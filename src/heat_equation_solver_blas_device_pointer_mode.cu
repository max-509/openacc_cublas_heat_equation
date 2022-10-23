#include "heat_equation_solver_impl.h"

#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#undef fabs
#undef fmax

#include <cublas_v2.h>

#include <type_traits>

#pragma acc routine seq
template <typename T>
constexpr T generic_abs(T v)
{
  if constexpr (std::is_same_v<T, float>)
  {
    return fabsf(v);
  }
  else
  {
    return fabs(v);
  }
}

#if (TARGET_DEVICE == GPU)

#define CUBLAS_CHECK(err)                                                           \
  do                                                                                \
  {                                                                                 \
    cublasStatus_t err_ = (err);                                                    \
    if (err_ != CUBLAS_STATUS_SUCCESS)                                              \
    {                                                                               \
      std::fprintf(stderr, "cublas error %d at %s:%d\n", err_, __FILE__, __LINE__); \
      throw std::runtime_error("cublas error");                                     \
    }                                                                               \
  } while (0)

template <typename T>
cublasStatus_t cublasCopy(cublasHandle_t handle,
                          int n,
                          const T *x,
                          int incx,
                          T *y,
                          int incy)
{
  if constexpr (std::is_same_v<float, T>)
  {
    return cublasScopy(handle, n, x, incx, y, incy);
  }
  else
  {
    return cublasDcopy(handle, n, x, incx, y, incy);
  }
}

template <typename T>
cublasStatus_t cublasAxpy(cublasHandle_t handle,
                          int n,
                          const T *alpha, /* host or device pointer */
                          const T *x,
                          int incx,
                          T *y,
                          int incy)
{
  if constexpr (std::is_same_v<float, T>)
  {
    return cublasSaxpy(handle, n, alpha, x, incx, y, incy);
  }
  else
  {
    return cublasDaxpy(handle, n, alpha, x, incx, y, incy);
  }
}

template <typename T>
cublasStatus_t cublasIamax(cublasHandle_t handle,
                           int n,
                           const T *x,
                           int incx,
                           int *result)
{
  if constexpr (std::is_same_v<float, T>)
  {
    return cublasIsamax(handle, n, x, incx, result);
  }
  else
  {
    return cublasIdamax(handle, n, x, incx, result);
  }
}

template <typename T>
struct err_computer
{
  T operator()(const T *__restrict__ buff_grid, T *__restrict__ diff_buff,
               const size_t grid_size, cublasHandle_t handle)
  {
    const size_t grid_sqr = grid_size * grid_size;

    T err;
    T a;
    int err_idx;

#pragma acc declare create(err) device_resident(a, err_idx)
    {

#pragma acc data present(a)
#pragma acc kernels num_gangs(1) num_workers(1)
      a = static_cast<T>(-1.0);

#pragma acc host_data use_device(buff_grid, diff_buff)
      {
        CUBLAS_CHECK(cublasCopy(handle, grid_sqr, buff_grid, 1, diff_buff, 1));
#pragma acc host_data use_device(a)
        CUBLAS_CHECK(cublasAxpy(handle, grid_sqr, &a, buff_grid + grid_sqr, 1, diff_buff, 1));
      }

#pragma acc host_data use_device(err_idx, diff_buff)
      CUBLAS_CHECK(cublasIamax(handle, grid_sqr, diff_buff, 1, &err_idx));
#pragma acc data present(err, err_idx, diff_buff)
#pragma acc kernels num_gangs(1) num_workers(1)
      err = diff_buff[err_idx];
#pragma acc update host(err)
    }

    return generic_abs(err);
  }
};

#else

template <typename T>
struct err_computer
{
  T operator()(T *__restrict__ buff_grid, T *__restrict__ diff_buff /*not used*/,
               const size_t grid_size, cublasHandle_t handle)
  {
    const size_t grid_sqr = grid_size * grid_size;
    T err = 0.0;
#pragma acc wait

#pragma acc data present(buff_grid [0:grid_sqr * 2])
    {
#pragma acc kernels
      {
#pragma acc loop independent collapse(2) reduction(max \
                                                   : err)
        for (size_t i = 1; i < grid_size - 1; ++i)
        {
          for (size_t j = 1; j < grid_size - 1; ++j)
          {
            const size_t grid_idx = i * grid_size + j;
            err = max(err, abs(buff_grid[grid_idx] - buff_grid[grid_sqr + grid_idx]));
          }
        }
      }
    }

    return err;
  }
};

#endif // TARGET_DEVICE

#ifndef N_ERR_COMPUTING_IN_DEVICE
#define N_ERR_COMPUTING_IN_DEVICE 1500
#endif // N_ERR_COMPUTING_IN_DEVICE

int solve_heat_equation(FLOAT_TYPE *__restrict__ init_grid, const size_t grid_size, const size_t max_iter, const FLOAT_TYPE etol, size_t *last_iter, FLOAT_TYPE *last_etol)
{
  const size_t grid_sqr = grid_size * grid_size;
  FLOAT_TYPE *__restrict__ buff_grid = (FLOAT_TYPE *)malloc(sizeof(FLOAT_TYPE) * (grid_sqr * 2));
  if (NULL == buff_grid)
  {
    return 1;
  }

  cublasHandle_t cublas_handle = NULL;
  if (CUBLAS_STATUS_SUCCESS != cublasCreate(&cublas_handle))
  {
    free(buff_grid);
    return 1;
  }

  cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_DEVICE);

#pragma acc data copy(init_grid [0:grid_sqr]) create(buff_grid [0:grid_sqr * 2])
  {
    FLOAT_TYPE err = (FLOAT_TYPE)INFINITY;

    size_t curr_iter;

#pragma acc data present(init_grid [0:grid_sqr], buff_grid [0:grid_sqr * 2])
#pragma acc parallel
    {
#pragma acc loop independent collapse(2)
      for (size_t i = 0; i < grid_size; ++i)
      {
        for (size_t j = 0; j < grid_size; ++j)
        {
          size_t grid_idx = i * grid_size + j;
          buff_grid[grid_idx] = init_grid[grid_idx];
          buff_grid[grid_sqr + grid_idx] = init_grid[grid_idx];
        }
      }
    }

    size_t n_err_iter;
    for (curr_iter = 0u; curr_iter < max_iter && err > etol; curr_iter += n_err_iter)
    {

      for (n_err_iter = 0; n_err_iter < N_ERR_COMPUTING_IN_DEVICE; n_err_iter += 2)
      {
#pragma acc data present(buff_grid [0:grid_sqr * 2])
#pragma acc kernels async
        {
#pragma acc loop independent collapse(2)
          for (size_t i = 1; i < grid_size - 1; ++i)
          {
            for (size_t j = 1; j < grid_size - 1; ++j)
            {
              const size_t grid_idx = i * grid_size + j;
              buff_grid[grid_sqr + grid_idx] = (FLOAT_TYPE)0.25 * (buff_grid[grid_idx - grid_size] +
                                                                   buff_grid[grid_idx + grid_size] +
                                                                   buff_grid[grid_idx - 1] +
                                                                   buff_grid[grid_idx + 1]);
            }
          }

#pragma acc loop independent collapse(2)
          for (size_t i = 1; i < grid_size - 1; ++i)
          {
            for (size_t j = 1; j < grid_size - 1; ++j)
            {
              const size_t grid_idx = i * grid_size + j;
              const size_t next_grid_idx = grid_sqr + grid_idx;
              buff_grid[grid_idx] = (FLOAT_TYPE)0.25 * (buff_grid[next_grid_idx - grid_size] +
                                                        buff_grid[next_grid_idx + grid_size] +
                                                        buff_grid[next_grid_idx - 1] +
                                                        buff_grid[next_grid_idx + 1]);
            }
          }
        }
      }

      err = err_computer<FLOAT_TYPE>{}(buff_grid, init_grid, grid_size, cublas_handle);
    }

#pragma acc wait

#pragma acc data present(init_grid [0:grid_sqr], buff_grid [0:grid_sqr * 2])
#pragma acc parallel
    {
#pragma acc loop independent collapse(2)
      for (size_t i = 0; i < grid_size; ++i)
      {
        for (size_t j = 0; j < grid_size; ++j)
        {
          size_t grid_idx = i * grid_size + j;
          init_grid[grid_idx] = buff_grid[grid_idx];
        }
      }
    }

    if (NULL != last_iter)
    {
      *last_iter = curr_iter;
    }
    if (NULL != last_etol)
    {
      *last_etol = err;
    }
  }

  cublasDestroy(cublas_handle);
  free(buff_grid);

  return 0;
}

const char *get_solver_version()
{
  return "BLAS device pointer mode";
}
