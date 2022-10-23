N_ERR_COMPUTING_IN_DEVICE?=100
FLOAT_TYPE?=double
COMMON_CC_FLAGS=-DFLOAT_TYPE=$(FLOAT_TYPE) -fopenmp -lm -O3
COMMON_CPP_FLAGS=-std=c++17
OPENACC_CC_FLAGS=$(COMMON_CC_FLAGS) -acc -Minfo=accel
GPU_CC_FLAGS=$(OPENACC_CC_FLAGS) -ta=tesla -DTARGET_DEVICE=GPU -Mcudalib=cublas
CPU_CC_FLAGS=$(OPENACC_CC_FLAGS) -ta=multicore -DTARGET_DEVICE=CPU -Mcudalib=cublas

all: build

rebuild: clean build

build: no_blas blas_naive blas_without_err_copy blas_device_pointer_mode
 
no_blas: src/heat_equation_solver_no_blas.c heat_equation_runner.o src/heat_equation_solver.h src/heat_equation_utils.h
	pgcc $(GPU_CC_FLAGS) -std=c11 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_gpu.out
	pgcc $(CPU_CC_FLAGS) -std=c11 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_cpu.out

blas_naive: src/heat_equation_solver_blas_naive.cu heat_equation_runner.o src/heat_equation_solver.h src/heat_equation_utils.h
	pgc++ $(GPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_gpu.out
	pgc++ $(CPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_cpu.out

blas_without_err_copy: src/heat_equation_solver_blas_without_err_copy.cu heat_equation_runner.o src/heat_equation_solver.h src/heat_equation_utils.h
	pgc++ $(GPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_gpu.out
	pgc++ $(CPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_cpu.out

blas_device_pointer_mode: src/heat_equation_solver_blas_device_pointer_mode.cu heat_equation_runner.o src/heat_equation_solver.h src/heat_equation_utils.h
	pgc++ $(GPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_gpu.out
	pgc++ $(CPU_CC_FLAGS) -std=c++17 -DN_ERR_COMPUTING_IN_DEVICE=$(N_ERR_COMPUTING_IN_DEVICE) $< heat_equation_runner.o -o $@_cpu.out

heat_equation_runner.o: src/heat_equation_runner.cpp
	g++ $(COMMON_CC_FLAGS) -c $< -o heat_equation_runner.o

clean:
	rm -rf *.o *.out *.so

.PHONY: run clean