#!/bin/bash

ETOL=${1:-1e-6}
GRID_SIZE=${2:-512}
MAX_ITER=${3:-1000000}

BENCHMARKS_TABLE='benchmarks.csv'

if [ ! -f $BENCHMARKS_TABLE ]; then
    echo "Target device;Algo ver;Grid size;Number of iters;Elapsed Time;Last error" > $BENCHMARKS_TABLE
fi

VERSIONS="no_blas blas_naive blas_without_err_copy blas_device_pointer_mode"

make clean
make

for version in $VERSIONS; do
    ./${version}_cpu.out $ETOL $GRID_SIZE $MAX_ITER >> $BENCHMARKS_TABLE
    ./${version}_gpu.out $ETOL $GRID_SIZE $MAX_ITER >> $BENCHMARKS_TABLE
done
