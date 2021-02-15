#!/usr/bin/bash

HEALPIX_ROOT=/path/to/healpix/root/healpix
COMM3_BIN=/path/to/commander/binary/commander3
PARAM_FILE=proc_BP8.1.txt
CHAINS_DIR=proc_chains
OUT_FILE=slurm.txt
NUM_PROC=64

export OMP_NUM_THREADS=1
export HEALPIX=$HEALPIX_ROOT
mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
