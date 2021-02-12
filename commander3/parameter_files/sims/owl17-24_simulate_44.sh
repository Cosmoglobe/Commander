#!/usr/bin/bash

HEALPIX_ROOT=/mn/stornext/u3/maksymb/commander/bp8.1_sims/build_owl17-24r/install/healpix
COMM3_BIN=/mn/stornext/u3/maksymb/commander/bp8.1_sims/build_owl17-24r/install/bin/commander3
PARAM_FILE=sims_BP8.1_44.txt
CHAINS_DIR=chains_44
OUT_FILE=slurm.txt
NUM_PROC=24

export OMP_NUM_THREADS=1
export HEALPIX=$HEALPIX_ROOT
mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
