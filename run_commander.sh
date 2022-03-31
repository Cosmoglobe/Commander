#!/usr/bin/bash

BUILD_DIR=/mn/stornext/u3/artemba/devel_commander/Commander/build_owl2528_intel
COMM3_BIN=$BUILD_DIR/install/bin/commander3
PARAM_FILE=param_BP10_final_v2.txt
CHAINS_DIR=chains_directory
OUT_FILE=slurm.txt
NUM_PROC=64

export LD_LIBRARY_PATH=$BUILD_DIR/install/lib:$BUILD_DIR/install/healpix/lib:$LD_LIBRARY_PATH
export HEALPIX=$BUILD_DIR/install/healpix
export COMMANDER_PARAMS_DEFAULT=$HOME/devel_commander/Commander/commander3/parameter_files/defaults
export OMP_NUM_THREADS=1
mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
