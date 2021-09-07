#!/usr/bin/bash
#================================================================================
# Global Configuration
#------------------------------------------------------------------------------
comm3_root_dir="$HOME/mycomm/stable"
PARAM_FILE=param_stable.txt #param_BP10_c48_test.txt
CHAINS_DIR=chains_stable #chains_BP10_c48_test
OUT_FILE=slurm.txt
NUM_PROC=48
COMM3_BIN=commander3
#------------------------------------------------------------------------------
owl_prefix="owl"
owl_suffix="\.uio\.no"
owl_1724="$owl_prefix+(1[7-9]|2[0-4])+$owl_suffix"
owl_2528="$owl_prefix+(2[5-8])+$owl_suffix"
owl_2930="$owl_prefix+(29|30)+$owl_suffix"
owl_3135="$owl_prefix+(3[1-5])+$owl_suffix"
owl_3637="$owl_prefix+(3[6-7])+$owl_suffix"
if [[ "${HOSTNAME}" =~ "owl"* ]]
then
	#------------------------------------------------------------------------------
	logicalCpuCount=$([ $(uname) = 'Darwin' ] &&
												sysctl -n hw.logicalcpu_max ||
												lscpu -p | egrep -v '^#' | wc -l)
  physicalCpuCount=$([ $(uname) = 'Darwin' ] &&
                         sysctl -n hw.physicalcpu_max ||
                         lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)

	#------------------------------------------------------------------------------
	echo "You are on ${HOSTNAME}"	
	#------------------------------------------------------------------------------
	# Choosing correct build directory to put CMake files into
	if [[ "${HOSTNAME}" =~ $owl_1724 ]]; then
		build_dir="build_owl1724"
	elif [[ "${HOSTNAME}" =~ $owl_2528 ]]; then
		build_dir="build_owl2528"
	elif [[ "${HOSTNAME}" =~ $owl_2930 ]]; then
		build_dir="build_owl2930"
	elif [[ "${HOSTNAME}" =~ $owl_3135 ]]; then
		build_dir="build_owl3135"
	elif [[ "${HOSTNAME}" =~ $owl_3637 ]]; then
		build_dir="build_owl3637"
	fi
	#------------------------------------------------------------------------------
	# (Re)loading necessary modules
	module purge
	module load gnu git/2.30.1 cmake/3.21.1 
	module load Intel_parallel_studio/2018/3.051
	# Ensuring there is no OpenMP involved
	export OMP_NUM_THREADS=1
	# Exporint the default part of paraneter file
	export COMMANDER_PARAMS_DEFAULT="$comm3_root_dir/commander3/parameter_files/defaults"
	# Adding compiled libs to correct paths
	export HEALPIX="$comm3_root_dir/$build_dir/install/healpix"
	export PATH="$comm3_root_dir/$build_dir/install/bin:$PATH"
	export LD_LIBRARY_PATH="$HEALPIX/lib:$comm3_root_dir/$build_dir/install/lib:$LD_LIBRARY_PATH"
	# Printing out the stack size and other info
	ulimit -a
	# Increasing the stack-size
	ulimit -s unlimited
	#------------------------------------------------------------------------------
	# Running Commander
	mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
fi
#------------------------------------------------------------------------------
#export PATH="$PATH"
#export CFITSIO="$HOME/mycomm/stable/build_owl3135/install/lib"
#export LD_LIBRARY_PATH="$CFITSIO:$HOME/mycomm/stable/build_owl3135/install/healpix/lib:$LD_LIBRARY_PATH"
#------------------------------------------------------------------------------
#HEALPIX_ROOT=/mn/stornext/u3/maksymb/mycomm/stable/build_owl3135/install/healpix
#COMM3_BIN=/mn/stornext/u3/maksymb/mycomm/stable/build_owl3135/install/bin/commander3
#------------------------------------------------------------------------------
#source /opt/rh/devtoolset-9/enable
#module load myopenmpi/4.0.3
#module load Intel_parallel_studio/2020/4.912 
#export HEALPIX=$HEALPIX_ROOT
