#!/usr/bin/bash
#================================================================================
# Global Configuration
#------------------------------------------------------------------------------
comm3_root_dir="$HOME/mycomm/maksymb"
PARAM_FILE=param_maksymb.txt 
CHAINS_DIR=chains_maksymb 
OUT_FILE=slurm.txt
NUM_PROC=8 #32 #48
COMM3_BIN=commander3
# Valid options are: gnu, intel, nvhpc
toolchain="intel"
#------------------------------------------------------------------------------
owl_prefix="owl"
owl_suffix="\.uio\.no"
owl_1724="$owl_prefix+(1[7-9]|2[0-4])+$owl_suffix"
owl_2528="$owl_prefix+(2[5-8])+$owl_suffix"
owl_2930="$owl_prefix+(29|30)+$owl_suffix"
owl_3135="$owl_prefix+(3[1-5])+$owl_suffix"
owl_3637="$owl_prefix+(3[6-7])+$owl_suffix"
# Using regex to figure out which beehive I am on.
prefix="beehive"
bee0123="$prefix+(\d{0}|[1-9](?!\d)|1[0-9]|2[0-3])+$suffix"
bee2631="$prefix+(2[6-9]|3[0-1])+$suffix"
bee3436="$prefix+(3[4-6])+$suffix"
bee3742="$prefix+(3[7-9]|4[0-2])+$suffix"
bee4345="$prefix+(4[3-5])+$suffix"
bee46="$prefix+(46)+$suffix"
bee47="$prefix+(47)+$suffix"
if [[ "${HOSTNAME}" =~ "owl"* ]] || [[ "${HOSTNAME}" =~ "beehive"* ]]
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
		build_dir="build_owl1724_$toolchain"
	elif [[ "${HOSTNAME}" =~ $owl_2528 ]]; then
		build_dir="build_owl2528_$toolchain"
	elif [[ "${HOSTNAME}" =~ $owl_2930 ]]; then
		build_dir="build_owl2930_$toolchain"
	elif [[ "${HOSTNAME}" =~ $owl_3135 ]]; then
		build_dir="build_owl3135_$toolchain"
	elif [[ "${HOSTNAME}" =~ $owl_3637 ]]; then
		build_dir="build_owl3637_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee0123 ]]; then
    build_dir="build_bee0123_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee2631 ]]; then
    build_dir="build_bee2631_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee3436 ]]; then
    build_dir="build_bee3436_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee3742 ]]; then
    build_dir="build_bee3742_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee4345 ]]; then
    build_dir="build_bee4345_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee46 ]]; then
    build_dir="build_bee46_$toolchain"
  elif [[ "${HOSTNAME}" =~ $bee47 ]]; then
    build_dir="build_bee47_$toolchain"
	fi
	echo "$build_dir"
	#------------------------------------------------------------------------------
	# (Re)loading necessary modules
	module purge
	module load gnu git/2.30.1 cmake/3.21.1 
  if [[ "$toolchain" =~ "intel" ]]
  then
		#module load Intel_parallel_studio/2018/3.051
		module load Intel_parallel_studio/2020/4.912
  elif [[ "$toolchain" =~ "gnu" ]]
  then
		source /opt/rh/devtoolset-9/enable
		module load myopenmpi/4.0.3
  elif [[ "$toolchain" =~ "flang" ]]
  then
		module load openmpi/aocc/4.0.5 AMD/aocc/3.0.0
  elif [[ "$toolchain" =~ "nvhpc" ]]
  then
		module load nvhpc/21.7
	fi
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
	#mpirun -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
	#mpirun -np $physicalCpuCount $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
	mpirun -env I_MPI_FABRICS=shm -env I_MPI_PIN_DOMAIN=numa -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
	#mpirun -env I_MPI_FABRICS=shm -env I_MPI_PIN_DOMAIN=numa -np $physicalCpuCount $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE
fi
#------------------------------------------------------------------------------
