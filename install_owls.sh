#!/bin/bash
#================================================================================
# Description: This script will install Commander3 on owls
#================================================================================
# Global configuration:
#------------------------------------------------------------------------------
# Compiler Toolchain to use
# Possible values: nvidia, flang, gnu, intel
toolchain="intel"
buildtype="RelWithDebInfo" #"Debug" #"Release" #"RelWithDebInfo"
#------------------------------------------------------------------------------
# Absolute path to Commander3 root directory
comm3_root_dir="$(pwd)"
# Using regex to figure out which owl I am on.
owl_prefix="owl"
owl_suffix="\.uio\.no"
owl_1724="$owl_prefix+(1[7-9]|2[0-4])+$owl_suffix"
owl_2528="$owl_prefix+(2[5-8])+$owl_suffix"
owl_2930="$owl_prefix+(29|30)+$owl_suffix"
owl_3135="$owl_prefix+(3[1-5])+$owl_suffix"
owl_3637="$owl_prefix+(3[6-7])+$owl_suffix"
#------------------------------------------------------------------------------
# Will compile commander only if on owl!
#------------------------------------------------------------------------------
if [[ "${HOSTNAME}" =~ "owl"* ]]
then
	#------------------------------------------------------------------------------
	# Getting the total number of CPUs, taken from this answer:
	# https://stackoverflow.com/questions/6481005/how-to-obtain-the-number-of-cpus-cores-in-linux-from-the-command-line
	# macOS:           Use `sysctl -n hw.*cpu_max`, which returns the values of
	#                  interest directly.
	#                  CAVEAT: Using the "_max" key suffixes means that the *maximum*
	#                          available number of CPUs is reported, whereas the
	#                          current power-management mode could make *fewer* CPUs
	#                          available; dropping the "_max" suffix would report the
	#                          number of *currently* available ones; see [1] below.
	#
	# Linux:           Parse output from `lscpu -p`, where each output line represents
	#                  a distinct (logical) CPU.
	#                  Note: Newer versions of `lscpu` support more flexible output
	#                        formats, but we stick with the parseable legacy format
	#                        generated by `-p` to support older distros, too.
	#                        `-p` reports *online* CPUs only - i.e., on hot-pluggable
	#                        systems, currently disabled (offline) CPUs are NOT
	#                        reported.
	#
	# Number of LOGICAL CPUs (includes those reported by hyper-threading cores)
	# Linux: Simply count the number of (non-comment) output lines from `lscpu -p`,
	# which tells us the number of *logical* CPUs.
	logicalCpuCount=$([ $(uname) = 'Darwin' ] &&
												 sysctl -n hw.logicalcpu_max ||
												 lscpu -p | egrep -v '^#' | wc -l)

	# Number of PHYSICAL CPUs (cores).
	# Linux: The 2nd column contains the core ID, with each core ID having 1 or
	#        - in the case of hyperthreading - more logical CPUs.
	#        Counting the *unique* cores across lines tells us the
	#        number of *physical* CPUs (cores).
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
	# Unloading any loaded module
	module purge
	# Loading GNU Autotools (autoconf, libtool, automake etc.), GIT and CMake
	module load gnu git/2.30.1 cmake/3.21.1
	# Choosing which compiler toolchain to use
	if [[ "$toolchain" =~ "intel" ]]
	then
		# Compilers
		fc="ifort"
		cc="icc"
		cxx="icpc"
		# MPI compilers
		mpifc="mpiifort" 
		mpicc="mpiicc"
		mpicxx="mpiicpc"
		printf "Using Intel:\nFC=$fc\nCC=$cc\nCXX=$cxx\nMPIF90=$mpifc\nMPICC=$mpicc\nMPICXX=$mpicxx"
		#module load Intel_parallel_studio/2020/4.912
		module load Intel_parallel_studio/2018/3.051
	elif [[ "$toolchain" =~ "gnu" ]]
	then
		# Compilers
		fc="gfortran"
		cc="gcc"
		cxx="g++"
		# MPI compilers
		mpifc="mpifort" 
		mpicc="mpicc"
		mpicxx="mpicxx"
		printf "Using GNU:\nFC=$fc\nCC=$cc\nCXX=$cxx\nMPIF90=$mpifc\nMPICC=$mpicc\nMPICXX=$mpicxx"
	  #module load foss/10.3.0 # custom GNU GCC + OpenMPI 
		#module load gcc/9.3.1 Mellanox/2.8.1/gcc/hpcx
		source /opt/rh/devtoolset-9/enable
		#export PATH="/usr/local/opt/openmpi-4.0.5/bin:$PATH"
		#export LD_LIBRARY_PATH="/usr/local/opt/openmpi-4.0.5/lib:$LD_LIBRARY_PATH"
		module load myopenmpi/4.0.3
	elif [[ "$toolchain" =~ "flang" ]]
	then
		# Compilers
		fc="flang"
		cc="clang"
		cxx="clang++"
		# MPI compilers
		mpifc="mpifort" 
		mpicc="mpicc"
		mpicxx="mpicxx"
		printf "Using AOCC:\nFC=$fc\nCC=$cc\nCXX=$cxx\nMPIF90=$mpifc\nMPICC=$mpicc\nMPICXX=$mpicxx"
		module load openmpi/aocc/4.0.5 AMD/aocc/3.0.0
	elif [[ "$toolchain" =~ "nvidia" ]]
	then
		# Compilers
		fc="nvfortran"
		cc="nvc"
		cxx="nvc++"
		# MPI compilers
		mpifc="mpifort" 
		mpicc="mpicc"
		mpicxx="mpicxx"
		printf "Using NVIDIA:\nFC=$fc\nCC=$cc\nCXX=$cxx\nMPIF90=$mpifc\nMPICC=$mpicc\nMPICXX=$mpicxx"
		module load nvhpc/21.7 
	fi
	# Printing Loaded modules
	printf "\n"
	printf "$(module list)"
	#------------------------------------------------------------------------------
	# Checking for existence of build directory
	echo "Checking for 'build' directory."
	abs_path_to_build="$comm3_root_dir/$build_dir"
	if [[ -d "$abs_path_to_build" ]]; 
	then
		echo "$abs_path_to_build exists! Proceeding..."
	else
		echo "$abs_path_to_build does not exist! Creating..."
		mkdir $abs_path_to_build 
	fi
	#------------------------------------------------------------------------------
	rm -rf $abs_path_to_build/CMakeCache.txt
	#------------------------------------------------------------------------------
	# Executing CMake commands for the first time
	#------------------------------------------------------------------------------
	cmake \
	-DCMAKE_INSTALL_PREFIX:PATH="$comm3_root_dir/$build_dir/install" \
	-DCMAKE_DOWNLOAD_DIRECTORY:PATH="$comm3_root_dir/downloads" \
	-DCMAKE_BUILD_TYPE="$buildtype" \
	-DCMAKE_Fortran_COMPILER=$fc \
	-DCMAKE_C_COMPILER=$cc \
	-DCMAKE_CXX_COMPILER=$cxx \
	-DMPI_C_COMPILER=$mpicc \
	-DMPI_CXX_COMPILER=$mpicxx \
	-DMPI_Fortran_COMPILER=$mpifc \
	-DCFITSIO_USE_CURL:BOOL=OFF \
	-DUSE_SYSTEM_FFTW:BOOL=OFF \
	-DUSE_SYSTEM_CFITSIO:BOOL=OFF \
	-DUSE_SYSTEM_HDF5:BOOL=ON \
	-DUSE_SYSTEM_HEALPIX:BOOL=OFF \
	-S $comm3_root_dir -B $abs_path_to_build
	#------------------------------------------------------------------------------
	# Build and install command
	#------------------------------------------------------------------------------
	cmake --build $comm3_root_dir/$build_dir --target install -j $physicalCpuCount #-v 
else
	printf "TERMINATING: NOT ON OWL!"
fi
