#!/bin/bash
#================================================================================
# Description: This script will install Commander3 on ITA systems: Owls/Beehives.
#================================================================================
# Global configuration:
#------------------------------------------------------------------------------
# Compiler Toolchain to use
# Possible values: nvidia, flang, gnu, intel
toolchain="gnu"
toolchain="intel"
buildtype="Release" #"Debug" #"Release" #"RelWithDebInfo"
#------------------------------------------------------------------------------
# Absolute path to Commander3 root directory
comm3_root_dir="$(pwd)"
# Address
suffix="\.uio\.no"
# Using regex to figure out which owl I am on.
prefix="owl"
owl1724="$prefix+(1[7-9]|2[0-4])+$suffix"
owl2528="$prefix+(2[5-8])+$suffix"
owl2930="$prefix+(29|30)+$suffix"
owl3135="$prefix+(3[1-5])+$suffix"
owl3637="$prefix+(3[6-7])+$suffix"
# Using regex to figure out which beehive I am on.
prefix="beehive"
bee0123="$prefix+(\d{0}|[1-9](?!\d)|1[0-9]|2[0-3])+$suffix"
bee2631="$prefix+(2[6-9]|3[0-1])+$suffix"
bee3436="$prefix+(3[4-6])+$suffix"
bee3742="$prefix+(3[7-9]|4[0-2])+$suffix"
bee4345="$prefix+(4[3-5])+$suffix"
bee46="$prefix+(46)+$suffix"
bee47="$prefix+(47)+$suffix"
#------------------------------------------------------------------------------
# Will compile commander only if on owl/beehive!
#------------------------------------------------------------------------------
if [[ "${HOSTNAME}" =~ "owl"* ]] || [[ "${HOSTNAME}" =~ "beehive"* ]]
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
  if [[ "${HOSTNAME}" =~ $owl1724 ]]; then
    build_dir="build_owl1724_$toolchain"
  elif [[ "${HOSTNAME}" =~ $owl2528 ]]; then
    build_dir="build_owl2528_$toolchain"
  elif [[ "${HOSTNAME}" =~ $owl2930 ]]; then
    build_dir="build_owl2930_$toolchain"
  elif [[ "${HOSTNAME}" =~ $owl3135 ]]; then
    build_dir="build_owl3135_$toolchain"
  elif [[ "${HOSTNAME}" =~ $owl3637 ]]; then
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
		module load Intel_parallel_studio/2020/4.912
		#module load Intel_parallel_studio/2018/3.051
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
		#source /opt/rh/devtoolset-9/enable
		#export PATH="/usr/local/opt/openmpi-4.0.5/bin:$PATH"
		#export LD_LIBRARY_PATH="/usr/local/opt/openmpi-4.0.5/lib:$LD_LIBRARY_PATH"
		module load gcc/10.2.1
		module load myopenmpi/4.0.3
		#module load gcc/9.3.1 Mellanox/2.8.1/gcc/hpcx
		printf "\n"
		$mpifc --version
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
		module load openmpi/aocc/4.1.0 AMD/aocc/3.0.0
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
	-DUSE_SYSTEM_BLAS:BOOL=ON \
	-S $comm3_root_dir -B $abs_path_to_build
	#------------------------------------------------------------------------------
	# Build and install command
	#------------------------------------------------------------------------------
	cmake --build $comm3_root_dir/$build_dir --target install -j $physicalCpuCount #-v 
else
	printf "TERMINATING: NOT ON ITA MACHINE!"
fi
