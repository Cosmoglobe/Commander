#!/bin/bash
#================================================================================
# Description: This script will install Commander3 on owls
#================================================================================
# Absolute path to Commander3 root directory
#comm3_root="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
comm3_root_dir="$(pwd)"
# Using regex to figure out which owl I am on.
owl_prefix="owl"
owl_suffix="\.uio\.no"
owl_1724="$owl_prefix+(1[7-9]|2[0-4])+$owl_suffix"
owl_2528="$owl_prefix+(2[5-8])+$owl_suffix"
owl_2930="$owl_prefix+(29|30)+$owl_suffix"
owl_3135="$owl_prefix+(3[1-5])+$owl_suffix"
owl_3637="$owl_prefix+(3[6-7])+$owl_suffix"
# Will compile commander only if on owl!
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
	module load gnu git/2.30.1 Intel_parallel_studio/2020/4.912 #Intel_parallel_studio/2018/3.051
	echo "(Re)loaded the following modules:"
	echo "$(module list)"
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
	# If directory is empty or doesn't have correct CMake produced structure, 
	# we simply remove its contents and start anew
	if [[ -e "$abs_path_to_build/CMakeCache.txt" && -d "$abs_path_to_build/CMakeFiles" 
		&& -e "$abs_path_to_build/Makefile" && -e "$abs_path_to_build/cmake_install.cmake" ]];
	then
		echo "Rebuilding the projects."
		cmake3 --build $comm3_root_dir/$build_dir --target install -j $physicalCpuCount 
	else
		echo "Building the project from scratch."
		rm -rf $abs_path_to_build/*
		# Executing CMake commands for the first time

		#using release build type
		cmake3 -DCMAKE_INSTALL_PREFIX:PATH="$comm3_root_dir/$build_dir/install" -DCMAKE_DOWNLOAD_DIRECTORY:PATH="$comm3_root_dir/downloads" -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort -DCFITSIO_USE_CURL:BOOL=OFF -DUSE_SYSTEM_FFTW:BOOL=OFF -DUSE_SYSTEM_CFITSIO:BOOL=OFF -DUSE_SYSTEM_HDF5:BOOL=OFF -DUSE_SYSTEM_HEALPIX:BOOL=OFF -S $comm3_root_dir -B $comm3_root_dir/$build_dir 

		#using release build type with som debugging flags
		#cmake3 -DCMAKE_INSTALL_PREFIX:PATH="$comm3_root_dir/$build_dir/install" -DCMAKE_DOWNLOAD_DIRECTORY:PATH="$comm3_root_dir/downloads" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DMPI_C_COMPILER=mpiicc -DMPI_CXX_COMPILER=mpiicpc -DMPI_Fortran_COMPILER=mpiifort -DCFITSIO_USE_CURL:BOOL=OFF -DUSE_SYSTEM_FFTW:BOOL=OFF -DUSE_SYSTEM_CFITSIO:BOOL=OFF -DUSE_SYSTEM_HDF5:BOOL=OFF -DUSE_SYSTEM_HEALPIX:BOOL=OFF -S $comm3_root_dir -B $comm3_root_dir/$build_dir 
		# Build and install command
		cmake3 --build $comm3_root_dir/$build_dir --target install -j $physicalCpuCount 

	fi


fi
