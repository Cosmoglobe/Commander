#!/bin/bash
#================================================================================
# Description: This script will install Commander3 on NOTUR system: Betzy.
#================================================================================
comm3_root_dir="$(pwd)"
build_dir="build_notur"
buildtype="RelWithDebInfo" #"Debug" #"Release" #"RelWithDebInfo"

#if [[ "${HOSTNAME}" =~ *"saga"* ]]
#then
module purge
module load intel/2020b 
module load FFTW/3.3.8-intel-2020b
module load HDF5/1.10.7-iimpi-2020b
#module load CMake/3.20.1-GCCcore-10.3.0 
module load git/2.28.0-GCCcore-10.2.0-nodocs
# Printing Loaded modules
printf "\n"
printf "$(module list)"

# Compilers
fc="ifort"
cc="icc"
cxx="icpc"
# MPI compilers
mpifc="mpiifort" 
mpicc="mpiicc"
mpicxx="mpiicpc"

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
-DUSE_SYSTEM_FFTW:BOOL=ON \
-DUSE_SYSTEM_CFITSIO:BOOL=OFF \
-DUSE_SYSTEM_HDF5:BOOL=ON \
-DUSE_SYSTEM_HEALPIX:BOOL=OFF \
-DUSE_SYSTEM_BLAS:BOOL=ON \
-S $comm3_root_dir -B $abs_path_to_build
#------------------------------------------------------------------------------
# Build and install command
#------------------------------------------------------------------------------
cmake --build $abs_path_to_build --target install -j 2
#fi
