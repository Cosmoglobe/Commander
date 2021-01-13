#==============================================================================
# This file contains general instructions how to
# fetch and build the Commander dependencies
# Author: Maksym Brilenkov
#==============================================================================

# TODO:
# [x] Split this file into several different files, containing corresponding instructions (versions, variables, toolchains etc.);
# [ ] Change URL_MD5 to URL_HASH of every project;
# [x] Change compiler variables from list to string (but leave APPEND);
# [ ] Remove include_directory() and use target_include_directory() instead (for commander3 target);
# [ ] Add one variable which will force all libraries to be recompiled;
# [x] Add ExternalProject_Add_Step() to copy files instead of COMMAND (and make it ALWAYS OFF);

#------------------------------------------------------------------------------
# including compiler definitions
include(project_toolchains)
# including project defined variables
include(project_variables)
#------------------------------------------------------------------------------

#==============================================================================
# MAIN PROJECT DEPENDENCIES
#==============================================================================
message(STATUS "---------------------------------------------------------------")
message(STATUS "Looking for packages...")
# use this to write your own find_package
find_package(PkgConfig)
# We will be using Git to download some dependencies, so we need to check if git available
find_package(Git REQUIRED)
# finding math library
message(STATUS "math (m) libraries are: ${LIBM_LIBRARIES}")
# printing out the dl libs, which are also required on some unix systems
message(STATUS "dl libs are: ${CMAKE_DL_LIBS}")

unset(projects)
# project names
list(APPEND projects 
	tempita
	blas # blas-lapack module 
	mpi
	openmp
	curl
	zlib
	##sharp2
	fftw
	cfitsio
	hdf5
	doxygen
	healpix
	commander3
	)
#==============================================================================
# PROJECTS' URL SOURCES, MD5 HASHES AND CONFIGURE COMMANDS
#==============================================================================
# cURL - needed by CFitsio and HEALPix
# need to specify command separately, othewise it won't work
set(curl_url "https://github.com/curl/curl/releases/download/curl-7_69_0/curl-7.69.0.zip")#"https://github.com/curl/curl/releases/download/curl-7_69_1/curl-7.69.1.tar.gz")
#------------------------------------------------------------------------------
# FFTW
set(fftw_url "http://fftw.org/fftw-3.3.8.tar.gz")
set(fftw_md5 "8aac833c943d8e90d51b697b27d4384d")
#------------------------------------------------------------------------------
# HDF5
set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz")
set(hdf5_md5 "e115eeb66e944fa7814482415dd21cc4")
#------------------------------------------------------------------------------
# LibSharp2
# Note: libsharp2 comes as a native part of Healpix 3.70 and thus, we do not
# need to compile it independently. But, I will leave this just in case we need
# to switch to older versions (for whatever reason).
set(sharp2_url "https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz")#"https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz") #"https://github.com/Libsharp/libsharp/archive/v1.0.0.tar.gz")
#------------------------------------------------------------------------------
# CFitsio
set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz")
#------------------------------------------------------------------------------
# HEALPix
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.tar.gz/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.zip/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")#"https://sourceforge.net/projects/healpix/files/latest/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")
#set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.tar.gz/download")
set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.70/Healpix_3.70_2020Jul23.tar.gz/download")
#set(healpix_md5 "ed7c9a3d7593577628ed1286fa7a9250")
#set(healpix_md5 "540b243406596205a7a82434d99af41e")
#set(healpix_md5 "9b51b2fc919f4e70076d296826eebee0")
set(healpix_md5 "bdcc2a4b1ede3ed5a07be57e4aec01d2")
# this command is for healpix 3.50 and below
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#set(healpix_configure_command "${CMAKE_COMMAND}" "-E" "env" "FC=${COMMANDER3_Fortran_COMPILER}" "CXX=${COMMANDER3_CXX_COMPILER}" "CPP=${COMMANDER3_CPP_COMPILER}" "CC=${COMMANDER3_C_COMPILER}" "./configure")
#------------------------------------------------------------------------------
# Doxygen
# there is some weird errors appearing for doxygen v1.8.17 and above, so will stick with this one
set(doxygen_url "https://github.com/doxygen/doxygen/archive/Release_1_8_16.tar.gz")#"https://github.com/doxygen/doxygen/archive/Release_1_8_18.tar.gz")#"https://github.com/doxygen/doxygen.git")
# compile doxygen with cmake itself, so no need to specify configure options here
set(flex_url "http://sourceforge.net/projects/flex/files/flex-2.5.39.tar.gz/download")#"https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.xz/download")
set(bison_url "http://ftp.gnu.org/gnu/bison/bison-3.6.tar.gz")#"http://ftp.gnu.org/gnu/bison/bison-3.6.2.tar.gz")
#------------------------------------------------------------------------------

# include all project configuration files
foreach(project ${projects})
	include("${project}")
endforeach()
# 
include_directories(${CMAKE_INSTALL_PREFIX}/include)
message(STATUS "---------------------------------------------------------------")
