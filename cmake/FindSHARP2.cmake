#==============================================================================
# Module to find libsharp2 on the system.
# Author: Maksym Brilenkov
#
# Libsharp2 is a library for efficient spherical harmonic transforms 
# at arbitrary spins, supporting CPU vectorization, OpenMP and MPI.
# To obtain libsharp2 visit main gitlab repo via: 
# https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master
#
# This file is part of Commander. For  details, see <LINK>. The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the distribution or at <LICENSE>.
#
# Usage:
#   find_package(SHARP2 [REQUIRED] [QUIET])
#
# The following variables are defined:
#   SHARP2_FOUND                 set to TRUE if SHARP2 is found on the system;
#   SHARP2_LIBRARIES             full path to SHARP2 libraries;
#   SHARP2_INCLUDE_DIR           SHARP2 include directory path; 
#==============================================================================
# Reference to CMake docs:
# https://cmake.org/cmake/help/v3.17/module/FindPackageHandleStandardArgs.html
include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)
#==============================================================================
# First check whether the SHARP2_ROOT variable was set by the user,
# as this variable can be used to guide detection of the libsharp2
# to a non-standard installation directory. The reference point will 
# be environment variable - SHARP2 
if(NOT SHARP2_ROOT)
	set(SHARP2_ROOT "$ENV{SHARP2}")
endif()
#==============================================================================
# Looking for headers and libraries. Libsharp2 doesn't have an enclude folder.
# Instead there is libsharp2 folder which contains all header files. Libsharp2
# also doesn't have an install command, but only compile one. It compiles all 
# the libraries into .libs folder. However, It is possible to compile test 
# programs. Keeping this in mind we do the following:

# First, we need to look for user installed (if any)
# version of libsharp2. The next will go Commander3
# compiled and copied version. And then we go to a 
# bunch of standard locations, where libsharp2, can be
# potentially stored.
find_path(SHARP2_INCLUDE_DIR 
	NAMES sharp.h
	PATH_SUFFIXES libsharp2
	HINTS ${SHARP2_ROOT} 
				${CMAKE_INSTALL_PREFIX}
				${CMAKE_INSTALL_PREFIX}/sharp2
				/usr
				/usr/src
				/usr/src/libsharp-master
				/usr/local
				/usr/local/src
				/usr/local/src/libsharp-master
				~/Downloads
				~/Downloads/libsharp-master
				~/Downloads/libsharp2
				~/Downloads/sharp2
				~/local
				~/local/src
				~/local/src/libsharp-master
	)

# This one is similar to the procedure for the 
# header file.
find_library(SHARP2_LIBRARIES
	NAMES "sharp2"
	PATH_SUFFIXES .libs lib libs lib64
	HINTS ${SHARP2_ROOT}
				${CMAKE_INSTALL_PREFIX}
				${CMAKE_INSTALL_PREFIX}/sharp2
				/usr
				/usr/src
				/usr/src/libsharp-master
				/usr/local
				/usr/local/src
				/usr/local/src/libsharp-master
				~/Downloads
				~/Downloads/libsharp-master
				~/Downloads/libsharp2
				~/Downloads/sharp2
				~/local
				~/local/src
				~/local/src/libsharp-master
	)
#==============================================================================
# Once libraries are found (or not found), we need to notify user about 
# the results and export a bunch of variables. 
find_package_handle_standard_args(SHARP2
	FOUND_VAR			SHARP2_FOUND 
	REQUIRED_VARS SHARP2_INCLUDE_DIR SHARP2_LIBRARIES
	)
