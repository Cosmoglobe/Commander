#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================
# Author: Maksym Brilenkov
#================================================================================
# Module to find libsharp2 on the system.
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
