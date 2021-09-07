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
# This file contains general instructions how to
# fetch and build the Commander dependencies
#==============================================================================

#==============================================================================
# PROJECT'S DOWNLOAD/INSTALL/OUTPUT DIRECTORIES 
#==============================================================================
# Download dir
# Setting default values if this variable wasn't defined
set(CMAKE_DOWNLOAD_DIRECTORY "${CMAKE_SOURCE_DIR}/build/downloads"
	CACHE STRING
	"Directory where to download commander dependencies' source files"
	)
# Where to output shared libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib"
	CACHE STRING
	"Directory where to install all the libraries."
	)
# Where to output static libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib"
	CACHE STRING
	"Directory where to install all the libraries."
	)
# Where to output executable(s)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/bin"
	CACHE STRING
	"Directory where to install all the executables."
	)
# setting the directory where to output all .mod and .o files
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_INSTALL_PREFIX}/mod"
	CACHE STRING
	"Directory where to install all .mod/.o files."
	)
# HEALPix install (root) dir - by default we will install it in "healpix" dir
set(HEALPIX_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/healpix"
	CACHE STRING
	"Directory where to install (copy compiled and data files of) HEALPix."
	)
# Where to store logs of each step (download, configure, build, install)
set(CMAKE_LOG_DIR "${CMAKE_INSTALL_PREFIX}/logs"
	CACHE STRING
	"Derictory where to output logs of each step (download, configure, build, install)."
	)
# To avoid errors such as "file doesn't exists" during
# configuration/build/installation of an external project,
# we need to manually recreate this folder (during configure time).
if(NOT EXISTS "${CMAKE_LOG_DIR}")
  file(MAKE_DIRECTORY "${CMAKE_LOG_DIR}")
endif()
set(DOXYGEN_BUILD_DOCS OFF 
	CACHE BOOL
	"Determine whether to use doxygen or not."
	)
#------------------------------------------------------------------------------
# If any problems with installation will occur, which cannot be fixed quickly,
# these variables will force a fresh installation for every specified library.
# If USE_SYSTEM_LIBS is set to FALSE, all libraries will be recompiled,
# whereas if set to TRUE, the libraries will first be searched on the
# system and only if not found, compiled from source. However, If USE_SYSTEM_LIBS
# is set to TRUE but, e.g. USE_SYSTEM_HDF5 is set to FALSE, then HDF5 will
# be compiled from source (i.e. the latter variable is having the advantage).
#------------------------------------------------------------------------------
option(USE_SYSTEM_LIBS    "Enables search for LIBS on the system."        ON)
# BLAS/LAPACK
option(USE_SYSTEM_BLAS    "Enables search for BLAS/LAPACK on the system." ON)
# HDF5
option(USE_SYSTEM_ZLIB    "Enables search for ZLIB on the system."        ON)
option(USE_SYSTEM_LIBAEC  "Enables search for LibAEC on the system."      ON)
option(USE_SYSTEM_HDF5    "Enables search for HDF5 on the system."        OFF)
# FFTW
option(USE_SYSTEM_FFTW    "Enables search for HDF5 on the system."        OFF)
# CFITSIO
option(USE_SYSTEM_MBEDTLS "Enables search for MbedTLS on the system."     ON)
option(USE_SYSTEM_LIBSSH2 "Enables search for LibSSH2 on the system."     ON)
option(USE_SYSTEM_CURL    "Enables search for cURL on the system."        ON)
option(USE_SYSTEM_CFITSIO "Enables search for CFITSIO on the system."     OFF)
# Can choose whether to compile CFITSIO with or without cURL support
option(CFITSIO_USE_CURL   "Installs CFITSIO with cURL support."           OFF)
# HEALPix
option(USE_SYSTEM_HEALPIX "Enables search for HEALPIX on the system."     OFF)
# Doxygen
option(USE_SYSTEM_FLEX    "Enables search for FLEX on the system."        ON)
option(USE_SYSTEM_BISON   "Enables search for BISON on the system."       ON)
option(USE_SYSTEM_DOXYGEN "Enables search for DOXYGEN on the system."     ON)
#------------------------------------------------------------------------------
# Commander3 source dir
set(COMMANDER3_SOURCE_DIR "${CMAKE_SOURCE_DIR}/commander3/src")
# tempita source dir
set(TEMPITA_DIR ${CMAKE_SOURCE_DIR}/commander3/python)
# Some projects are build with configure and others are build with CMake. 
# To avoid errors when compiling on different owls (or other machines),
# we put special variable prefix for all subprojects' builds to be located 
# inside build directory.
set(LIBS_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/subbuilds")
# The following variables' purpose is to keep source directory separate from
# build directory etc. In such way we are able to safely remove the build dir
# without removing sources; thus, can easily recompile, no need to redownload.
set(ZLIB_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/zlib")
set(LIBAEC_SOURCE_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}/libaec")
set(HDF5_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/hdf5")
set(MBEDTLS_SOURCE_DIR		"${CMAKE_DOWNLOAD_DIRECTORY}/mbedtls")
set(LIBSSH2_SOURCE_DIR		"${CMAKE_DOWNLOAD_DIRECTORY}/libssh2")
set(CURL_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/curl")
set(CFITSIO_SOURCE_DIR		"${CMAKE_DOWNLOAD_DIRECTORY}/cfitsio")
set(HEALPIX_SOURCE_DIR		"${CMAKE_DOWNLOAD_DIRECTORY}/healpix")
set(CAMB_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/camb")
set(FFTW_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/fftw")
set(BLAS_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/blas")
#
#------------------------------------------------------------------------------
# output of the summary into the screen
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON INSTALLATION:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Projects will be downloaded into: ${CMAKE_DOWNLOAD_DIRECTORY}")
message(STATUS "Projects will be installed into: ${CMAKE_INSTALL_PREFIX}")
