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
# PROJECT'S DOWNLOAD/INSTALL/OUTPUT DIRECTORIES 
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
# Subprojects' source directories. 
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
# Open Source alternatives 
set(FFTW_SOURCE_DIR				"${CMAKE_DOWNLOAD_DIRECTORY}/fftw")
set(OPENBLAS_SOURCE_DIR	  "${CMAKE_DOWNLOAD_DIRECTORY}/openblas")
# AMD AOCL
set(BLIS_SOURCE_DIR	      "${CMAKE_DOWNLOAD_DIRECTORY}/blis")
set(FLAME_SOURCE_DIR	    "${CMAKE_DOWNLOAD_DIRECTORY}/flame")
set(AMDFFTW_SOURCE_DIR  	"${CMAKE_DOWNLOAD_DIRECTORY}/amdfftw")
#------------------------------------------------------------------------------
# To have a summary of the installation, we need to have a summary of the Host 
# system. These variables serve this purpose.
#------------------------------------------------------------------------------
# CPU 
cmake_host_system_information(RESULT CPU_DESCRIPTION       QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT N_LOGICAL_CORES       QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT N_PHYSICAL_CORES      QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT CPU_HAS_SERIAL_NUMBER QUERY HAS_SERIAL_NUMBER)
cmake_host_system_information(RESULT CPU_SERIAL_NUMBER     QUERY PROCESSOR_SERIAL_NUMBER)
cmake_host_system_information(RESULT CPU_IS_64BIT          QUERY IS_64BIT)
cmake_host_system_information(RESULT CPU_HAS_FPU           QUERY HAS_FPU)
cmake_host_system_information(RESULT CPU_HAS_MMX           QUERY HAS_MMX)
# One if processor supports Ext. MMX instructions
cmake_host_system_information(RESULT CPU_HAS_MMX_PLUS      QUERY HAS_MMX_PLUS) 
cmake_host_system_information(RESULT CPU_HAS_SSE           QUERY HAS_SSE)
cmake_host_system_information(RESULT CPU_HAS_SSE2          QUERY HAS_SSE2)
cmake_host_system_information(RESULT CPU_HAS_SSE_FP        QUERY HAS_SSE_FP)
# One if processor supports SSE MMX instructions
cmake_host_system_information(RESULT CPU_HAS_SSE_MMX       QUERY HAS_SSE_MMX)
# Third-party package to identify the presence of SSE, AVX etc.
# Returns: 
# SSE2_FOUND, SSE3_FOUND, SSSE3_FOUND, SSE4_1_FOUND, SSE4_2_FOUND, 
# AVX_FOUND, AVX2_FOUND, AVX512_FOUND
find_package(SSE)
# To conform to the style above, creating new variables
if(SSE3_FOUND)
  set(CPU_HAS_SSE3   1)
else()
  set(CPU_HAS_SSE3   0)
endif()
if(SSSE3_FOUND)
  set(CPU_HAS_SSSE3  1)
else()
  set(CPU_HAS_SSSE3  0)
endif()
if(SSE4_1_FOUND)
  set(CPU_HAS_SSE4_1 1)
else()
  set(CPU_HAS_SSE4_1 0)
endif()
if(SSE4_2_FOUND)
  set(CPU_HAS_SSE4_2 1)
else()
  set(CPU_HAS_SSE4_2 0)
endif()
if(AVX_FOUND)
  set(CPU_HAS_AVX    1)
else()
  set(CPU_HAS_AVX    0)
endif()
if(AVX2_FOUND)
  set(CPU_HAS_AVX2   1)
else()
  set(CPU_HAS_AVX2   0)
endif()
if(AVX512_FOUND)
  set(CPU_HAS_AVX512 1)
else()
  set(CPU_HAS_AVX512 0)
endif()

message("SSE2_FOUND   ${SSE2_FOUND}")
message("SSE3_FOUND   ${SSE3_FOUND}")
message("SSSE3_FOUND  ${SSSE3_FOUND}")
message("SSE4_1_FOUND ${SSE4_1_FOUND}")
message("SSE4_2_FOUND ${SSE4_2_FOUND}")
message("AVX_FOUND    ${AVX_FOUND}")
message("AVX2_FOUND   ${AVX2_FOUND}")
message("AVX512_FOUND ${AVX512_FOUND}")
# OS information
cmake_host_system_information(RESULT HOST_OS_RELEASE       QUERY OS_RELEASE)
cmake_host_system_information(RESULT HOST_OS_VERSION       QUERY OS_VERSION)
cmake_host_system_information(RESULT HOST_OS_PLATFORM      QUERY OS_PLATFORM)
# RAM information
cmake_host_system_information(RESULT TOT_VIRTUAL_MEMORY    QUERY TOTAL_VIRTUAL_MEMORY)
cmake_host_system_information(RESULT AVAIL_VIRTUAL_MEMORY  QUERY AVAILABLE_VIRTUAL_MEMORY)
cmake_host_system_information(RESULT TOT_PHYSICAL_MEMORY   QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT AVAIL_PHYSICAL_MEMORY QUERY AVAILABLE_PHYSICAL_MEMORY)
#------------------------------------------------------------------------------
# If any problems with installation will occur, which cannot be fixed quickly,
# these variables will force a fresh installation for every specified library.
# If USE_SYSTEM_LIBS is set to FALSE, all libraries will be recompiled,
# whereas if set to TRUE, the libraries will first be searched on the
# system and only if not found, compiled from source. However, If USE_SYSTEM_LIBS
# is set to TRUE but, e.g. USE_SYSTEM_HDF5 is set to FALSE, then HDF5 will
# be compiled from source (i.e. the latter variable is having the advantage).
#------------------------------------------------------------------------------
option(USE_SYSTEM_LIBS    "Enables search for all LIBS on the system."    ON)
# BLAS/LAPACK
option(USE_SYSTEM_BLAS    "Enables search for BLAS/LAPACK on the system." ON)
# FFTW
option(USE_SYSTEM_FFTW    "Enables search for FFTW on the system."        ON) #OFF)
if(NOT FFTW_ENABLE_AVX)
  if(CPU_HAS_AVX)
    option(FFTW_ENABLE_AVX    "Enables AVX support for FFTW library"      ON)
  else()
    option(FFTW_ENABLE_AVX    "Enables AVX support for FFTW library"      OFF)
  endif()
endif()
if(NOT FFTW_ENABLE_AVX2)
  if(CPU_HAS_AVX2)
    option(FFTW_ENABLE_AVX2   "Enables AVX2 support for FFTW library"     ON)
  else()
    option(FFTW_ENABLE_AVX2   "Enables AVX2 support for FFTW library"     OFF)
  endif()
endif()
if(NOT FFTW_ENABLE_SSE)
  if(CPU_HAS_SSE)
    option(FFTW_ENABLE_SSE    "Enables SSE support for FFTW library"      ON)
  else()
    option(FFTW_ENABLE_SSE    "Enables SSE support for FFTW library"      OFF)
  endif()
endif()
if(NOT FFTW_ENABLE_SSE2)
  if(CPU_HAS_SSE2)
    option(FFTW_ENABLE_SSE2   "Enables SSE2 support for FFTW library"     ON)
  else()
    option(FFTW_ENABLE_SSE2   "Enables SSE2 support for FFTW library"     OFF)
  endif()
endif()
# HDF5
option(USE_SYSTEM_ZLIB    "Enables search for ZLIB on the system."        ON)
option(USE_SYSTEM_LIBAEC  "Enables search for LibAEC on the system."      ON)
option(USE_SYSTEM_HDF5    "Enables search for HDF5 on the system."        ON) #OFF)
# CFITSIO
option(USE_SYSTEM_MBEDTLS "Enables search for MbedTLS on the system."     ON)
option(USE_SYSTEM_LIBSSH2 "Enables search for LibSSH2 on the system."     ON)
option(USE_SYSTEM_CURL    "Enables search for cURL on the system."        ON)
option(USE_SYSTEM_CFITSIO "Enables search for CFITSIO on the system."     ON) #OFF)
# Can choose whether to compile CFITSIO with or without cURL support
option(CFITSIO_USE_CURL   "Installs CFITSIO with cURL support."           OFF)
# HEALPix
option(USE_SYSTEM_HEALPIX "Enables search for HEALPIX on the system."     ON) #OFF)
#------------------------------------------------------------------------------
# Commander3 can use various BLAS/LAPACK & FFT implementations from different 
# vendors (AMD, Intel etc.). The following variables define the particular 
# implementation to use -- the backend.
#------------------------------------------------------------------------------
set(COMM3_BACKEND "any"
  CACHE STRING
  "Defines which BLAS/LAPACK & FFTW implementation to use. 
  Possible values are: aocl, mkl, opensrc, any. Default: any."
  )
#------------------------------------------------------------------------------
