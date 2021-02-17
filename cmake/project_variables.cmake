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
set(CMAKE_LIBRARY64_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib64"
	CACHE STRING
	"Directory where to install all the libraries for 64."
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
set(DOXYGEN_BUILD_DOCS OFF 
	CACHE BOOL
	"Determine whether to use doxygen or not."
	)
#------------------------------------------------------------------------------
# If any problems with installation will occur, which cannot be fixed quickly,
# these variables will force a fresh installation for every specified library.
# If LIBSALL_FORCE_COMPILE is set to TRUE, all libraries will be recompiled,
# whereas if set to FALSE, the libraries will first be searched on the
# system and only if not found, compiled from source. If LIBSALL_FORCE_COMPILE
# is set to FALSE but, e.g. HDF5_FORCE_COMPILE is set to TRUE, then HDF5 will
# be compiled from source (it will be given the advantage).
#------------------------------------------------------------------------------
#set(ALL_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of all Commander3 dependencies."
#  )
#set(BLAS_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of OpenBLAS."
#  )
#set(HDF5_FORCE_COMPILE FALSE
#  CACHE BOOL
#  "Forces fresh installation of HDF5."
#  )
#set(FFTW_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of FFTW."
#  )
#set(CFITSIO_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of CFITSIO."
#  )
#set(HEALPIX_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of HEALPIX."
#  )
#set(CURL_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of CURL."
#  )
#set(DOXYGEN_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of DOXYGEN."
#  )
#set(FLEX_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of FLEX."
#  )
#set(BISON_FORCE_COMPILE FALSE
#  CACHE BOOL
#	"Forces fresh installation of BISON."
#  )
# All options default to OFF
#option(ALL_FORCE_COMPILE "Forces fresh installation of all Commander3 dependencies.")
#option(BLAS_FORCE_COMPILE "Forces fresh installation of OpenBLAS.")
#
#option(ZLIB_FORCE_COMPILE "Forces fresh installation of ZLib.")
#option(LIBAEC_FORCE_COMPILE "Forces fresh installation of LibAEC.")
#option(HDF5_FORCE_COMPILE "Forces fresh installation of HDF5.")
#
#option(FFTW_FORCE_COMPILE "Forces fresh installation of FFTW.")
# TODO: Add option "CFITSIO_USE_CURL" to give user an opportunity 
# to compile CFitsIO with cURL support. In this case, will download,
# compile and install ZLib -> MbedTLS -> LibSSH2 -> cURL.
#option(MBEDTLS_FORCE_COMPILE "Forces fresh installation of MbedTLS.")
#option(LIBSSH2_FORCE_COMPILE "Forces fresh installation of LibSSH2.")
#option(CURL_FORCE_COMPILE "Forces fresh installation of cURL.")
#option(CFITSIO_FORCE_COMPILE "Forces fresh installation of CFITSIO.")
#option(CFITSIO_USE_CURL "Installs CFITSIO with cURL support.")
#option(HEALPIX_FORCE_COMPILE "Forces fresh installation of HEALPix.")
#
#option(FLEX_FORCE_COMPILE "Forces fresh installation of FLEX.")
#option(BISON_FORCE_COMPILE "Forces fresh installation of BISON.")
#option(DOXYGEN_FORCE_COMPILE "Forces fresh installation of Doxygen.")
#
#cmake_dependent_option(BLAS_FORCE_COMPILE "Forces fresh installation of OpenBLAS." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(ZLIB_FORCE_COMPILE "Forces fresh installation of ZLib." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(HDF5_FORCE_COMPILE "Forces fresh installation of HDF5." ON "ALL_FORCE_COMPILE" ON)
#
#cmake_dependent_option(FFTW_FORCE_COMPILE "Forces fresh installation of FFTW." ON "ALL_FORCE_COMPILE" ON)
#
#cmake_dependent_option(CURL_FORCE_COMPILE "Forces fresh installation of cURL." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(CFITSIO_FORCE_COMPILE "Forces fresh installation of CFITSIO." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(HEALPIX_FORCE_COMPILE "Forces fresh installation of HEALPix." ON "ALL_FORCE_COMPILE" ON)
#
#cmake_dependent_option(FLEX_FORCE_COMPILE "Forces fresh installation of FLEX." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(BISON_FORCE_COMPILE "Forces fresh installation of BISON." ON "ALL_FORCE_COMPILE" ON)
#cmake_dependent_option(DOXYGEN_FORCE_COMPILE "Forces fresh installation of Doxygen." ON "ALL_FORCE_COMPILE" ON)
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
# Commander source dir
set(COMMANDER3_SOURCE_DIR "${CMAKE_SOURCE_DIR}/commander3/src")
# tempita source dir
set(TEMPITA_DIR ${CMAKE_SOURCE_DIR}/commander3/python)
# adding custom cmake modules directory, e.g. for FindSomething.cmake
# Note: It should be already inside root CmakeLists.txt, so 
# don't need to include in here
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
#------------------------------------------------------------------------------
# output of the summary into the screen
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON INSTALLATION:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Projects will be downloaded into: ${CMAKE_DOWNLOAD_DIRECTORY}")
message(STATUS "Projects will be installed into: ${CMAKE_INSTALL_PREFIX}")
