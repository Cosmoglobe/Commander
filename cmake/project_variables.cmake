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
# Where to output libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/lib"
	CACHE STRING
	"Directory where to install all the libraries."
	)
# Where to output executable(s)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_PREFIX}/bin"
	CACHE STRING
	"Directory where to install all the executables."
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
# whereas iof i set to FALSE, the libraries will first be searched on the
# system and only if not found, compile dfrom source. If LIBSALL_FORCE_COMPILE
# is set to FALSE but, e.g. HDF5_FORCE_COMPILE is set to TRUE, then HDF5 will
# be compile dfrom source (it will be given the advantage).
#------------------------------------------------------------------------------
set(LIBSALL_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of Commander3 dependencies."
  )
set(HDF5_FORCE_COMPILE FALSE
  CACHE BOOL
  "Forces fresh installation of HDF5."
  )
set(FFTW_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of FFTW."
  )
set(CFITSIO_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of CFITSIO."
  )
set(HEALPIX_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of HEALPIX."
  )
set(CURL_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of CURL."
  )
set(DOXYGEN_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of DOXYGEN."
  )
set(FLEX_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of FLEX."
  )
set(BISON_FORCE_COMPILE FALSE
  CACHE BOOL
	"Forces fresh installation of BISON."
  )

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
