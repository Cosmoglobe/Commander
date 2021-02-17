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
# Description: This script determines the location of ZLib on the host system.
# If it fails to do so, it will download, compile and install ZLib from source.
# ZLib is (not strictly) required by OpenSSL, LibSSH2, cURL, and HDF5.
#================================================================================

message(STATUS "---------------------------------------------------------------")

#if(NOT (ZLIB_FORCE_COMPILE OR HDF5_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#if(USE_SYSTEM_ZLIB AND NOT USE_SYSTEM_HDF5 AND USE_SYSTEM_LIBS)
if(USE_SYSTEM_ZLIB AND USE_SYSTEM_LIBS)
	#set(zlib_minimal_accepted_version "1.2.11")
	find_package(ZLIB 1.2.11)
	# Require ZLib to be of the most recent version
	#if(ZLIB_VERSION_STRING VERSION_LESS_EQUAL zlib_minimal_accepted_version)
	#	message(STATUS "Required version -- ${zlib_minimal_accepted_version} -- will be compiled from source.")
	#endif()
endif()

if(NOT ZLIB_FOUND) #OR (ZLIB_VERSION_STRING VERSION_LESS_EQUAL zlib_minimal_accepted_version)) 
	message(STATUS "Required version -- 1.2.11 -- will be compiled from source.")
	#------------------------------------------------------------------------------
	# Getting ZLib from source.
	#------------------------------------------------------------------------------
	ExternalProject_Add(zlib
		DEPENDS required_libraries
		URL "${zlib_url}"
		URL_MD5 "${zlib_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/zlib"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/zlib/src/zlib"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# commands how to build the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			# Specifying compilers
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		)
	#------------------------------------------------------------------------------
	# Defining the variable which will show the path to the compiled libraries
	set(ZLIB_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include"
		)
	include_directories(${ZLIB_INCLUDE_DIRS}) 
	set(ZLIB_LIBRARIES
		#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}z${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}z${CMAKE_SHARED_LIBRARY_SUFFIX}" 
		)
	# Setting ZLIB_ROOT environment variable to point future HDF5 installation to zlib libraries.
	# One of the requirements to correctly detect ZLib with CMake implementation of HDF5.
	set($ENV{ZLIB_ROOT}
		"${CMAKE_INSTALL_PREFIX}"
		)
	#------------------------------------------------------------------------------
	message(STATUS "ZLIB LIBRARIES will be: ${ZLIB_LIBRARIES}")
	message(STATUS "ZLIB INCLUDE DIRS will be: ${ZLIB_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	add_custom_target(zlib ALL "")
	#------------------------------------------------------------------------------
	message(STATUS "ZLIB LIBRARIES are: ${ZLIB_LIBRARIES}")
	message(STATUS "ZLIB INCLUDE DIRS are: ${ZLIB_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
