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
# Description: This script determines the location of ZLIB on the host system.
# If it fails to do so, it will download, compile and install ZLIB from source.
# ZLIB is required to successfully compile HDF5 and other libraries.
#================================================================================

message(STATUS "---------------------------------------------------------------")
if(NOT (ZLIB_FORCE_COMPILE OR ALL_FORCE_COMPILE))
	find_package(ZLIB)
endif()

if(NOT ZLIB_FOUND)
	# Creating configure command for cURL
	set(curl_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>"
		)
	#------------------------------------------------------------------------------
	# Getting ZLIB from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}
		URL "${${project}_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		CONFIGURE_COMMAND "${${project}_configure_command}"
		)
	# getting curl directories
	ExternalProject_Get_Property(${project} source_dir)
  ExternalProject_Get_Property(${project} install_dir)
	# specifying curl libraries and binaries
	set(CURL_SOURCE_DIR ${source_dir})
  set(CURL_BINARY_DIR ${install_dir}/bin)
	set(CURL_INCLUDE_DIR ${install_dir}/include)#/${project})
  set(CURL_LIBRARIES ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	# including curl into a project
	#include_directories(${CURL_SOURCE_DIR})
	include_directories(${CURL_BINARY_DIR})
	include_directories(${CURL_INCLUDE_DIR})
	message(STATUS "Curl INCLUDE DIR will be ${CURL_INCLUDE_DIR}")
	#message(STATUS "Curl BINARY DIR will be ${CURL_BINARY_DIR}")
	#message(STATUS "Curl SOURCE DIR will be ${CURL_SOURCE_DIR}")
	message(STATUS "Curl LIBRARIES will be ${CURL_LIBRARIES}")

	# setting an environment variable for cfitsio to find curl library
	set(ENV{PATH} 
		#$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		$ENV{PATH}:${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
		)
	#message(STATUS "ENV PATH: " $ENV{PATH})
else()
	add_custom_target(${project} ALL "")
	include_directories(${ZLIB_INCLUDE_DIRS})
	message(STATUS "ZLIB LIBRARIES ARE: ${ZLIB_LIBRARIES}")
	message(STATUS "ZLIB INCLUDE DIRS ARE: ${ZLIB_INCLUDE_DIRS}")
	# setting an environment variable for cfitsio to find curl library
	set(ENV{PATH} 
		#$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		$ENV{PATH}
		)
endif()

