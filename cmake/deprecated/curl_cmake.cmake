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
# Description: This script determines the location of cURL on the host system.
# If it fails to do so, it will download, compile and install cURL from source.
# Together with CFitsio, cURL is required to successfully compile HEALPix.
#================================================================================

message(STATUS "---------------------------------------------------------------")
# looking for curl in the system and download it if it is not present
if(NOT (CURL_FORCE_COMPILE OR ALL_FORCE_COMPILE))
	set(CURL_NO_CURL_CMAKE ON)
	find_package(CURL)
endif()

if(NOT CURL_FOUND)
	#------------------------------------------------------------------------------
	# Getting cURL from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}
		DEPENDS required_libraries zlib 
		URL "${${project}_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		#LOG_CONFIGURE ON
		#LOG_BUILD ON
		#LOG_INSTALL ON
		# commands how to build the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_INSTALL_LIBDIR:PATH=lib
			-DCMAKE_INSTALL_INCLUDEDIR:PATH=include
			-DBUILD_CURL_EXE:BOOL=ON
			-DBUILD_SHARED_LIBS:BOOL=ON
			-DCMAKE_USE_OPENSSL:BOOL=ON
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${COMMANDER3_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${COMMANDER3_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${COMMANDER3_C_COMPILER}
			# Specifying the location of OpenSSL library
			-DOPENSSL_INCLUDE_DIR:PATH=${OPENSSL_INCLUDE_DIR}
			-DOPENSSL_SSL_LIBRARY:FILEPATH=${OPENSSL_SSL_LIBRARIES}
			-DOPENSSL_CRYPTO_LIBRARY:FILEPATH=${OPENSSL_CRYPTO_LIBRARIES}
			# Specyfing the location of ZLIB library
			-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
			-DZLIB_LIBRARY_RELEASE:FILEPATH=${ZLIB_LIBRARIES}
			-DCURL_ZLIB:BOOL=ON
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
	include_directories(${CURL_INCLUDE_DIR})
	include_directories(${CURL_BINARY_DIR})
	message(STATUS "Curl INCLUDE DIR is ${CURL_INCLUDE_DIR}")
	#message(STATUS "Curl BINARY DIR is ${CURL_BINARY_DIR}")
	#message(STATUS "Curl SOURCE DIR is ${CURL_SOURCE_DIR}")
	message(STATUS "Curl LIBRARIES are ${CURL_LIBRARIES}")
	# setting an environment variable for cfitsio to find curl library
	set(ENV{PATH} 
		#$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		$ENV{PATH}
		)
	#message(STATUS "ENV PATH: " $ENV{PATH})
	# Healpix complains about curl library - it needs to be in the same location as cfitsio
	add_custom_command(
		TARGET ${project} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CURL_LIBRARIES} 
				${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}curl${CMAKE_SHARED_LIBRARY_SUFFIX}
				)
endif()

