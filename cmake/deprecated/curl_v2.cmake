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

if(CFITSIO_USE_CURL)
	message(STATUS "---------------------------------------------------------------")
	# looking for curl in the system and download it if it is not present
	if(NOT (CURL_FORCE_COMPILE OR ALL_FORCE_COMPILE))
		find_package(CURL)
	endif()

	# TODO: to properly link zlib to curl (cmake build) need to export it in .bashrc (or in $ENV{ZLIB_ROOT})
	# export ZLIB_ROOT=/path/to/zlib/root, e.g. ${CMAKE_INSTALL_PREFIX}
	# It seems that CMake build is still lacking compared to ./configure
	# even though it is pretty good for version 74. The reason for it is
	# that cURL doesn't have option to specify paths to OpenSSl, ZLib and
	# SSH2 manually, but it uses its own find_package(s) which is good.
	if(NOT CURL_FOUND)
		# Creating configure command for cURL
		# Note: if checked downloaded fron git it needs autoreconf -fi
		# to create configure first
		set(ENV{PATH} 
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}:$ENV{PATH}
			)
		set(ENV{LD_LIBRARY_PATH}
			${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$ENV{LD_LIBRARY_PATH}
			)
		set(curl_configure_command 
			"autoreconf" "-fi" "&&"
			"${CMAKE_COMMAND}" "-E" "env" 
			"PATH=$ENV{PATH}"
			"LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}"
			"CPPFLAGS="-I/path/to/ssl/include"
			"FC=${COMMANDER3_Fortran_COMPILER}" 
			"CXX=${COMMANDER3_CXX_COMPILER}" 
			"CC=${COMMANDER3_C_COMPILER}" 
			"CPP=${COMMANDER3_CPP_COMPILER}" 
			#"CPPFLAGS=-I${OPENSSL_INCLUDE_DIR}"
			#"LDFLAGS=-L${OPENSSL_LIBRARIES}"
			"./configure" 
			"--prefix=<INSTALL_DIR>"
			"--with-ssl"
			"--with-zlib" #=PATH_VAR"
			"--with-libssh2=${LIBSSH2_LIBRARY}"
			)
		#------------------------------------------------------------------------------
		# Getting cURL from source
		#------------------------------------------------------------------------------
		ExternalProject_Add(${project}_src
			DEPENDS zlib
			GIT_REPOSITORY "${${project}_git_url}"
			GIT_TAG "${${project}_git_tag}"
			# PREFIX should be present, otherwise it will pull it into "build" dir
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
			LOG_DIR "${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD ON
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
			)
		#------------------------------------------------------------------------------
		ExternalProject_Add(${project}
			DEPENDS zlib libssh2 ${project}_src
			#URL "${${project}_url}"
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
			#DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}" #"${download_dir}"
			SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_src"
			# Binary dir is for ./configure
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_src"
			INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" #"${out_install_dir}"
			LOG_DIR "${CMAKE_LOG_DIR}"
			#LOG_DOWNLOAD ON
			LOG_CONFIGURE OFF
			LOG_BUILD ON
			LOG_INSTALL ON
			# Disabling download
			DOWNLOAD_COMMAND ""
			CONFIGURE_COMMAND "${${project}_configure_command}"
			BUILD_ALWAYS FALSE
			)
		# getting curl directories
		ExternalProject_Get_Property(${project} source_dir)
		ExternalProject_Get_Property(${project} install_dir)
		# specifying curl libraries and binaries
		set(CURL_SOURCE_DIR ${source_dir})
		set(CURL_BINARY_DIR ${install_dir}/bin)
		set(CURL_INCLUDE_DIR #${install_dir}/include)#/${project})
			"${CMAKE_INSTALL_PREFIX}/include"
			)
		#set(CURL_LIBRARIES ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
		set(CURL_LIBRARIES
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}curl${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			)
		# including curl into a project
		#include_directories(${CURL_SOURCE_DIR})
		include_directories(${CURL_BINARY_DIR})
		include_directories(${CURL_INCLUDE_DIR})
		message(STATUS "Curl LIBRARIES will be ${CURL_LIBRARIES}")
		message(STATUS "Curl INCLUDE DIR will be ${CURL_INCLUDE_DIR}")
		#message(STATUS "Curl BINARY DIR will be ${CURL_BINARY_DIR}")
		#message(STATUS "Curl SOURCE DIR will be ${CURL_SOURCE_DIR}")

		# setting an environment variable for cfitsio to find curl library
		#set(ENV{PATH} 
		#	#$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		#	$ENV{PATH}:${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
		#	)
		#message(STATUS "ENV PATH: " $ENV{PATH})
	else()
		add_custom_target(${project} ALL "")
		include_directories(${CURL_INCLUDE_DIR})
		include_directories(${CURL_BINARY_DIR})
		message(STATUS "Curl LIBRARIES are ${CURL_LIBRARIES}")
		message(STATUS "Curl INCLUDE DIR is ${CURL_INCLUDE_DIR}")
		#message(STATUS "Curl BINARY DIR is ${CURL_BINARY_DIR}")
		#message(STATUS "Curl SOURCE DIR is ${CURL_SOURCE_DIR}")
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
else()
	# making an empty target so the project will compile regardless of libssh2
	add_custom_target(${project} ALL "")
endif()
