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
# Description: This script determines the location of cURL on the host system.
# If it fails to do so, it will download, compile and install cURL from source.
# Together with CFitsio, cURL is required to successfully compile HEALPix.
#================================================================================

if(NOT CFITSIO_FOUND AND CFITSIO_USE_CURL)
	#message(STATUS "---------------------------------------------------------------")
	## looking for cURL in the system. 
	#if(USE_SYSTEM_CURL AND USE_SYSTEM_LIBS)
	#	# CMake configure scripts (in versions 7.69-7.74) doesn't work properly,
	#	# so we look for cURL in a standard manner.
	#	set(CURL_NO_CURL_CMAKE ON)
	#	find_package(CURL 7.54)
	#endif()

	# TODO: to properly link zlib to curl (cmake build) need to export it in .bashrc (or in $ENV{ZLIB_ROOT})
	# export ZLIB_ROOT=/path/to/zlib/root, e.g. ${CMAKE_INSTALL_PREFIX}
	# It seems that CMake build is still lacking compared to ./configure
	# even though it is pretty good for version 74. The reason for it is
	# that cURL doesn't have option to specify paths to OpenSSl, ZLib and
	# SSH2 manually, but it uses its own find_package(s) which is good.
	# Note2: Switched to version 7.78, which (should) have even better
	# CMake support.
	if(COMPILE_CURL)
		# Creating configure command for cURL
		# Note: if checked downloaded fron git it needs autoreconf -fi
		# to create configure first
		set(ENV{PATH} 
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}:$ENV{PATH}
			)
		set(ENV{LD_LIBRARY_PATH}
			${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$ENV{LD_LIBRARY_PATH}
			)
		#------------------------------------------------------------------------------
		# Note: the explicit splitting for download and install step is done on purpose
		# to avoid errors when you want to recompile libraries for different owls etc.
		# In addition, this will allow us to download sources only once and then just 
		# reuse it whenever possible.
		#------------------------------------------------------------------------------
		# Getting cURL from source
		#------------------------------------------------------------------------------
		# Checking whether we have source directory and this directory is not empty.
		if(NOT EXISTS "${CURL_SOURCE_DIR}/CMakeLists.txt")
			message(STATUS "No CURL sources were found; thus, will download it from source:\n${curl_git_url}")
			ExternalProject_Add(
				curl_src
				GIT_REPOSITORY		"${curl_git_url}"
				GIT_TAG						"${curl_git_tag}"
				PREFIX						"${LIBS_BUILD_DIR}"
				DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
				SOURCE_DIR				"${CURL_SOURCE_DIR}"
				LOG_DIR						"${CMAKE_LOG_DIR}"
				LOG_DOWNLOAD			ON
				# commands how to build the project
				CONFIGURE_COMMAND ""
				BUILD_COMMAND			""
				INSTALL_COMMAND		""
				)
		else()
			message(STATUS "Found an existing CURL sources inside:\n${CURL_SOURCE_DIR}")
			add_custom_target(curl_src
				ALL ""
				)
		endif()
		#------------------------------------------------------------------------------
		# Compiling and Installing Static and Shared cURL 
		#------------------------------------------------------------------------------
		list(APPEND _CURL_LIB_TYPES_ static shared)
		list(APPEND _CURL_LIB_BOOL_VALS_ -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON)
		foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _CURL_LIB_TYPES_ _CURL_LIB_BOOL_VALS_)
			ExternalProject_Add(
				curl_${_lib_type_}
				DEPENDS						required_libraries
													zlib 
													mbedtls
													libssh2 
													curl_src
				PREFIX						"${LIBS_BUILD_DIR}"
				SOURCE_DIR				"${CURL_SOURCE_DIR}"
				INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
				LOG_DIR						"${CMAKE_LOG_DIR}"
				LOG_CONFIGURE			ON
				LOG_BUILD					ON
				LOG_INSTALL				ON
				# Disabling download
				DOWNLOAD_COMMAND	""
				# commands how to build the project
				CMAKE_ARGS
          -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
					# Specifying installations paths for binaries and libraries
					-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
					#-DCMAKE_INSTALL_LIBDIR:PATH=lib
					-DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
					-DCMAKE_INSTALL_INCLUDEDIR:PATH=include
					-DBUILD_CURL_EXE:BOOL=ON
					# Building both static and shared libraries
					${_bool_val_}
					# Specifying compilers
					-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
					-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
					-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
					# Specifying the location of MbedTLS library
					-DCMAKE_USE_MBEDTLS:BOOL=ON
					-DMBEDTLS_LIBRARY:FILEPATH=${MBEDTLS_LIBRARY}
					-DMBEDX509_LIBRARY:FILEPATH=${MBEDX509_LIBRARY}
					-DMBEDCRYPTO_LIBRARY:FILEPATH=${MBEDCRYPTO_LIBRARY}
					# Specifying the location of OpenSSL library
					#-DCMAKE_USE_OPENSSL:BOOL=ON
					#-DOPENSSL_INCLUDE_DIR:PATH=${OPENSSL_INCLUDE_DIR}
					#-DOPENSSL_SSL_LIBRARY:FILEPATH=${OPENSSL_SSL_LIBRARIES}
					#-DOPENSSL_CRYPTO_LIBRARY:FILEPATH=${OPENSSL_CRYPTO_LIBRARIES}
					# Specyfing the location of ZLIB library
					-DCURL_ZLIB:BOOL=ON
					-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
					-DZLIB_LIBRARY_RELEASE:FILEPATH=${ZLIB_LIBRARIES}
					# Specyfing the location of LibSSH2
					-DLIBSSH2_LIBRARY:PATH=${LIBSSH2_LIBRARY}
					-DLIBSSH2_INCLUDE_DIR:FILEPATH=${LIBSSH2_INCLUDE_DIR}
			)
		endforeach()
		#------------------------------------------------------------------------------
		# Creating Unified Target
		#------------------------------------------------------------------------------
		add_custom_target(curl 
			ALL ""
			DEPENDS curl_static
							curl_shared
			)
		#------------------------------------------------------------------------------
		# getting curl directories
		set(CURL_INCLUDE_DIR
			"${CMAKE_INSTALL_PREFIX}/include"
			)
		# Using static linking otherwise we need to put things into LD_LIBRARY_PATH
		# libcurl.a(mime.c.o): relocation R_X86_64_32S against `base64' can not be 
		# used when making a shared object; recompile with -fPIC
		# Note: Use shared linking otherwise will get a lot of linking problems when 
		# compiling commander
		set(CURL_LIBRARIES
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}curl${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}curl${CMAKE_STATIC_LIBRARY_SUFFIX}" 
			"${MBEDTLS_LIBRARIES}"
			)
		# including curl into a project
		include_directories(${CURL_INCLUDE_DIR})
		#------------------------------------------------------------------------------
		#message(STATUS "cURL LIBRARIES will be ${CURL_LIBRARIES}")
		#message(STATUS "cURL INCLUDE DIR will be ${CURL_INCLUDE_DIR}")
		#------------------------------------------------------------------------------
	else()
		add_custom_target(curl ALL "")
		include_directories(${CURL_INCLUDE_DIR})
		include_directories(${CURL_BINARY_DIR})
		#------------------------------------------------------------------------------
		#message(STATUS "cURL LIBRARIES are ${CURL_LIBRARIES}")
		#message(STATUS "cURL INCLUDE DIR is ${CURL_INCLUDE_DIR}")
		#------------------------------------------------------------------------------
	endif()
else()
	# making an empty target so the project will compile regardless of cURL 
	add_custom_target(curl ALL "")
endif()
