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
# Description: This script determines the location of MbedTLS on the host system.
# If it fails to do so, it will download, compile and install MbedTLS from source.
# MbedTLS is (not strictly) required by cURL. If either CFITSIO and/or cURL were 
# found, the script will go with system library.
#================================================================================

if(NOT (CFITSIO_FOUND AND CURL_FOUND) AND CFITSIO_USE_CURL)
	#message(STATUS "---------------------------------------------------------------")
	#if(USE_SYSTEM_MBEDTLS AND USE_SYSTEM_LIBS)
	#	find_package(MBEDTLS)
	#endif()

	if(COMPILE_MBEDTLS) 
		#------------------------------------------------------------------------------
		# Note: the explicit splitting for download and install step is done on purpose
		# to avoid errors when you want to recompile libraries for different owls etc.
		# In addition, this will allow us to download sources only once and then just 
		# reuse it whenever possible.
		#------------------------------------------------------------------------------
		# Getting MbedTLS from source.
		#------------------------------------------------------------------------------
		# Checking whether we have source directory and this directory is not empty.
		if(NOT EXISTS "${MBEDTLS_SOURCE_DIR}/CMakeLists.txt")
			message(STATUS "No MBEDTLS sources were found; thus, will download it from source:\n${mbedtls_git_url}")
			ExternalProject_Add(
				mbedtls_src
				GIT_REPOSITORY		"${mbedtls_git_url}"
				GIT_TAG						"${mbedtls_git_tag}"
				PREFIX						"${LIBS_BUILD_DIR}"
				DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
				SOURCE_DIR				"${MBEDTLS_SOURCE_DIR}"
				LOG_DIR						"${CMAKE_LOG_DIR}"
				LOG_DOWNLOAD			ON
				CONFIGURE_COMMAND ""
				BUILD_COMMAND			""
				INSTALL_COMMAND		""
				)
		else()
			message(STATUS "Found an existing MBEDTLS sources inside:\n${MBEDTLS_SOURCE_DIR}")
			add_custom_target(mbedtls_src
				ALL ""
				)
		endif()
		#------------------------------------------------------------------------------
		# Getting MbedTLS from source and compiling both static and shared libraries.
		#------------------------------------------------------------------------------
		ExternalProject_Add(
			mbedtls
			DEPENDS						required_libraries
												zlib
												mbedtls_src
			PREFIX						"${LIBS_BUILD_DIR}"
			SOURCE_DIR				"${MBEDTLS_SOURCE_DIR}"
			INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_CONFIGURE			ON
			LOG_BUILD					ON
			LOG_INSTALL				ON
			# commands how to build the project
			DOWNLOAD_COMMAND ""
			CMAKE_ARGS
				-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
				# Specifying installations paths for binaries and libraries
				-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
				# ensuring that the lib name would be lib and not lib64
				#-DLIB_INSTALL_DIR=lib
				-DLIB_INSTALL_DIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
				# Specifying compilers
				-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
				# There are some problems with generated programs
				# everythign works, but if user wants to use MbedTLS
				# he/she needs to add *.so into their PATH.	
				-DENABLE_PROGRAMS:BOOL=ON
				-DENABLE_TESTING:BOOL=ON
				-DINSTALL_MBEDTLS_HEADERS:BOOL=ON
				-DUSE_STATIC_MBEDTLS_LIBRARY:BOOL=ON
				-DUSE_SHARED_MBEDTLS_LIBRARY:BOOL=ON
				# Enabling ZLIB support
				-DENABLE_ZLIB_SUPPORT:BOOL=ON
				-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
				-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}
			)
		#------------------------------------------------------------------------------
		# Note: we need to define all these libraries separately to sucessfully compile
		# LibSSH2, as it fails to detect MbedTLS.
		# Use static library linking otherwise we need to add these to LD_LIBRARY_PATH.
		# Note2: In this case getitng libmbedcrypto.a(bignum.c.o): relocation R_X86_64_32S 
		# against `.rodata' can not be used when making a shared object; recompile with -fPIC
		set(MBEDTLS_LIBRARY
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}mbedtls${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}mbedtls${CMAKE_STATIC_LIBRARY_SUFFIX}" 
			)
		set(MBEDX509_LIBRARY
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}mbedx509${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}mbedx509${CMAKE_STATIC_LIBRARY_SUFFIX}" 
			)
		set(MBEDCRYPTO_LIBRARY
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}mbedcrypto${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}mbedcrypto${CMAKE_STATIC_LIBRARY_SUFFIX}" 
			)
		set(MBEDTLS_LIBRARIES
			"${MBEDTLS_LIBRARY}"
			"${MBEDX509_LIBRARY}"
			"${MBEDCRYPTO_LIBRARY}"
			)
		set(MBEDTLS_INCLUDE_DIRS
			"${CMAKE_INSTALL_PREFIX}/include"
			)
		include_directories(${MBEDTLS_INCLUDE_DIRS})
		#------------------------------------------------------------------------------
		#message(STATUS "MbedTLS LIBRARIES will be: ${MBEDTLS_LIBRARIES}")
		#message(STATUS "MbedTLS INCLUDE DIRS will be: ${MBEDTLS_INCLUDE_DIRS}")
		#------------------------------------------------------------------------------
	else()
		# If mbedtls exists on the system, we just use this version instead.
		add_custom_target(mbedtls ALL "")
		#------------------------------------------------------------------------------
		#message(STATUS "MbedTLS LIBRARIES are: ${MBEDTLS_LIBRARIES}")
		#message(STATUS "MbedTLS INCLUDE DIRS are: ${MBEDTLS_INCLUDE_DIRS}")
		#------------------------------------------------------------------------------
	endif()
else()
	# making an empty target so the project will compile regardless of mbedtls
	add_custom_target(${project} ALL "")
endif()
