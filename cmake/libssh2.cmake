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
# Description: This script determines the location of LibSSH2 on the host system.
# If it fails to do so, it will download, compile and install LibSSH2 from source.
# LibSSH2 is (not strictly) required by cURL.
#================================================================================

if(NOT (CFITSIO_FOUND AND CURL_FOUND) AND CFITSIO_USE_CURL)
	message(STATUS "---------------------------------------------------------------")
	#if(NOT (LIBSSH2_FORCE_COMPILE OR ALL_FORCE_COMPILE))
	if(USE_SYSTEM_LIBSSH2 AND USE_SYSTEM_LIBS)
		find_package(LibSSH2)
	endif()

	if(NOT LIBSSH2_FOUND) 
		#------------------------------------------------------------------------------
		# Getting LibSSH2 from source.
		#------------------------------------------------------------------------------
		# Checking whether we have source directory and this directory is not empty.
		if(NOT EXISTS "${LIBSSH2_SOURCE_DIR}/CMakeLists.txt")
			message(STATUS "No LIBSSH2 sources were found; thus, will download it from source:\n${libssh2_git_url}")
			ExternalProject_Add(
				libssh2_src
				GIT_REPOSITORY		"${libssh2_git_url}"
				GIT_TAG						"${libssh2_git_tag}"
				PREFIX						"${LIBS_BUILD_DIR}"
				SOURCE_DIR				"${LIBSSH2_SOURCE_DIR}"
				DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
				LOG_DIR						"${CMAKE_LOG_DIR}"
				LOG_DOWNLOAD			ON
				CONFIGURE_COMMAND ""
				BUILD_COMMAND			""
				INSTALL_COMMAND		""
				)
		else()
			message(STATUS "Found an existing LIBSSH2 sources inside:\n${LIBSSH2_SOURCE_DIR}")
			add_custom_target(libssh2_src
				ALL ""
				)
		endif()
		#------------------------------------------------------------------------------
		# Building both -- static and shared -- versions of the library.
		# Building is done consequently, i.e. first it will install static
		# and then shared libraries. It is done in this way, because there
		# are some error occuring from time to time while compiling/installing
		# some parts of the library.	
		#------------------------------------------------------------------------------
		# Building Static LibSSH2
		#------------------------------------------------------------------------------
		ExternalProject_Add(
			libssh2_static
			DEPENDS						zlib 
												mbedtls
												libssh2_src
			PREFIX						"${LIBS_BUILD_DIR}"
			SOURCE_DIR				"${LIBSSH2_SOURCE_DIR}"
			INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_CONFIGURE			ON
			LOG_BUILD					ON
			LOG_INSTALL				ON
			# Disabling download
			DOWNLOAD_COMMAND	""
			# commands how to build the project
			CMAKE_ARGS
				-DCMAKE_BUILD_TYPE=Release
				# Specifying installations paths for binaries and libraries
				-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
				-DCMAKE_INSTALL_LIBDIR=lib
				# Specifying compilers
				-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
				# Building static library with MbedTLS backend
				-DBUILD_SHARED_LIBS=OFF
				#-D$ENV{mbedTLS_ROOT}=$ENV{MBEDTLS_ROOT}
				-DCRYPTO_BACKEND:STRING=mbedTLS
				# MbedTLS comes with its own FindmbedTLS.cmake file
				# and if this module fails to find MbedTLS on the
				# system, it will crush. So, we must artificially
				# provide below variables to mimick FindmbedTLS
				# behavior.	
				-DMBEDTLS_FOUND=TRUE
				#-DMBEDTLS_LIBRARIES:FILEPATH=${MBEDTLS_LIBRARIES}
				-DMBEDTLS_LIBRARY:FILEPATH=${MBEDTLS_LIBRARY}
				-DMBEDX509_LIBRARY:FILEPATH=${MBEDX509_LIBRARY}
				-DMBEDCRYPTO_LIBRARY:FILEPATH=${MBEDCRYPTO_LIBRARY}
				-DMBEDTLS_INCLUDE_DIR:PATH=${MBEDTLS_INCLUDE_DIRS}
				#-DMBEDTLS_LIBRARY_DIR:FILEPATH=${MBEDTLS_LIBRARIES}
				# Enabling ZLIB support
				-DENABLE_ZLIB_COMPRESSION:BOOL=ON
				-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
				-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}
			)
		#------------------------------------------------------------------------------
		# Building Shared LibSSH2
		#------------------------------------------------------------------------------
		ExternalProject_Add(
			libssh2_shared
			DEPENDS						zlib 
												mbedtls
												libssh2_src
												libssh2_static
			PREFIX						"${LIBS_BUILD_DIR}"
			SOURCE_DIR				"${LIBSSH2_SOURCE_DIR}"
			INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_CONFIGURE			ON
			LOG_BUILD					ON
			LOG_INSTALL				ON
			# Disabling download
			DOWNLOAD_COMMAND	""
			# commands how to build the project
			CMAKE_ARGS
				-DCMAKE_BUILD_TYPE=Release
				# Specifying installations paths for binaries and libraries
				-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
				-DCMAKE_INSTALL_LIBDIR=lib
				# Specifying compilers
				-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
				# Building shared library with MbedTLS backend
				-DBUILD_SHARED_LIBS=ON
				-DCRYPTO_BACKEND:STRING=mbedTLS
				# MbedTLS comes with its own FindmbedTLS.cmake file
				# and if this module fails to find MbedTLS on the
				# system, it will crush. So, we must artificially
				# provide below variables to mimick FinmbedTLS.cmake
				# behavior.	
				-DMBEDTLS_FOUND=TRUE
				#-DMBEDTLS_LIBRARIES:FILEPATH=${MBEDTLS_LIBRARIES}
				-DMBEDTLS_LIBRARY:FILEPATH=${MBEDTLS_LIBRARY}
				-DMBEDX509_LIBRARY:FILEPATH=${MBEDX509_LIBRARY}
				-DMBEDCRYPTO_LIBRARY:FILEPATH=${MBEDCRYPTO_LIBRARY}
				-DMBEDTLS_INCLUDE_DIR:PATH=${MBEDTLS_INCLUDE_DIRS}
				# Enabling ZLIB support
				-DENABLE_ZLIB_COMPRESSION:BOOL=ON
				-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
				-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}
			)
		#------------------------------------------------------------------------------
		set(LIBSSH2_LIBRARY
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}ssh2${CMAKE_SHARED_LIBRARY_SUFFIX}" 
			)
		set(LIBSSH2_INCLUDE_DIR
			"${CMAKE_INSTALL_PREFIX}/include"
			)
		include_directories(${LIBSSH2_INCLUDE_DIR})
		# Adding new target -- libssh2 -- to ensure that only after all libraries built
		# the project can use this target.
		add_custom_target(libssh2 
			ALL ""
			DEPENDS libssh2_static
							libssh2_shared
			)
		#------------------------------------------------------------------------------
		message(STATUS "LibSSH2 LIBRARY will be: ${LIBSSH2_LIBRARY}")
		message(STATUS "LibSSH2 INCLUDE DIR will be: ${LIBSSH2_INCLUDE_DIR}")
		#------------------------------------------------------------------------------
	else()
		# If libssh2 exists on the system, we just use this version instead.
		add_custom_target(${project} ALL "")
		#------------------------------------------------------------------------------
		message(STATUS "LibSSH2 LIBRARY are: ${LIBSSH2_LIBRARY}")
		message(STATUS "LibSSH2 INCLUDE DIR are: ${LIBSSH2_INCLUDE_DIR}")
		#------------------------------------------------------------------------------
	endif()
else()
	# making an empty target so the project will compile regardless of libssh2
	add_custom_target(libssh2 ALL "")
endif()
