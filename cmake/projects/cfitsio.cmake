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
# Description: This script determines the location of CFitsio on the host system.
# If it fails to do so, it will download, compile and install CFitsio from source.
# Together with cURL, CFitsio is required to successfully compile HEALPix.
#================================================================================

#message(STATUS "---------------------------------------------------------------")
#if(USE_SYSTEM_CFITSIO AND USE_SYSTEM_LIBS)
#	find_package(CFITSIO 3.470)
#endif()

#if(NOT CFITSIO_FOUND)
if(COMPILE_CFITSIO)
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting CFITSIO from source
	#------------------------------------------------------------------------------

	# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP
  #if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
	#	cmake_policy(SET CMP0135 NEW)
	#endif()

	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${CFITSIO_SOURCE_DIR}/CMakeLists.txt")
		message(STATUS "No CFITSIO sources were found; thus, will download it from source:\n${cfitsio_url}")
		ExternalProject_Add(
			cfitsio_src
			URL								"${cfitsio_url}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${CFITSIO_SOURCE_DIR}"
			#BINARY_DIR				"${CFITSIO_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing CFITSIO sources inside:\n${CFITSIO_SOURCE_DIR}")
		add_custom_target(cfitsio_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Creating CMake configure command for CFitsIO
	#------------------------------------------------------------------------------
	# List of arguments to apply to CFitsIO build
	list(APPEND _CFITSIO_ARGS_ 
			# Build type
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			# Specifying compilers
			-DCMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${MPI_C_COMPILER}
			# Specifying the location of ZLIB library (required from version 4.0.0)
			-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
			-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}
		)
	if(CFITSIO_USE_CURL)
		list(APPEND _CFITSIO_ARGS_ 
			-DUSE_CURL:BOOL=ON	
			# Specyfying location of cURL library
			-DCURL_INCLUDE_DIR:PATH=${CURL_INCLUDE_DIR}
			-DCURL_LIBRARY:FILEPATH=${CURL_LIBRARIES}
			)
	else()
		list(APPEND _CFITSIO_ARGS_ 
			-DUSE_CURL:BOOL=OFF	
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and Installing Static and Shared CFitsIO
	#------------------------------------------------------------------------------
	list(APPEND _CFITSIO_LIB_TYPES_ static shared)
	list(APPEND _CFITSIO_LIB_BOOL_VALS_ -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON)
	foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _CFITSIO_LIB_TYPES_ _CFITSIO_LIB_BOOL_VALS_)
		ExternalProject_Add(
			cfitsio_${_lib_type_}
			# Specifying that cfitsio depends on the curl project and should be built after it
			DEPENDS						required_libraries
												zlib
												curl
												cfitsio_src
			PREFIX						"${LIBS_BUILD_DIR}"
			SOURCE_DIR				"${CFITSIO_SOURCE_DIR}"
			INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_CONFIGURE			ON
			LOG_BUILD					ON
			LOG_INSTALL				ON
			# commands how to build the project
			DOWNLOAD_COMMAND	""
			CMAKE_ARGS
				${_CFITSIO_ARGS_}
				# Building both static and shared libraries
				${_bool_val_}
			)
	endforeach()
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(cfitsio
		ALL ""
		DEPENDS cfitsio_static
						cfitsio_shared
		)
	#------------------------------------------------------------------------------
	# setting cfitsio library and include variables
	# Note: for some reason, only set() works but not list(APPEND)
	# Note2: After switching to 4.0.0 version it turns out you need
	# to add Zlib as well, because it is now one of the CFitsIO dependencies
	if(CFITSIO_USE_CURL)
		set(CFITSIO_LIBRARIES 
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}cfitsio${CMAKE_STATIC_LIBRARY_SUFFIX}"
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}cfitsio${CMAKE_SHARED_LIBRARY_SUFFIX}"
			"${CURL_LIBRARIES}"
			"${ZLIB_LIBRARIES}"
			)
	else()
		set(CFITSIO_LIBRARIES 
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}cfitsio${CMAKE_STATIC_LIBRARY_SUFFIX}"
			#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}cfitsio${CMAKE_SHARED_LIBRARY_SUFFIX}"
			"${ZLIB_LIBRARIES}"
			)
	endif()
	set(CFITSIO_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include/"
		)
	include_directories(${CFITSIO_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "CFITSIO LIBRARIES will be: ${CFITSIO_LIBRARIES}")
	#message(STATUS "CFITSIO INCLUDE DIRS will be: ${CFITSIO_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	# adding empty targets in case CFITSIO was found on the system
	add_custom_target(cfitsio ALL "")
	include_directories(${CFITSIO_INCLUDE_DIRS})
	if(CFITSIO_USE_CURL)
		set(CFITSIO_LIBRARIES
			"${CFITSIO_LIBRARY}"
			"${CURL_LIBRARIES}"
			"${ZLIB_LIBRARIES}"
			)
	else()
		set(CFITSIO_LIBRARIES
			"${CFITSIO_LIBRARY}"
			"${ZLIB_LIBRARIES}"
			)
	endif()
	#------------------------------------------------------------------------------
	#message(STATUS "CFITSIO LIBRARIES are: ${CFITSIO_LIBRARIES}")
	#message(STATUS "CFITSIO INCLUDE DIRS are: ${CFITSIO_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
