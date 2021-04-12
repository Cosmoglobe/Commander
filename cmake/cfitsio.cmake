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
#if(NOT (CFITSIO_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#	find_package(CFITSIO 3.470)
#endif()

# TODO: change cfitsio version to 3.49 and figure out how to link cURL correctly.
# Also, try again to switch to cmake build and try to install healpix on top of it.
message(STATUS "---------------------------------------------------------------")
#if(NOT (CFITSIO_FORCE_COMPILE OR ALL_FORCE_COMPILE))
if(USE_SYSTEM_CFITSIO AND USE_SYSTEM_LIBS)
	find_package(CFITSIO 3.470)
endif()

if(NOT CFITSIO_FOUND)
	# Creating configure command for CFITSIO
	list(APPEND cfitsio_configure_command
		"${CMAKE_COMMAND}" "-E" "env" 
		)
	# If user wanted to enable cURL in CFITSIO
	# Adding cURL and other installed libraries to PATH variable
	# because cfitsio is looking for curl-config during configuration.
	if(CFITSIO_USE_CURL)
		set(ENV{PATH} 
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}:$ENV{PATH}
			)
		list(APPEND cfitsio_configure_command
			"PATH=$ENV{PATH}"
			)	
	endif()
	# Adding the rest of commands to the compilation
	list(APPEND cfitsio_configure_command 
		#"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>" 
		)
	if(CFITSIO_USE_CURL)
		list(APPEND cfitsio_configure_command 
			"--enable-curl"
			)
	else()
		list(APPEND cfitsio_configure_command 
			"--disable-curl"
			)
	endif()
	#message("${cfitsio_configure_command}")
	#------------------------------------------------------------------------------
	# Getting CFITSIO from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(cfitsio
		# specifying that cfitsio depends on the curl project and should be built after it
		# Note: I have removed curl support from cfitsio, but curl is left here for now
		DEPENDS required_libraries
						curl
		URL "${cfitsio_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/cfitsio"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/cfitsio/src/cfitsio"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# commands how to build the project
		CONFIGURE_COMMAND "${cfitsio_configure_command}"
		)
	#------------------------------------------------------------------------------
	# setting cfitsio library and include variables
	# Note: for some reason, only set() works but not list(APPEND)
	if(CFITSIO_USE_CURL)
		set(CFITSIO_LIBRARIES 
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX}"
			"${CURL_LIBRARIES}"
			)
	else()
		set(CFITSIO_LIBRARIES 
			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX}"
			)
	endif()
	set(CFITSIO_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include/"
		)
	include_directories(${CFITSIO_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	message(STATUS "CFITSIO LIBRARIES will be: ${CFITSIO_LIBRARIES}")
	message(STATUS "CFITSIO INCLUDE DIRS will be: ${CFITSIO_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	# adding empty targets in case CFITSIO was found on the system
	add_custom_target(cfitsio ALL "")
	include_directories(${CFITSIO_INCLUDE_DIRS})
	if(CFITSIO_USE_CURL)
		set(CFITSIO_LIBRARIES
			"${CFITSIO_LIBRARY}"
			"${CURL_LIBRARIES}"
			)
	else()
		set(CFITSIO_LIBRARIES
			"${CFITSIO_LIBRARY}"
			)
	endif()
	#------------------------------------------------------------------------------
	message(STATUS "CFITSIO LIBRARIES are: ${CFITSIO_LIBRARIES}")
	message(STATUS "CFITSIO INCLUDE DIRS are: ${CFITSIO_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
