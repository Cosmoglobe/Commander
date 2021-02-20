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
# Description: This script determines the location of CFitsio on the host system.
# If it fails to do so, it will download, compile and install CFitsio from source.
# Together with cURL, CFitsio is required to successfully compile HEALPix.
#================================================================================

message(STATUS "---------------------------------------------------------------")
if(NOT (CFITSIO_FORCE_COMPILE OR ALL_FORCE_COMPILE))
	find_package(CFITSIO 3.47)
endif()

if(NOT CFITSIO_FOUND)
	#------------------------------------------------------------------------------
	# Getting CFITSIO from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}
		# specifying that cfitsio depends on the curl project and should be built after it
		# Note: I have removed curl support from cfitsio, but curl is left here for now
		DEPENDS required_libraries curl
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
			-DBUILD_SHARED_LIBS:BOOL=OFF # there are problems when building shared library
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${COMMANDER3_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${COMMANDER3_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${COMMANDER3_C_COMPILER}
			# Specyfying location of cURL library
			-DCURL_INCLUDE_DIR:PATH=${CURL_INCLUDE_DIR}
			#-DCURL_LIBRARY_RELEASE:FILEPATH=${CURL_LIBRARIES}
			-DCURL_LIBRARY:FILEPATH=${CURL_LIBRARIES}
			# Specyfing the location of ZLIB library
			-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
			-DZLIB_LIBRARY_RELEASE:FILEPATH=${ZLIB_LIBRARIES}
		)

	# setting cfitsio library and include variables
	set(CFITSIO_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(CFITSIO_INCLUDES "${CMAKE_INSTALL_PREFIX}/include/")
	include_directories(${CFITSIO_INCLUDES})

	message(STATUS "CFITSIO LIBRARIES will be: ${CFITSIO_LIBRARIES}")
	message(STATUS "CFITSIO INCLUDE DIR will be: ${CFITSIO_INCLUDES}")
else()
	# adding empty targets in case CFITSIO was found on the system
	add_custom_target(${project} 
		ALL ""
		DEPENDS required_libraries curl
		)
	include_directories(${CFITSIO_INCLUDES})
endif()
