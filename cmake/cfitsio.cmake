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
	# Creating configure command for CFITSIO
	set(cfitsio_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>" 
		"--enable-curl=${CURL_LIBRARIES}"
		"--disable-curl"
		)
	#------------------------------------------------------------------------------
	# Getting CFITSIO from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}
		# specifying that cfitsio depends on the curl project and should be built after it
		# Note: I have removed curl support from cfitsio, but curl is left here for now
		DEPENDS curl
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
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		#COMMAND ${CMAKE_COMMAND} -E env --unset=PATH PATH=$ENV{PATH} ./configure --prefix=<INSTALL_DIR> 
		#COMMAND export PATH=${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl-config #"${${project}_configure_command}"
		#COMMAND export PATH=$PATH:/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild2/build/install/bin #"${${project}_configure_command}"
		#COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} MPFC=${CMAKE_Fortran_COMPILER} MPCXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR> --disable-curl # <= if specified manually, the cmake args will not work
		#COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --disable-curl # <= if specified manually, the cmake args will not work
		#CMAKE_ARGS
		# specifying where to find curl library
		#-DCURL_INCLUDE_DIR:PATH=${CURL_INCLUDE_DIR}
		#-DCURL_BINARY_DIR:PATH=${CURL_BINARY_DIR}
		#-DCURL_SOURCE_DIR:PATH=${CURL_SOURCE_DIR}
		#-DCURL_LIBRARIES:PATH=${CURL_LIBRARIES}
		# specifying where to install CFitsio
		# (and how to install it)
		#-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		#-DBUILD_SHARED_LIBS:BOOL=OFF
		# specifying build command
		#BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Debug #--target INSTALL
		#BUILD_IN_SOURCE 1	
		)

	# setting cfitsio library and include variables
	set(CFITSIO_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(CFITSIO_INCLUDES "${CMAKE_INSTALL_PREFIX}/include/")
	include_directories(${CFITSIO_INCLUDES})

	message(STATUS "CFITSIO LIBRARIES will be: ${CFITSIO_LIBRARIES}")
	message(STATUS "CFITSIO INCLUDE DIR will be: ${CFITSIO_INCLUDES}")
else()
	# adding empty targets in case CFITSIO was found on the system
	add_custom_target(${project} ALL "")
	include_directories(${CFITSIO_INCLUDES})
endif()
