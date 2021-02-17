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
# Description: This script determines the location of FFTW on the host system.
# If it fails to do so, it will download, compile and install FFTW from source.
#================================================================================

message(STATUS "---------------------------------------------------------------")
# TODO: make it so components will matter because now it install everything because 
# I gave the command to add appropriate configure suboptions to configure command
#if(NOT (FFTW_FORCE_COMPILE OR ALL_FORCE_COMPILE))
if(USE_SYSTEM_FFTW AND USE_SYSTEM_LIBS)
	find_package(FFTW 
		COMPONENTS 
		DOUBLE 
		#DOUBLE_THREADS
		FLOAT 
		#FLOAT_MPI 
		FLOAT_OPENMP
		#FLOAT_THREADS
		)
endif()
# Is TRUE if one of the components were missing
if(NOT FFTW_FOUND)
	# First, we determine which component were found, so we can link them
	# others will be compiled from source
	# double component is default one, so no configuration option required
	if(NOT FFTW_DOUBLE_FOUND)
		message(STATUS "Missing component - DOUBLE - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_DOUBLE_LIB: ${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_DOUBLE_OPENMP_FOUND)
		message(STATUS "Missing component - DOUBLE_OPENMP - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_DOUBLE_OPENMP_LIB: ${FFTW_DOUBLE_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		message(STATUS "Missing component - FLOAT - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_FLOAT_LIB: ${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		message(STATUS "Missing component - FLOAT_OPENMP - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_FLOAT_OPENMP_LIB: ${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	#------------------------------------------------------------------------------
	# Getting FFTW from source
	#------------------------------------------------------------------------------
	# Splitting external project add into 3 steps:
	# 1. To download the project
	# 2. To compile with single and double precision - requires by GNU compilers
	ExternalProject_Add(fftw_src
		DEPENDS required_libraries
		URL "${fftw_url}"
		URL_MD5 "${fftw_md5}"
		# PREFIX should be present, otherwise it will pull it into "build" dir
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		# Ommiting Configuration, build and install steps
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		)

	#------------------------------------------------------------------------------
	# FFTW Float precision
	#------------------------------------------------------------------------------
	# Configure scripts usually compile both static and shared libs, but in this case
	# we will have only of them, so we need to call compilations script twice -- for
	# once for static libs and once for shared ones.	
	ExternalProject_Add(fftw_float_static
		DEPENDS fftw_src 
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_INSTALL_LIBDIR=lib
			# FFTW configuration
			-DBUILD_SHARED_LIBS:BOOL=OFF
			-DENABLE_FLOAT:BOOL=ON
			-DENABLE_OPENMP:BOOL=ON
			-DENABLE_AVX:BOOL=ON
			-DENABLE_AVX2:BOOL=ON
			-DENABLE_SSE2:BOOL=ON
			-DENABLE_SSE:BOOL=ON
			-DBUILD_TESTS:BOOL=OFF
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		)
	
	ExternalProject_Add(fftw_float_shared
		DEPENDS fftw_src 
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON 
		LOG_BUILD ON 
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_INSTALL_LIBDIR=lib
			# FFTW configuration
			-DBUILD_SHARED_LIBS:BOOL=ON
			-DENABLE_FLOAT:BOOL=ON
			-DENABLE_OPENMP:BOOL=ON
			-DENABLE_AVX:BOOL=ON
			-DENABLE_AVX2:BOOL=ON
			-DENABLE_SSE2:BOOL=ON
			-DENABLE_SSE:BOOL=ON
			-DBUILD_TESTS:BOOL=OFF
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		)
	add_custom_target(fftw_float 
		ALL ""
		DEPENDS fftw_float_static
						fftw_float_shared
			)
	#------------------------------------------------------------------------------
	# FFTW Double precision
	#------------------------------------------------------------------------------
	ExternalProject_Add(fftw_double_static
		DEPENDS fftw_src 
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_INSTALL_LIBDIR=lib
			# FFTW configuration
			-DBUILD_SHARED_LIBS:BOOL=OFF
			-DENABLE_OPENMP:BOOL=ON
			-DENABLE_AVX:BOOL=ON
			-DENABLE_AVX2:BOOL=ON
			-DENABLE_SSE2:BOOL=ON
			-DENABLE_SSE:BOOL=ON
			-DBUILD_TESTS:BOOL=OFF
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		)

	ExternalProject_Add(fftw_double_shared
		DEPENDS fftw_src 
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/fftw"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/fftw/src/fftw_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		BUILD_ALWAYS FALSE
		# Commands to configure, build and install the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_INSTALL_LIBDIR=lib
			# FFTW configuration
			-DBUILD_SHARED_LIBS:BOOL=ON
			-DENABLE_OPENMP:BOOL=ON
			-DENABLE_AVX:BOOL=ON
			-DENABLE_AVX2:BOOL=ON
			-DENABLE_SSE2:BOOL=ON
			-DENABLE_SSE:BOOL=ON
			-DBUILD_TESTS:BOOL=OFF
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		)

	add_custom_target(fftw_double
		ALL ""
		DEPENDS fftw_double_static
						fftw_double_shared
			)
	#------------------------------------------------------------------------------
	# adding fftw3, fftw3_omp, and fftw3f, fftws3f_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	set(FFTW3_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3_omp${CMAKE_SHARED_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f_omp${CMAKE_SHARED_LIBRARY_SUFFIX}"		
		)

	#------------------------------------------------------------------------------
	add_custom_target(fftw
		ALL ""
		DEPENDS fftw_double
						fftw_float
			)
	#------------------------------------------------------------------------------
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(fftw ALL "")
	add_custom_target(fftw_double ALL "")
	add_custom_target(fftw_float ALL "")
	set(FFTW3_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_DOUBLE_OPENMP_LIB}
		#${FFTW_DOUBLE_THREADS_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		#${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
endif()
