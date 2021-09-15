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
# Description: This script determines the location of FFTW on the host system.
# If it fails to do so, it will download, compile and install FFTW from source.
#================================================================================

message(STATUS "---------------------------------------------------------------")
if(USE_SYSTEM_FFTW AND USE_SYSTEM_LIBS)
	find_package(FFTW 
		COMPONENTS 
		DOUBLE 
		DOUBLE_THREADS
		FLOAT 
		#FLOAT_MPI 
		FLOAT_OPENMP
		FLOAT_THREADS
		)
endif()
# Is TRUE if one of the components were missing
if(NOT FFTW_FOUND)
	# Configure command to compile FFTW from source
	# DOUBLE (default one)
	list(APPEND fftw_double_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
		"${FFTW_SOURCE_DIR}/configure" 
		"--prefix=<INSTALL_DIR>")
	# FLOAT
	list(APPEND fftw_float_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
		"${FFTW_SOURCE_DIR}/configure" 
		"--prefix=<INSTALL_DIR>")
	# First, we determine which component were found, so we can link them
	# others will be compiled from source
	# double component is default one, so no configuration option required
	if(NOT FFTW_DOUBLE_FOUND)
		message(STATUS "Missing component - DOUBLE - will be compiled from source")	
	else()
		message(STATUS "Found FFTW_DOUBLE_LIB: ${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_DOUBLE_THREADS_FOUND)
		message(STATUS "Missing component - DOUBLE_THREADS - will be compiled from source")	
		list(APPEND fftw_double_configure_command "--enable-threads")
	else()
		message(STATUS "Found FFTW_DOUBLE_THREADS_LIB: ${FFTW_DOUBLE_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		message(STATUS "Missing component - FLOAT - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-float")
	else()
		message(STATUS "Found FFTW_FLOAT_LIB: ${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		message(STATUS "Missing component - FLOAT_THREADS - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-threads")
	else()
		message(STATUS "Found FFTW_FLOAT_THREADS_LIB: ${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		message(STATUS "Missing component - FLOAT_OPENMP - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-openmp")
	else()
		message(STATUS "Found FFTW_FLOAT_OPENMP_LIB: ${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		message(STATUS "Missing component - FLOAT_MPI - will be compiled from source")	
		list(APPEND fftw_float_configure_command "--enable-mpi")
	else()
		message(STATUS "Found FFTW_FLOAT_MPI_LIB: ${FFTW_FLOAT_MPI_LIB}")
	endif()
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting FFTW from source
	#------------------------------------------------------------------------------
	# Splitting external project add into 3 steps:
	# 1. To download the project
	# 2. To compile with single 
	# 3. and double precision - requires by GNU compilers
	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${FFTW_SOURCE_DIR}/CMakeLists.txt")
		message(STATUS "No FFTW sources were found; thus, will download it from source:\n${fftw_url}")
		ExternalProject_Add(
			fftw_src
			URL								"${fftw_url}"
			URL_MD5						"${fftw_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${FFTW_SOURCE_DIR}"
			BINARY_DIR				"${FFTW_SOURCE_DIR}" 
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# Ommiting Configuration, build and install steps
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing FFTW sources inside:\n${FFTW_SOURCE_DIR}")
		add_custom_target(fftw_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# FFTW Float precision
	#------------------------------------------------------------------------------
	ExternalProject_Add(
		fftw_float
		DEPENDS						fftw_src
		PREFIX						"${LIBS_BUILD_DIR}"
		SOURCE_DIR				"${FFTW_SOURCE_DIR}"
		BINARY_DIR				"${FFTW_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_CONFIGURE			ON 
		LOG_BUILD					ON 
		LOG_INSTALL				ON
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${fftw_float_configure_command}"
		)
	#------------------------------------------------------------------------------
	# FFTW Double precision
	#------------------------------------------------------------------------------
	ExternalProject_Add(
		fftw_double
		DEPENDS						fftw_src
											fftw_float
		PREFIX						"${LIBS_BUILD_DIR}"
		SOURCE_DIR				"${FFTW_SOURCE_DIR}"
		BINARY_DIR				"${FFTW_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_CONFIGURE			ON 
		LOG_BUILD					ON 
		LOG_INSTALL				ON
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${fftw_double_configure_command}"
		)
	# Adding fftw3, fftw3_threads, fftw3_mpi and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	if(NOT FFTW_DOUBLE_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_DOUBLE_LIB}")
	endif()
	if(NOT FFTW_DOUBLE_THREADS_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3_threads${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_LIB}")
	endif()
	if(NOT FFTW_FLOAT_THREADS_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_THREADS_LIB}")
	endif()
	if(NOT FFTW_FLOAT_OPENMP_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_OPENMP_LIB}")
	endif()
	if(NOT FFTW_FLOAT_MPI_FOUND)
		list(APPEND FFTW3_LIBRARIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3f_mpi${CMAKE_STATIC_LIBRARY_SUFFIX}")
	else()
		list(APPEND FFTW3_LIBRARIES "${FFTW_FLOAT_MPI_LIB}")
	endif()
	#------------------------------------------------------------------------------
	add_custom_target(fftw 
		ALL ""
		DEPENDS fftw_float
						fftw_double
		)
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(fftw ALL "")
	#add_custom_target(fftw_double ALL "")
	#add_custom_target(fftw_float ALL "")
	set(FFTW3_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_DOUBLE_THREADS_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
endif()
