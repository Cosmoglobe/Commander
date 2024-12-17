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

# TODO: I am using FFTW3_* instead of FFTW_* variables (because MPI is treated 
# separetely), that is why in the ned I have FFTW3_LIBRARIES instead of 
# FFTW_LIBRARIES. Perhaps this is incorrect? Need to take a look into FindFFTW.cmake
# and rewrite it if necessary.
#message(STATUS "---------------------------------------------------------------")
#if(USE_SYSTEM_FFTW AND USE_SYSTEM_LIBS)
#	find_package(FFTW 
#		COMPONENTS 
#		DOUBLE 
#		DOUBLE_THREADS
#		FLOAT 
#		#FLOAT_MPI 
#		FLOAT_OPENMP
#		FLOAT_THREADS
#		)
#endif()
# Is TRUE if one of the components were missing
if(COMPILE_FFTW)
	#------------------------------------------------------------------------------
	# Splitting the project into 5 steps:
	# 1. To download the project
	# 2. To compile with double static 
	# 3. To compile with double shared 
	# 4. To compile with float static 
	# 5. To compile with float shared 
	#
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting FFTW from source
	#------------------------------------------------------------------------------

	# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP
  #if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
	#	cmake_policy(SET CMP0135 NEW)
	#endif()

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
	# Compiling and Installing Static and Shared FFTW
	#------------------------------------------------------------------------------
	# Looping over libraries we need to compile
	list(APPEND _FFTW_NAMES_ double float)
	list(APPEND _FFTW_ARGS_ -DENABLE_FLOAT:BOOL=OFF -DENABLE_FLOAT:BOOL=ON)
	list(APPEND _FFTW_LIB_TYPE_ shared static)
	list(APPEND _FFTW_LIB_BOOL_VAL_ -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF)
	foreach(_fftw_component_ _fftw_arg_ IN ZIP_LISTS _FFTW_NAMES_ _FFTW_ARGS_)
		foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _FFTW_LIB_TYPE_ _FFTW_LIB_BOOL_VAL_)
			ExternalProject_Add(
				fftw_${_fftw_component_}_${_lib_type_}
				DEPENDS           fftw_src
				PREFIX            "${LIBS_BUILD_DIR}"
				SOURCE_DIR        "${FFTW_SOURCE_DIR}"
				#BINARY_DIR        "${FFTW_SOURCE_DIR}"
				INSTALL_DIR       "${CMAKE_INSTALL_PREFIX}"
				LOG_DIR           "${CMAKE_LOG_DIR}"
				LOG_CONFIGURE     ON 
				LOG_BUILD         ON 
				LOG_INSTALL       ON 
				# Disabling download
				DOWNLOAD_COMMAND  ""
				CMAKE_ARGS
				  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
					# Specifying installations paths for binaries and libraries
					-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
					# Specifying compilers
					-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
					-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
					-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
					# Building both static and shared libraries
					${_bool_val_}
					# Which libraries to produce
					-DENABLE_OPENMP:BOOL=ON
					-DENABLE_THREADS:BOOL=ON
					#-DENABLE_AVX2:BOOL=${FFTW_ENABLE_AVX2}
					-DENABLE_AVX2:BOOL=OFF
					${_fftw_arg_}
					# ensuring it will be installed inside `lib` and not `lib64`
					-DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
				)
		endforeach()
	endforeach()
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(fftw 
		ALL ""
		DEPENDS fftw_float_static
						fftw_double_static
		        fftw_float_shared
						fftw_double_shared
		)
	#------------------------------------------------------------------------------
	# Adding fftw3, fftw3_threads, and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	set(FFTW_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"	
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		)
	set(FFTW_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include"	
		)
	include_directories(${FFTW_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "FFTW LIBRARIES will be: ${FFTW_LIBRARIES}")
	#message(STATUS "FFTW INCLUDE DIRS will be: ${FFTW_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(fftw ALL "")
	set(FFTW_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_DOUBLE_THREADS_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
	include_directories(${FFTW_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "FFTW LIBRARIES will be: ${FFTW_LIBRARIES}")
	#message(STATUS "FFTW INCLUDE DIRS will be: ${FFTW_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
