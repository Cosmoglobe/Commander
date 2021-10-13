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
# Description: This script determines the location of HEALPix on the host system.
# If it fails to do so, it will download, compile and install HEALPix from source.
#================================================================================

#message(STATUS "---------------------------------------------------------------")
#if(USE_SYSTEM_HEALPIX AND USE_SYSTEM_LIBS)
#	#find_package(HEALPIX 3.70 COMPONENTS SHARP Fortran)
#	find_package(HEALPIX COMPONENTS SHARP Fortran)
#endif()

if(COMPILE_HEALPIX)
	# Writing this to be consistent with fftw.cmake, otherwise 
	# the if statement is unnecessary.
	if(NOT HEALPIX_Fortran_FOUND)
		message(STATUS "Missing component - Fortran - will be compiled from source")	
	endif()
	#------------------------------------------------------------------------------
	# Below flags used to configure Libsharp as part of HEALPix
	if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
		#set(healpix_sharp2_C_FLAGS "-static-intel -O3 -ffast-math -march=native -std=c99 -DUSE_MPI -qopenmp -D__PURE_INTEL_C99_HEADERS__")
		#set(healpix_sharp2_C_FLAGS "-static-intel -O3 -ffast-math -std=c99 -DUSE_MPI -qopenmp -D__PURE_INTEL_C99_HEADERS__")
		#set(healpix_sharp2_C_FLAGS "-static-intel -O3 -ffast-math -mavx2 -std=c99 -DUSE_MPI -qopenmp -D__PURE_INTEL_C99_HEADERS__")
		#MPI support, OpenMP, portable binary:
		set(healpix_sharp2_C_FLAGS "-DUSE_MPI -DMULTIARCH -std=c99 -O3 -ffast-math")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
		#set(healpix_sharp2_C_FLAGS "-O3 -ffast-math -march=native -std=c99 -DUSE_MPI -fopenmp")
		set(healpix_sharp2_C_FLAGS "-DUSE_MPI -DMULTIARCH -std=c99 -O3 -ffast-math")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
		set(healpix_sharp2_C_FLAGS "-O4 -fast -Mipa=fast,inline -Msmartalloc -std=c99 -DUSE_MPI -mp")
	#elseif(CMAKE_Fortran_COMPILER_ID MATCHES NVIDIA)
		#set(healpix_sharp2_C_FLAGS "-O4 -fast -Mipa=fast,inline -Msmartalloc -std=c99 -DUSE_MPI -mp")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
		#set(healpix_sharp2_C_FLAGS "-O4 -fast -Mipa=fast,inline -Msmartalloc -std=c99 -DUSE_MPI -mp")
		set(healpix_sharp2_C_FLAGS "-DUSE_MPI -DMULTIARCH -std=c99 -O3 -ffast-math")
	endif()
	#------------------------------------------------------------------------------
	# Copying modyfied configure script to healpix root
	#list(APPEND healpix_copy_configure_script 
	#	"${CMAKE_COMMAND}" "-E" "copy"
	#	"${CMAKE_SOURCE_DIR}/cmake/third_party/healpix/hpxconfig_functions.sh"
	#	#"${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}/hpxconfig_functions.sh" 
	#	"${HEALPIX_SOURCE_DIR}/hpxconfig_functions.sh" 
	#	#"&&"
	#	)
	# Creating configure command for HEALPix
	list(APPEND healpix_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		)
	if(CFITSIO_USE_CURL)
		set(ENV{PATH} 
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}:$ENV{PATH}
			)
		set(ENV{LD_LIBRARY_PATH} 
			${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$ENV{LD_LIBRARY_PATH}
			)
		list(APPEND healpix_configure_command
			"PATH=$ENV{PATH}"
			"LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}"
			)	
	endif()
	# TODO: make this work with and without CFITSIO_ROOT
	# This method is not optimal, need to figure somethin else.
	set(CFITSIO_ROOT
		$ENV{CFITSIO_ROOT}
		)
	if(CFITSIO_ROOT AND NOT(HEALPIX_FORCE_COMPILE OR ALL_FORCE_COMPILE))
		list(APPEND healpix_configure_command 
			"FITSDIR=$ENV{CFITSIO_ROOT}/lib"#${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
			)
	else()
		list(APPEND healpix_configure_command 
			"FITSDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
			)
	endif()
	list(APPEND healpix_configure_command 
		#"${CMAKE_COMMAND}" "-E" "env" 
		#"PATH=$ENV{PATH}"
		#"FITSDIR=$ENV{CFITSIO_ROOT}/lib"#${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
		#"FITSDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
		#"FITSINC=${CMAKE_INSTALL_PREFIX}/include"
		#"FITSDIR=${CFITSIO_LIBRARY}"
		"FITSINC=${CFITSIO_INCLUDE_DIRS}"
		"F_SHARED=0"
		#"F_SHARED=1"
		"FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"SHARP_COPT=${healpix_sharp2_C_FLAGS}"
		"./configure" 
		"--auto=sharp,f90" #${healpix_components}" #profile,f90,c,cxx;" 
		#"--prefix=<INSTALL_DIR>" 
		)
	#------------------------------------------------------------------------------
	# Getting HEALPix from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${HEALPIX_SOURCE_DIR}/configure")
		message(STATUS "No HEALPIX sources were found; thus, will download it from source:\n${healpix_url}")
		ExternalProject_Add(
			healpix_src
			URL								"${healpix_url}"
			URL_MD5						"${healpix_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${HEALPIX_SOURCE_DIR}"
			BINARY_DIR				"${HEALPIX_SOURCE_DIR}" 
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing HEALPIX sources inside:\n${HEALPIX_SOURCE_DIR}")
		add_custom_target(healpix_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and installing HEALPix
	#------------------------------------------------------------------------------
	# Note: HEALPix doesn't have an install command
	ExternalProject_Add(
		healpix
		DEPENDS						required_libraries
											curl
											cfitsio 
											healpix_src
		PREFIX						"${LIBS_BUILD_DIR}"
		SOURCE_DIR				"${HEALPIX_SOURCE_DIR}"
		BINARY_DIR				"${HEALPIX_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_CONFIGURE			ON
		LOG_BUILD					ON
		# commands how to build the project
		DOWNLOAD_COMMAND	""
		#CONFIGURE_COMMAND "${healpix_copy_configure_script}"
		#COMMAND						"${healpix_configure_command}"
		CONFIGURE_COMMAND	"${healpix_configure_command}"
		# HEALPix doesn't have an install command 
		INSTALL_COMMAND		""
		# copying Healpix and all its files (src and compiled) into CMAKE_INSTALL_PREFIX directory
		COMMAND ${CMAKE_COMMAND} -E copy_directory "${HEALPIX_SOURCE_DIR}" "${HEALPIX_INSTALL_PREFIX}"
		)
	#------------------------------------------------------------------------------
	set(HEALPIX_LIBRARIES 
		${HEALPIX_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sharp${CMAKE_STATIC_LIBRARY_SUFFIX}
		${HEALPIX_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}healpix${CMAKE_STATIC_LIBRARY_SUFFIX}
		)
	set(HEALPIX_INCLUDE_DIRS
		${HEALPIX_INSTALL_PREFIX}/include
		${HEALPIX_INSTALL_PREFIX}/include/libsharp
		)
	#include_directories("${CMAKE_INSTALL_PREFIX}/healpix/include")
	#include_directories("${HEALPIX_INSTALL_PREFIX}/include")
	#include_directories("${HEALPIX_INSTALL_PREFIX}/include/libsharp")
	include_directories("${HEALPIX_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
	message(STATUS "HEALPIX LIBRARIES will be: ${HEALPIX_LIBRARIES}")
	message(STATUS "HEALPix INCLUDES will be: ${HEALPIX_INCLUDE_DIRS}")
else()
	add_custom_target(healpix ALL "")
	message(STATUS "HEALPix LIBRARIES are: ${HEALPIX_LIBRARIES}")
	message(STATUS "HEALPix INCLUDES are: ${HEALPIX_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
	include_directories("${HEALPIX_INCLUDE_DIRS}")
endif()
