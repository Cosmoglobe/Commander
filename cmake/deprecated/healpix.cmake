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
# Description: This script determines the location of HEALPix on the host system.
# If it fails to do so, it will download, compile and install HEALPix from source.
#================================================================================

message(STATUS "---------------------------------------------------------------")
if(NOT (HEALPIX_FORCE_COMPILE OR ALL_FORCE_COMPILE))
	find_package(HEALPIX 3.70 COMPONENTS SHARP Fortran)
endif()

if(NOT HEALPIX_FOUND)
	# Writing this to be consistent with fftw.cmake, otherwise 
	# the if statement is unnecessary.
	if(NOT HEALPIX_Fortran_FOUND)
		message(STATUS "Missing component - Fortran - will be compiled from source")	
	endif()
	#------------------------------------------------------------------------------
	# Below flags used to configure Libsharp as part of HEALPix
	if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
		set(healpix_sharp2_C_FLAGS "-O3 -ffast-math -march=native -std=c99 -DUSE_MPI -qopenmp")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
		set(healpix_sharp2_C_FLAGS "-O3 -ffast-math -march=native -std=c99 -DUSE_MPI -fopenmp")
	elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
		set(healpix_sharp2_C_FLAGS "-O4 -fast -Mipa=fast,inline -Msmartalloc -std=c99 -DUSE_MPI -mp")
	endif()
	#------------------------------------------------------------------------------
	# Creating configure command for HEALPix
	#message("my message ${CFITSIO_ROOT_DIR}")
	#message("my message ${CFITSIO_ROOT}")
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
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"SHARP_COPT=${healpix_sharp2_C_FLAGS}"
		"./configure" 
		"--auto=f90" #${healpix_components}" #profile,f90,c,cxx;" 
		#"--prefix=<INSTALL_DIR>" 
		)
	#------------------------------------------------------------------------------
	# Getting HEALPix from source
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}
		URL "${${project}_url}"
		URL_MD5 "${${project}_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		#SOURCE_DIR "${download_dir}/${project}/src/${project}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" 
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		#INSTALL_DIR "${HEALPIX_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		#BUILD_ALWAYS FALSE
		BUILD_ALWAYS TRUE 
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		#BUILD_COMMAND ""
		# making healpix to be installed the last before commander3
		DEPENDS cfitsio 
						#hdf5 
						#sharp2 
						#fftw 
						#fftw_double 
						#fftw_float 
						#doxygen 
						#tempita 
						#blas 
						#openmp 
						#curl 
						#mpi 
						#zlib
		# HEALPix doesn't have an install command 
		INSTALL_COMMAND ""
		# copying Healpix and all its files (src and compiled) into CMAKE_INSTALL_PREFIX directory
		#COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_INSTALL_PREFIX}/healpix"
		COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${HEALPIX_INSTALL_PREFIX}"
		)

	#ExternalProject_Add_Step(${project} ${project}_copy_step
	#		COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${HEALPIX_INSTALL_PREFIX}"
	#	ALWAYS FALSE
	#		)
	#ExternalProject_Add_StepTargets(${project} ${project}_copy_step)

	#set(HEALPIX_LIBRARIES 
	#	${CMAKE_INSTALL_PREFIX}/healpix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sharp${CMAKE_STATIC_LIBRARY_SUFFIX}
	#	${CMAKE_INSTALL_PREFIX}/healpix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX}
	#	)
	set(HEALPIX_LIBRARIES 
		${HEALPIX_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}sharp${CMAKE_STATIC_LIBRARY_SUFFIX}
		${HEALPIX_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX}
		)
	#include_directories("${CMAKE_INSTALL_PREFIX}/healpix/include")
	include_directories("${HEALPIX_INSTALL_PREFIX}/include")
	#------------------------------------------------------------------------------
	message(STATUS "HEALPIX LIBRARIES will be: ${HEALPIX_LIBRARIES}")
else()
	add_custom_target(${project} ALL "")
	message(STATUS "HEALPIX LIBRARIES are: ${HEALPIX_LIBRARIES}")
	#------------------------------------------------------------------------------
	include_directories("${HEALPIX_INCLUDE_DIRS}")
endif()
