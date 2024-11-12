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
# Description: This file determines the location of all required (the ones which 
# will not be compiled by CMake script) and all dependent (the ones which will be
# compiled via the CMake script) libraries. Below is the full list of both of them
# reproduced here for convenience.
# Required:
# - Git
# - MPI 
# - OpenMP
# Dependent:
# - BLAS/LAPACK
# - CFITSIO
# - HEALPix + LIBSHARP2
# - HDF5 
# - FFTW
# The search for the latter is performed depending on user preferences. If lib-
# rary was not found, the special variable COMPILE_<LIBNAME> will be set to TRUE
# and thus this library will be compiled from source. This is done to avoid com-
# pilation of some libraries on which Commander3 doesn't depend on but its depen-
# dencies do (e.g. ZLIB and LIBAEC which were compiled all the time before).
# The search logic is performed (roughly) as follows:
# - if not HDF5, CFITSIO => install ZLIB
# - HDF5 => LIBAEC
# - if not CFITSIO => install CURL => install ZLIB, MbedTLS, LIBSSH2
# - HEALPix
# - BLAS/LAPACK;
# - FFTW;
#------------------------------------------------------------------------------
# Required Libraries
#------------------------------------------------------------------------------
# Looking for Linux Math Library
#------------------------------------------------------------------------------
find_library(LIBM_LIBRARY m)
message(STATUS "math (m) libraries are: ${LIBM_LIBRARY}")
# Printing out the dl libs, which are also required on some unix systems
message(STATUS "dl libs are: ${CMAKE_DL_LIBS}")
#------------------------------------------------------------------------------
# Looking for Git.
#------------------------------------------------------------------------------
find_package(Git REQUIRED)
#------------------------------------------------------------------------------
# Looking for MPI with C, CXX and Fortran components. 
#------------------------------------------------------------------------------
message(STATUS "---------------------------------------------------------------")
find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
#find_package(Threads)
message(STATUS "MPI Libs Path: ${MPI_Fortran_LIBRARIES}")
message(STATUS "MPI Include Path: ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "MPI Compile Options are: ${MPI_Fortran_COMPILE_OPTIONS}") 
message(STATUS "MPI Link options are: ${MPI_Fortran_LINK_FLAGS}")
# to avoid cmake errors we create and empty target
add_custom_target(mpi ALL "")
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${MPI_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES}) #Threads::Threads)
#------------------------------------------------------------------------------
# Looking for OpenMP. 
#------------------------------------------------------------------------------
message(STATUS "---------------------------------------------------------------")
find_package(OpenMP REQUIRED)
add_custom_target(openmp ALL "")
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
message(STATUS "OPENMP Fortran LIBRARIES are: ${OpenMP_Fortran_LIBRARIES}")
#------------------------------------------------------------------------------
# Creating comm_hdf_mod.f90 with Tempita language. Python is required. 
#------------------------------------------------------------------------------
add_custom_target(tempita ALL "")
set(comm_hdf_mod "${COMMANDER3_SOURCE_DIR}/comm_hdf_mod.f90")
#------------------------------------------------------------------------------
# Creating a Unified target
#------------------------------------------------------------------------------
add_custom_target(required_libraries ALL "" 
	DEPENDS tempita 
					mpi
					openmp
					)
#------------------------------------------------------------------------------
# Dependent/Compiled libraries
#------------------------------------------------------------------------------
if(USE_SYSTEM_LIBS)
	#------------------------------------------------------------------------------
	# Performing search for BLAS and LAPACK
	if(USE_SYSTEM_BLAS)
		# TODO: Need to add FindMKL.cmake to separately search for MKL, which
		# has both BLAS/LAPACK and FFTW3. If it exists will go with that, and
		# if not then seacrh for OpenBLAS and FFTW3 and compile those (if necessary)
		#
		# Note: Sometimes this doesn't work, i.e. it cannot detect MKL/OpenBLAS 
		# for some weird reason. In this case it is a good idea to logout and login
		# to refresh terminal.
		set($ENV{BLA_VENDOR} 
				OpenBLAS
				Intel10_32
				Intel10_64lp
				Intel10_64lp_seq
				Intel10_64ilp
				Intel10_64ilp_seq
				Intel10_64_dyn
				)
		find_package(BLAS)
		find_package(LAPACK)
		if(NOT (BLAS_FOUND OR LAPACK_FOUND))
			set(COMPILE_BLAS TRUE)
		endif()
	else()
		set(COMPILE_BLAS TRUE)
	endif()
	#------------------------------------------------------------------------------
	# Performing search for HDF5 and its dependencies
	if(USE_SYSTEM_HDF5)
		message(STATUS "---------------------------------------------------------------")
		# Using static linking instead of dynamic
		set(HDF5_USE_STATIC_LIBRARIES FALSE)
		# Using parallel build instead of serial
		set(HDF5_PREFER_PARALLEL TRUE)
		#find_package(HDF5 1.12.0 COMPONENTS Fortran Fortran_HL)
		find_package(HDF5 1.10.5 COMPONENTS Fortran Fortran_HL)
		if(NOT HDF5_FOUND)
			set(COMPILE_HDF5 TRUE)
		endif()
	else()
		set(COMPILE_HDF5 TRUE)
	endif()
	#------------------------------------------------------------------------------
	# Performing search for CFITSIO
	if(USE_SYSTEM_CFITSIO)
		message(STATUS "---------------------------------------------------------------")
		find_package(CFITSIO 3.470)
		if(NOT CFITSIO_FOUND)
			set(COMPILE_CFITSIO TRUE)
		endif()
	else()
		set(COMPILE_CFITSIO TRUE)
	endif()
	#------------------------------------------------------------------------------
	# Looking for HDF5 and CFITSIO mutual dependency -- ZLIB.
	if(COMPILE_HDF5 OR COMPILE_CFITSIO)
		if(USE_SYSTEM_ZLIB)
			message(STATUS "---------------------------------------------------------------")
			find_package(ZLIB 1.2.11)
			if(NOT ZLIB_FOUND)
				set(COMPILE_ZLIB TRUE)
			endif()
		else()
			set(COMPILE_ZLIB TRUE)
		endif()
	endif()
	# Other HDF5 dependencies
	if(COMPILE_HDF5)	
		if(USE_SYSTEM_LIBAEC)
			message(STATUS "---------------------------------------------------------------")
			# Placeholder for FindLIBAEC.cmake
			message(STATUS "No LibAEC, will compile from source.")
			# find_package(LIBAEC)
			if(NOT LIBAEC_FOUND)
				set(COMPILE_LIBAEC TRUE)
			endif()
		else()
			set(COMPILE_LIBAEC TRUE)
		endif()
	endif()
	# Other CFITSIO dependencies. Compile only if we want to have CURL support
	if(COMPILE_CFITSIO AND CFITSIO_USE_CURL)
		message(STATUS "---------------------------------------------------------------")
		if(USE_SYSTEM_CURL)
			# CMake configure scripts (in versions 7.69-7.74) doesn't work properly,
			# so we look for cURL in a standard manner.
			set(CURL_NO_CURL_CMAKE ON)
			find_package(CURL 7.54)
			if(NOT CURL_FOUND)
				set(COMPILE_CURL TRUE)
			endif()
		else()
			set(COMPILE_CURL TRUE)
		endif()
		# Looking for CURL dependencies
		if(COMPILE_CURL)
			# MbedTLS (depends on ZLIB)
			if(USE_SYSTEM_MBEDTLS)
				find_package(MBEDTLS)
				if(NOT MBEDTLS_FOUND)
					set(COMPILE_MBEDTLS TRUE)
				endif()
			else()
				set(COMPILE_MBEDTLS TRUE)
			endif()
			# LIBSSH2 (depends on ZLIB and MbedTLS)
			if(USE_SYSTEM_LIBSSH2)
				find_package(LibSSH2)
				if(NOT LIBSSH2_FOUND)
					set(COMPILE_LIBSSH2 TRUE)
				endif()
			else()
				set(COMPILE_LIBSSH2 TRUE)
			endif()
		endif()
	endif()
	#------------------------------------------------------------------------------
	# Performing search for HEALPix
	if(USE_SYSTEM_HEALPIX)
		message(STATUS "---------------------------------------------------------------")
		find_package(HEALPIX COMPONENTS SHARP Fortran)
		if(NOT HEALPIX_FOUND)
			set(COMPILE_HEALPIX TRUE)
		endif()
	else()
		set(COMPILE_HEALPIX TRUE)
	endif()
	#------------------------------------------------------------------------------
	# Performing search for FFTW
	if(USE_SYSTEM_FFTW)
		message(STATUS "---------------------------------------------------------------")
		find_package(FFTW
			COMPONENTS
			DOUBLE
			DOUBLE_THREADS
			FLOAT
			FLOAT_OPENMP
			FLOAT_THREADS
			)
		if(NOT FFTW_FOUND)
			set(COMPILE_FFTW TRUE)
		endif()
	else()
		set(COMPILE_FFTW TRUE)
	endif()
endif()
#------------------------------------------------------------------------------
# Including the projects to compile 
unset(projects)
# project names <= order matters
list(APPEND projects 
	zlib
	libaec
	mbedtls
	libssh2
	curl
	cfitsio
	blas # blas-lapack module 
	##sharp2
	fftw
	hdf5
	healpix
	#camb
	commander3
	)
# Include all project configuration files
foreach(project ${projects})
	include("${project}")
endforeach()
