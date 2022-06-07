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
# Description: This file determines the location of all `required` (the ones which 
# will not be compiled by CMake script) and all `dependent` (the ones which will be
# compiled via the CMake script) libraries. Below is the full list of both of them
# reproduced here for convenience.
#
# Required:
# - Linux Math Library
# - Git
# - MPI 
# - OpenMP
# 
# Dependent:
# - BLAS/LAPACK: MKL, AOCL or OpenBLAS
# - CFITSIO + (if needed) cURL, MbedTLS, and SSH2 
# - HEALPix + LIBSHARP2
# - HDF5 
# - FFTW
#
# The below code defines boolean variables COMPILE_<LIBNAME> to identify which 
# library should be compiled from source. Since the search order of the libraries 
# matters and is reversed comparing to the compilation order, we first look for 
# libs and then include all the subsequent project files at the bottom of this 
# script (these describe how individual libraries should be build). For instance, 
# HDF5 depends on ZLIB and LIBAEC, so the search needs to be done in the following
# order: first HDF5 and then LIBAEC and ZLIB (if HDF5 wasn't found). However, when 
# compiling, ZLIB goes first, then LIBAEC and only then HDF5. Thus, the search and  
# install logic performed (roughly) as follows:
# - if not HDF5, CFITSIO found => install ZLIB if needed => install HDF5, CFITSIO
# - if not HDF5 found  => install LIBAEC if needed 
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
#message(STATUS "math (m) libraries are: ${LIBM_LIBRARY}")
# Printing out the dl libs, which are also required on some unix systems
#message(STATUS "dl libs are: ${CMAKE_DL_LIBS}")
#------------------------------------------------------------------------------
# Looking for Git.
#------------------------------------------------------------------------------
find_package(Git REQUIRED)
#------------------------------------------------------------------------------
# Looking for MPI with C, CXX and Fortran components. 
#------------------------------------------------------------------------------
#message(STATUS "---------------------------------------------------------------")
find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
#find_package(Threads)
#message(STATUS "MPI Libs Path: ${MPI_Fortran_LIBRARIES}")
#message(STATUS "MPI Include Path: ${MPI_Fortran_INCLUDE_PATH}")
#message(STATUS "MPI Compile Options are: ${MPI_Fortran_COMPILE_OPTIONS}") 
#message(STATUS "MPI Link options are: ${MPI_Fortran_LINK_FLAGS}")
# to avoid cmake errors we create and empty target
add_custom_target(mpi ALL "")
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${MPI_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES}) #Threads::Threads)
#------------------------------------------------------------------------------
# Looking for OpenMP. 
#------------------------------------------------------------------------------
#message(STATUS "---------------------------------------------------------------")
find_package(OpenMP REQUIRED)
add_custom_target(openmp ALL "")
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
#message(STATUS "OPENMP Fortran LIBRARIES are: ${OpenMP_Fortran_LIBRARIES}")
#------------------------------------------------------------------------------
# Creating comm_hdf_mod.f90 with Tempita language. Python is required. 
#------------------------------------------------------------------------------
add_custom_target(tempita ALL "")
set(comm_hdf_mod "${COMMANDER3_SOURCE_DIR}/comm_hdf_mod.f90")
# running python command at configure time
execute_process(
	COMMAND ${TEMPITA_DIR}/tempita_proc.py < $< > $@
	INPUT_FILE ${comm_hdf_mod}.in
	OUTPUT_FILE ${comm_hdf_mod}
	)
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
message(STATUS "COMM3_BACKEND was defined as '${COMM3_BACKEND}'")
if(USE_SYSTEM_LIBS)
	#------------------------------------------------------------------------------
	# Performing search for BLAS and LAPACK
	#------------------------------------------------------------------------------
  # From docs: Note C, CXX or Fortran must be enabled to detect a BLAS/LAPACK 
  # library. C or CXX must be enabled to use Intel Math Kernel Library (MKL).
  # Note: Because native (shipped with Linux) BLAS/LAPACK implementations are not 
  # optimized, we require usage of AOCL, MKL or OpenBLAS. 
	#------------------------------------------------------------------------------
	if(USE_SYSTEM_BLAS)
    if(COMM3_BACKEND MATCHES "any")

      get_cpu_vendor(${CPU_DESCRIPTION} CPU_VENDOR)

      if(CPU_VENDOR MATCHES "Intel")
        message(STATUS "Looking for MKL...")
        # Note: Sometimes this doesn't work, i.e. it cannot detect MKL/OpenBLAS 
        # for some weird reason. In this case it is a good idea to logout and login
        # to refresh terminal. In addition, defining two variables, since it may also
        # affect finding MKL on some systems (sometimes it just ignores BLA_VENDOR & 
        # sometimes it needs $ENV{BLA_VENDOR} for MKL)
        set(BLA_VENDOR
            Intel10_32
            Intel10_64lp 
            Intel10_64lp_seq
            Intel10_64ilp
            Intel10_64ilp_seq
            Intel10_64_dyn
          )
        set($ENV{BLA_VENDOR}
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
          set(COMPILE_OPENBLAS TRUE)
        endif()

      elseif(CPU_VENDOR MATCHES "AMD")
        message(STATUS "Looking for AOCL...")
        set(BLA_VENDOR
            FLAME
          )
        # Finds both BLIS and FLAME
        find_package(BLAS)
        find_package(LAPACK)
        if(NOT (BLAS_FOUND OR LAPACK_FOUND))
          set(COMPILE_FLAME TRUE)
        endif()

      elseif(CPU_VENDOR MATCHES "Unknown")
        message(STATUS "Looking for OpenBLAS...")
        set(BLA_VENDOR
            OpenBLAS
          )
        find_package(BLAS)
        find_package(LAPACK)
        if(NOT (BLAS_FOUND OR LAPACK_FOUND))
          set(COMPILE_OPENBLAS TRUE)
        endif()

      else(CPU_VENDOR MATCHES "") #<= just a check, it should be 'Unknown' in this case
        message(FATAL_ERROR 
          "Something went terribly wrong while identifying CPU for BLAS & FFTW3..."
          )
      endif()
    
    elseif(COMM3_BACKEND MATCHES "aocl")
      message(STATUS "Looking for AOCL...")
      set(BLA_VENDOR
          FLAME
        )
      # Finds both BLIS and FLAME
      find_package(BLAS)
      find_package(LAPACK)
      if(NOT (BLAS_FOUND OR LAPACK_FOUND))
        message(STATUS "AOCL is not found => will compile it from source")
        set(COMPILE_FLAME TRUE)
      endif()

    elseif(COMM3_BACKEND MATCHES "mkl")
      message(STATUS "Looking for MKL...")
      set(BLA_VENDOR
          Intel10_32
          Intel10_64lp 
          Intel10_64lp_seq
          Intel10_64ilp
          Intel10_64ilp_seq
          Intel10_64_dyn
        )
      set($ENV{BLA_VENDOR}
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
        message(STATUS "MKL is not found => looking for OpenBLAS instead")#& FFTW3 instead")
        set(BLA_VENDOR
            OpenBLAS
          )
        find_package(BLAS)
        find_package(LAPACK)
        if(NOT (BLAS_FOUND OR LAPACK_FOUND))
          message(STATUS "OpenBLAS is not found => will compile it from source")
          set(COMPILE_OPENBLAS TRUE)
        endif()
      endif()
    
    elseif(COMM3_BACKEND MATCHES "opensrc")
      message(STATUS "Looking for OpenBLAS...")
      set(BLA_VENDOR
          OpenBLAS
        )
      find_package(BLAS)
      find_package(LAPACK)
      if(NOT (BLAS_FOUND OR LAPACK_FOUND))
        message(STATUS "OpenBLAS is not found => will compile it from source")
        set(COMPILE_OPENBLAS TRUE)
      endif()
    else()
      message(FATAL_ERROR 
        "COMM3_BACKEND was defined as ${COMM3_BACKEND}.\n"
        "Possible values are: aocl, mkl, opensrc, any."
        )
    endif()
  endif()
# Since MKL and/or AOCL will be handled above this chunk will be simpler then the one 
# above (?). However, there is no FindFFTW for AMD FFT, so may need to write it myself,
# or do a workaround.
	#------------------------------------------------------------------------------
	# Performing search for FFTW
	#------------------------------------------------------------------------------
	if(USE_SYSTEM_FFTW)
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
	#------------------------------------------------------------------------------
	# Performing search for HDF5 and its dependencies
	#------------------------------------------------------------------------------
	if(USE_SYSTEM_HDF5)
    #message(STATUS "---------------------------------------------------------------")
		# Using dynamic linking instead of static 
		set(HDF5_USE_STATIC_LIBRARIES FALSE)
		# Using serial build instead of parallel
    set(HDF5_PREFER_PARALLEL FALSE)
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
	#------------------------------------------------------------------------------
	if(USE_SYSTEM_CFITSIO)
    #message(STATUS "---------------------------------------------------------------")
		find_package(CFITSIO 3.470)
		if(NOT CFITSIO_FOUND)
			set(COMPILE_CFITSIO TRUE)
		endif()
	else()
		set(COMPILE_CFITSIO TRUE)
	endif()
	#------------------------------------------------------------------------------
	# Looking for HDF5 and CFITSIO mutual dependency -- ZLIB.
	#------------------------------------------------------------------------------
	if(COMPILE_HDF5 OR COMPILE_CFITSIO)
		if(USE_SYSTEM_ZLIB)
      #message(STATUS "---------------------------------------------------------------")
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
      #message(STATUS "---------------------------------------------------------------")
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
    #message(STATUS "---------------------------------------------------------------")
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
	#------------------------------------------------------------------------------
	if(USE_SYSTEM_HEALPIX)
    #message(STATUS "---------------------------------------------------------------")
		find_package(HEALPIX COMPONENTS SHARP Fortran)
		if(NOT HEALPIX_FOUND)
			set(COMPILE_HEALPIX TRUE)
		endif()
	else()
		set(COMPILE_HEALPIX TRUE)
	endif()
else()
  # Identifying which BLAS to compile
  if(COMM3_BACKEND MATCHES "any")
    get_cpu_vendor(${CPU_DESCRIPTION} CPU_VENDOR)
    if(CPU_VENDOR MATCHES "Intel")
      message(STATUS "Cannot compile MKL from source => will compile OpenBLAS & FFTW3 instead")
      set(COMPILE_OPENBLAS TRUE)
    elseif(CPU_VENDOR MATCHES "AMD")
      set(COMPILE_FLAME TRUE)
    elseif(CPU_VENDOR MATCHES "Unknown")
      set(COMPILE_OPENBLAS TRUE)
    else(CPU_VENDOR MATCHES "") #<= just a check, it should be 'Unknown' in this case
      message(FATAL_ERROR 
        "Something went terribly wrong while identifying CPU for BLAS & FFTW3..."
        )
    endif()
  elseif(COMM3_BACKEND MATCHES "aocl")
    set(COMPILE_FLAME TRUE)
  elseif(COMM3_BACKEND MATCHES "opensrc")
    set(COMPILE_OPENBLAS TRUE)
  elseif(COMM3_BACKEND MATCHES "mkl")
    message(STATUS "Cannot compile MKL from source => will compile OpenBLAS & FFTW3 instead")
    set(COMPILE_OPENBLAS TRUE)
  else()
    message(FATAL_ERROR 
      "COMM3_BACKEND was defined as ${COMM3_BACKEND}.\n"
      "Possible values are: aocl, mkl, opensrc, any."
      )
  endif()

#			set(COMPILE_FFTW TRUE)
  set(COMPILE_HDF5 TRUE)
  set(COMPILE_ZLIB TRUE)
  set(COMPILE_LIBAEC TRUE)
  set(COMPILE_CFITSIO TRUE)
	if(CFITSIO_USE_CURL)
    set(COMPILE_CURL TRUE)
    set(COMPILE_MBEDTLS TRUE)
    set(COMPILE_LIBSSH2 TRUE)
  endif()
  set(COMPILE_HEALPIX TRUE)
endif()
#------------------------------------------------------------------------------
# Including the projects to compile 
unset(projects)
# project names <= order matters!!!!!
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
