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
# Description: This script determines the location of several libraries/projects 
# required to compile Commander3. The full list of them is:
# - Git
# - MPI
# - OpenMP
# - ZLIB
# - BLAS/LAPACK
# To avoid different cmake errors, we create an empty targets for each of the projects.
# We are also using tempita language to generate comm_hdf_mod.f90 from comm_hdf_mod.f90.in
#================================================================================

message(STATUS "---------------------------------------------------------------")
message(STATUS "Looking for packages...")
# use this to write your own find_package
find_package(PkgConfig)
#------------------------------------------------------------------------------
# Looking for Linux Math Library. 
#------------------------------------------------------------------------------
# finding math library
find_library(LIBM_LIBRARY m)
message(STATUS "math (m) libraries are: ${LIBM_LIBRARY}")
# printing out the dl libs, which are also required on some unix systems
message(STATUS "dl libs are: ${CMAKE_DL_LIBS}")
#------------------------------------------------------------------------------
# Looking for Git.
#------------------------------------------------------------------------------
# We will be using Git to download some dependencies, so we need to check if git available
find_package(Git REQUIRED)


#message(STATUS "---------------------------------------------------------------")
#find_package(OpenSSL REQUIRED)

#message(STATUS "OpenSSL INCLUDE DIR is: ${OPENSSL_INCLUDE_DIR}")
#message(STATUS "${OPENSSL_SSL_LIBRARIES}")
#message(STATUS "${OPENSSL_CRYPTO_LIBRARIES}")
#message(STATUS "OpenSSL LIBRARIES are: ${OPENSSL_LIBRARIES}")

#message(STATUS "---------------------------------------------------------------")
#find_package(LibSSH2)
#find_package(ZLIB REQUIRED)

#message(STATUS "ZLIB LIBRARIES are: ${ZLIB_LIBRARIES}")
#message(STATUS "ZLIB INCLUDE DIRS are: ${ZLIB_INCLUDE_DIRS}")

#add_custom_target(zlib ALL "")
#include_directories(${ZLIB_INCLUDE_DIRS})
#add_library(zlib_lib SHARED IMPORTED GLOBAL) 
#set_target_properties(zlib_lib PROPERTIES IMPORTED_LOCATION ${ZLIB_LIBRARIES})

#------------------------------------------------------------------------------
# Looking for MPI with C, CXX and Fortran components. 
#------------------------------------------------------------------------------
message(STATUS "---------------------------------------------------------------")
find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
find_package(Threads)
# printing out status of the search
message(STATUS "MPI Libs Path: ${MPI_Fortran_LIBRARIES}")
message(STATUS "MPI Include Path: ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "MPI Compile Options are: ${MPI_Fortran_COMPILE_OPTIONS}") 
message(STATUS "MPI Link options are: ${MPI_Fortran_LINK_FLAGS}")
message(STATUS "THREADS Libs Path: ${CMAKE_THREAD_LIBS_INIT}")#MPIexec: ${MPIEXEC_EXECUTABLE}"

# to avoid cmake errors we create and empty target
add_custom_target(mpi ALL "")
set(mpi_lib MPI::MPI_Fortran)
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${MPI_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES} Threads::Threads)

#------------------------------------------------------------------------------
# Looking for OpenMP. 
#------------------------------------------------------------------------------
message(STATUS "---------------------------------------------------------------")
find_package(OpenMP REQUIRED)

add_custom_target(openmp ALL "")
set(openmp_lib OpenMP::OpenMP_Fortran)
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
# running python command at configure time
execute_process(
	COMMAND ${TEMPITA_DIR}/tempita_proc.py < $< > $@
	INPUT_FILE ${comm_hdf_mod}.in
	OUTPUT_FILE ${comm_hdf_mod}
	)

add_custom_target(required_libraries ALL "" 
	DEPENDS tempita 
					mpi
					openmp
					#blas
					#zlib
					)

#------------------------------------------------------------------------------
# Looking for CFITSIO and (optionally) its dependencies.
#------------------------------------------------------------------------------
#message(STATUS "---------------------------------------------------------------")
#if(NOT (CFITSIO_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#	find_package(CFITSIO 3.470)
#endif()
#if(NOT CFITSIO_FOUND AND CFITSIO_USE_CURL)
#	#------------------------------------------------------------------------------
#	# cURL
#	#------------------------------------------------------------------------------
#	message(STATUS "---------------------------------------------------------------")
#	# looking for cURL in the system. 
#	if(NOT (CURL_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#		# CMake configure scripts foesn't work properly,
#		# so we look for cURL in a standard manner.
#		set(CURL_NO_CURL_CMAKE ON)
#		find_package(CURL)
#	endif()
#
#	# If there is no cURL, we install MbedTLS, LibSSH2 and cURL
#	if(NOT CURL_FOUND)
#		set(CURL_INCLUDE_DIR
#			"${CMAKE_INSTALL_PREFIX}/include"
#			)
#		set(CURL_LIBRARIES
#			"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}curl${CMAKE_SHARED_LIBRARY_SUFFIX}" 
#			)
#		#------------------------------------------------------------------------------
#		message(STATUS "cURL LIBRARIES will be ${CURL_LIBRARIES}")
#		message(STATUS "cURL INCLUDE DIR will be ${CURL_INCLUDE_DIR}")
#		#------------------------------------------------------------------------------
#		# MbedTLS
#		#------------------------------------------------------------------------------
#		message(STATUS "---------------------------------------------------------------")
#		if(NOT (MDEBTLS_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#			find_package(MbedTLS)
#		endif()
#		if(NOT MBEDTLS_FOUND) 
#			message(STATUS "MbedTLS will be compiled from source.")
#			#------------------------------------------------------------------------------
#			message(STATUS "MbedTLS LIBRARIES will be: ${MBEDTLS_LIBRARIES}")
#			message(STATUS "MbedTLS INCLUDE DIRS will be: ${MBEDTLS_INCLUDE_DIRS}")
#			#------------------------------------------------------------------------------
#		else()
#			#------------------------------------------------------------------------------
#			message(STATUS "MbedTLS LIBRARIES are: ${MBEDTLS_LIBRARIES}")
#			message(STATUS "MbedTLS INCLUDE DIRS are: ${MBEDTLS_INCLUDE_DIRS}")
#			#------------------------------------------------------------------------------
#		endif()
#		#------------------------------------------------------------------------------
#		# LibSSH2
#		#------------------------------------------------------------------------------
#		message(STATUS "---------------------------------------------------------------")
#		if(NOT (LIBSSH2_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#			find_package(LibSSH2)
#		endif()
#		if(NOT LIBSSH2_FOUND) 
#			#------------------------------------------------------------------------------
#			message(STATUS "LibSSH2 LIBRARY will be: ${LIBSSH2_LIBRARY}")
#			message(STATUS "LibSSH2 INCLUDE DIR will be: ${LIBSSH2_INCLUDE_DIR}")
#			#------------------------------------------------------------------------------
#		else()
#			#------------------------------------------------------------------------------
#			message(STATUS "LibSSH2 LIBRARY are: ${LIBSSH2_LIBRARY}")
#			message(STATUS "LibSSH2 INCLUDE DIR are: ${LIBSSH2_INCLUDE_DIR}")
#			#------------------------------------------------------------------------------
#		endif()
#	else()
#		#------------------------------------------------------------------------------
#		message(STATUS "cURL LIBRARIES are ${CURL_LIBRARIES}")
#		message(STATUS "cURL INCLUDE DIR is ${CURL_INCLUDE_DIR}")
#		#------------------------------------------------------------------------------
#	endif()
#endif()
