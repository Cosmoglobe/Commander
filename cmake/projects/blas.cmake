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
# Description: This script determines the location of MKL/OpenBLAS on the host system.
# If it fails to do so, it will download, compile and install OpenBLAS from source.
#================================================================================
#message(STATUS "---------------------------------------------------------------")
# require BLAS and LAPACK
# From docs: Note C, CXX or Fortran must be enabled
# to detect a BLAS/LAPACK library. C or CXX must be
# enabled to use Intel Math Kernel Library (MKL).
# Note: Because native (shipped with Linux) BLAS/LAPACK
# implementations are not optimized, we require usage of
# either MKL or OpenBLAS. 

#if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
# This works for MKL
#	set($ENV{BLA_VENDOR}
#		#OpenBLAS
#			Intel10_32
#			Intel10_64lp
#			Intel10_64lp_seq
#			Intel10_64ilp
#			Intel10_64ilp_seq
#			Intel10_64_dyn
#			)
#else()
#if(NOT (BLAS_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#if(USE_SYSTEM_BLAS AND USE_SYSTEM_LIBS)
	# This works for OpenBLAS
	# Note: Sometimes this doesn't work, i.e. it cannot detect MKL/OpenBLAS 
	# for some weird reason. In this case it is a good idea to logout and login
	# to refresh terminal.
#	set($ENV{BLA_VENDOR} 
#			OpenBLAS
#			Intel10_32
#			Intel10_64lp
#			Intel10_64lp_seq
#			Intel10_64ilp
#			Intel10_64ilp_seq
#			Intel10_64_dyn
#			)
#	find_package(BLAS) #REQUIRED)
#	find_package(LAPACK) #REQUIRED)
#endif()

# Download OpenBLAS from source if neither MKL or OpenBLAS
# were found on the host system.
#if(NOT (BLAS_FOUND OR LAPACK_FOUND))
if(COMPILE_BLAS)
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting OpenBLAS from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${BLAS_SOURCE_DIR}/CMakeLists.txt")
		message(STATUS "No BLAS sources were found; thus, will download it from source:\n${blas_url}")
		ExternalProject_Add(
			blas_src
			DEPENDS						required_libraries
			URL								"${blas_url}"
			URL_MD5						"${blas_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${BLAS_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing BLAS sources inside:\n${BLAS_SOURCE_DIR}")
		add_custom_target(blas_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and installing OpenBLAS
	#------------------------------------------------------------------------------
	ExternalProject_Add(
		blas
		DEPENDS						required_libraries
											blas_src
		URL								"${blas_url}"
		URL_MD5						"${blas_md5}"
		PREFIX						"${LIBS_BUILD_DIR}"
		SOURCE_DIR				"${BLAS_SOURCE_DIR}"
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_CONFIGURE 		ON
		LOG_BUILD					ON
		LOG_INSTALL				ON
		# commands how to build the project
		DOWNLOAD_COMMAND	""
		CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DUSE_OPENMP=ON
			-DCMAKE_Fortran_COMPILER=${MPI_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${MPI_C_COMPILER}
			#-DCMAKE_INSTALL_LIBDIR=lib
			-DCMAKE_INSTALL_LIBDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
		)
	#------------------------------------------------------------------------------
	# In case of static linking, we do not need to specify linker flags.
	# NOTE: If OpenBLAS 0.3.12 is installed with CMake, there is a special
	# variable, which controls the installation of libraries and it is 
	# CMAKE_INSTALL_LIBDIR. We need to set it manually as above to "lib"
	# otherwise on some platforms the dir will be named "lib" and on others
	# "lib64".
	set(BLAS_LINKER_FLAGS "")
	set(BLAS_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	set(LAPACK_LINKER_FLAGS "")
	set(LAPACK_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	#------------------------------------------------------------------------------
	#message(STATUS "BLAS LINKER FLAGS will be:   ${BLAS_LINKER_FLAGS}")
	#message(STATUS "BLAS LIBRARIES will be:      ${BLAS_LIBRARIES}")
	#message(STATUS "LAPACK LINKER FLAGS will be: ${LAPACK_LINKER_FLAGS}")
	#message(STATUS "LAPACK LIBRARIES will be:    ${LAPACK_LIBRARIES}")
	#------------------------------------------------------------------------------
else()
	# to avoid cmake errors we create and empty target
	add_custom_target(blas 
		ALL ""
		DEPENDS required_libraries
		)
	set(blas_lib ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
	#------------------------------------------------------------------------------
	#message(STATUS "BLAS LINKER FLAGS:   ${BLAS_LINKER_FLAGS}")
	#message(STATUS "BLAS LIBRARIES:      ${BLAS_LIBRARIES}")
	#message(STATUS "LAPACK LINKER FLAGS: ${LAPACK_LINKER_FLAGS}")
	#message(STATUS "LAPACK LIBRARIES:    ${LAPACK_LIBRARIES}")
	#------------------------------------------------------------------------------
endif()
