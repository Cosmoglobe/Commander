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
# Description: This script determines the location of BLAS/LAPACK on the host system.
# If it fails to do so, it will download, compile and install OpenBLAS from source.
#================================================================================
message(STATUS "---------------------------------------------------------------")
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
# This works for OpenBLAS
# Note: Sometimes this doesn't work, i.e. it cannot detect MKL/OpenBLAS 
# for some weird reason. Int his case it is a good idea to logout an login
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
#endif()
find_package(BLAS) #REQUIRED)
find_package(LAPACK) #REQUIRED)

# Download OpenBLAS from source if neither MKL or OpenBLAS
# were found on the host system.
if(NOT (BLAS_FOUND OR LAPACK_FOUND))
	#list(APPEND blas_configure_command
		#"${CMAKE_COMMAND}" "-E" "env" 
		#"FC=${COMMANDER3_Fortran_COMPILER}" 
		#"CXX=${COMMANDER3_CXX_COMPILER}" 
		#"CPP=${COMMANDER3_CPP_COMPILER}" 
		#"CC=${COMMANDER3_C_COMPILER}" 
		#"export USE_OPENMP=1"
		#"export USE_THREADS=1"
		#"PREFIX=${CMAKE_INSTALL_PREFIX}"
		#)
	#list(APPEND blas_install_command
		#"make"
		#"PREFIX=${CMAKE_INSTALL_PREFIX}"
		#"install"
		#)
	ExternalProject_Add(${project}
		DEPENDS required_libraries
		URL "${${project}_url}"
		URL_MD5 "${${project}_md5}"
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
		#COMMAND "${${project}_configure_command}" 
		#CONFIGURE_COMMAND "${${project}_configure_command}"
		#INSTALL_COMMAND "${blas_install_command}"	
		CMAKE_ARGS
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_BUILD_TYPE=Release
			-DUSE_OPENMP=ON
			-DCMAKE_Fortran_COMPILER=${COMMANDER3_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${COMMANDER3_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${COMMANDER3_C_COMPILER}
			#-DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
			#-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
		)
	#message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	# In case of static linking, we do not need to specify linker flags
	set(BLAS_LINKER_FLAGS "")
	set(BLAS_LIBRARIES
		"${CMAKE_LIBRARY64_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		)
	set(LAPACK_LINKER_FLAGS "")
	set(LAPACK_LIBRARIES
		"${CMAKE_LIBRARY64_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		)
	#message(${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX})
	message(STATUS "BLAS LINKER FLAGS will be:   ${BLAS_LINKER_FLAGS}")
	message(STATUS "BLAS LIBRARIES will be:      ${BLAS_LIBRARIES}")
	message(STATUS "LAPACK LINKER FLAGS will be: ${LAPACK_LINKER_FLAGS}")
	message(STATUS "LAPACK LIBRARIES will be:    ${LAPACK_LIBRARIES}")
else()
	# to avoid cmake errors we create and empty target
	add_custom_target(blas ALL "")
	set(blas_lib ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
	message(STATUS "BLAS LINKER FLAGS:   ${BLAS_LINKER_FLAGS}")
	message(STATUS "BLAS LIBRARIES:      ${BLAS_LIBRARIES}")
	message(STATUS "LAPACK LINKER FLAGS: ${LAPACK_LINKER_FLAGS}")
	message(STATUS "LAPACK LIBRARIES:    ${LAPACK_LIBRARIES}")
endif()
