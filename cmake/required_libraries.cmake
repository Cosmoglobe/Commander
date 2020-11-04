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

# We will be using Git to download some dependencies, so we need to check if git available
find_package(Git REQUIRED)

message(STATUS "---------------------------------------------------------------")
find_package(ZLIB REQUIRED)

message(STATUS "ZLIB LIBRARIES ARE: ${ZLIB_LIBRARIES}")
message(STATUS "ZLIB INCLUDE DIRS ARE: ${ZLIB_INCLUDE_DIRS}")

add_custom_target(zlib ALL "")
include_directories(${ZLIB_INCLUDE_DIRS})
add_library(zlib_lib SHARED IMPORTED GLOBAL) 
set_target_properties(zlib_lib PROPERTIES IMPORTED_LOCATION ${ZLIB_LIBRARIES})

message(STATUS "---------------------------------------------------------------")
# require BLAS and LAPACK
# From docs: Note C, CXX or Fortran must be enabled
# to detect a BLAS/LAPACK library. C or CXX must be
# enabled to use Intel Math Kernel Library (MKL).
#set(BLA_VENDOR Intel10_64lp)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# to avoid cmake errors we create and empty target
add_custom_target(blas ALL "")
set(blas_lib ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
message(STATUS "BLAS LINKER FLAGS:   ${BLAS_LINKER_FLAGS}")
message(STATUS "BLAS LIBRARIES:      ${BLAS_LIBRARIES}")
message(STATUS "LAPACK LINKER FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK LIBRARIES:    ${LAPACK_LIBRARIES}")

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

message(STATUS "---------------------------------------------------------------")
find_package(OpenMP REQUIRED)

add_custom_target(openmp ALL "")
set(openmp_lib OpenMP::OpenMP_Fortran)
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})

message(STATUS "OPENMP Fortran LIBRARIES are: ${OpenMP_Fortran_LIBRARIES}")

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
					blas
					zlib
					)
