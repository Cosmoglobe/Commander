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
# Description: This script compiles/installs Commander3 on host system and it also 
# links Commander3 to external dependencies (such as HDF5, CFitsio, HEALPix etc.)
#================================================================================

message(STATUS "---------------------------------------------------------------")

# To resolve conflict with the flags, we need to compile *.cpp
# into a stand alone library and then just link it to rest of 
# the files. However, with CMake it is possible to compile the
# file into object, which will not create *.a files. We go with
# this approach.
#add_library(comm_system_backend 
#	STATIC 
#	${COMMANDER3_SOURCE_DIR}/comm_system_backend.cpp
#	)
#target_compile_options(comm_system_backend
#	PRIVATE
#	${COMMANDER3_CXX_COMPILER_FLAGS}
#	)
## Installing `comm_system_backend` as a library
#install(TARGETS comm_system_backend ARCHIVE DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#	)

# TODO: add all sources manually instead of this command, as 
# there seems to be a problem with tempita language
#file(GLOB_RECURSE sources *.f90 *.cpp *.f)
set(sources
	${COMMANDER3_SOURCE_DIR}/commander.f90
	${COMMANDER3_SOURCE_DIR}/comm_map_mod.f90
	${COMMANDER3_SOURCE_DIR}/sharp.f90
	)

# Setting executable name
set(commander3 commander3)
add_executable(${commander3} "")
# make sure that commander executable will be built last
add_dependencies(${commander3} ${projects}) #fftw_float fftw_double)
target_sources(${commander3}
	PUBLIC	
	${sources}
	)
set_property(TARGET ${commander3} PROPERTY ENABLE_EXPORTS TRUE)
set_property(TARGET ${commander3} PROPERTY LINKER_LANGUAGE Fortran)
# Resolving Preprocessor statements for a given Compiler Toolchain
set_source_files_properties(
	${sources}	
	PROPERTIES Fortran_PREPROCESS ON
)
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	target_compile_definitions(${commander3}
		PUBLIC
		USE_INTEL
		)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	target_compile_definitions(${commander3}
		PUBLIC
		USE_GNU
		)
endif()
# adding compiler flags to commander3 target
target_compile_options(${commander3}
	PUBLIC
	# setting flags depending on configuration
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_COMPILER_FLAGS}
	)
# adding linker flags to commander3 target
target_link_options(${commander3}
	PUBLIC
	# setting flags depending on configuration
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_LINKER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_LINKER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_LINKER_FLAGS}
	)

# LINKING ORDER IN LIBRARIES IS IMPORTANT!
# Order is:
# MPI => OpenMP => Blas => LAPACK => HEALPix => 
# CFITSIO => cURL => libm => dl => HDF5 => ZLib
# => FFTW => comm_system_backend
target_link_libraries(${commander3} 
	PRIVATE
	MPI::MPI_Fortran
	OpenMP::OpenMP_Fortran
	${BLAS_LINKER_FLAGS} 
	${BLAS_LIBRARIES}
	${LAPACK_LINKER_FLAGS} 
	${LAPACK_LIBRARIES}
	${HEALPIX_LIBRARIES}
	${CFITSIO_LIBRARIES}
	${LIBM_LIBRARY}
	${CMAKE_DL_LIBS}
	#comm_system_backend
	)

# Installing Commander3 into appropriate folder
#install(TARGETS ${commander3} RUNTIME DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${commander3} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
