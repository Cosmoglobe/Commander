# Project: HDF5 
# File which contains setup for current project 
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
# asking for an exact hdf5 version
#find_package(HDF5 1.10.5 EXACT COMPONENTS Fortran) #Fortran_HL)

if(NOT HDF5_FORCE_COMPILE)
	find_package(HDF5 1.10.0 COMPONENTS Fortran) #Fortran_HL)
endif()

if(NOT HDF5_FOUND)
	if(NOT HDF5_Fortran_FOUND)
		message(STATUS "Missing component - Fortran - will be compiled from source")	
	endif()
	# Creating configure command for HDF5 autocofig (can be compiled with CMake as well) 
	list(APPEND hdf5_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
		)
	# HDF5 has problems with PGI compilers (couldn't identify correct flags)
	# so we include them here manually, e.g. -fPIC
	if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
		set(hdf5_CFLAGS "-fPIC")
		set(hdf5_FCFLAGS "-fPIC")
		set(hdf5_CXXFLAGS "-fPIC")
		list(APPEND hdf5_configure_command
			"CFLAGS=${hdf5_CFLAGS}"
			"CXXFLAGS=${hdf5_CXXFLAGS}"
			"FCFLAGS=${hdf5_FCFLAGS}"
			)
	endif()
	list(APPEND hdf5_configure_command
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"./configure" 
		"--prefix=<INSTALL_DIR>" 
		"--enable-fortran"
		"--enable-parallel"
		)

	# Getting HDF5 from source
	ExternalProject_Add(${project}
		DEPENDS zlib
		URL "${${project}_url}"
		URL_MD5 "${${project}_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		)
	# adding hdf5_fortran and hdf5 into a variable (to be consistent with cmake docs)
	set(HDF5_Fortran_LIBRARIES 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_fortran${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${ZLIB_LIBRARIES}"
		)
	set(HDF5_Fortran_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include/")
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	message(STATUS "HDF5 Fortran LIBRARIES will be: ${HDF5_Fortran_LIBRARIES}")
	message(STATUS "HDF5 Fortran INCLUDE DIRS will be: ${HDF5_Fortran_INCLUDE_DIRS}")
	#set($ENV{PATH} "${out_lib_dir}/")
else()
	add_custom_target(hdf5 ""
		ALL
		DEPENDS zlib
		)
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	message(STATUS "HDF5 Fortran INCLUDE DIRS are: ${HDF5_Fortran_INCLUDE_DIRS}")
	message(STATUS "HDF5 Fortran LIBRARIES are: ${HDF5_Fortran_LIBRARIES}")
endif()
