# Project: HDF5 
# File which contains setup for current project 
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
# asking for an exact hdf5 version
#find_package(HDF5 1.10.5 EXACT COMPONENTS Fortran) #Fortran_HL)

if(USE_EXISTING_HDF5)
find_package(HDF5 1.10.0 COMPONENTS Fortran) #Fortran_HL)
endif()

if(NOT HDF5_FOUND)
	message(STATUS "Will download HDF5v1.10.5 from source.")
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
		#COMMAND ${CMAKE_COMMAND} -E env FC=${CMAKE_Fortran_COMPILER} CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} MPCC=${CMAKE_C_COMPILER} MPFC=${CMAKE_Fortran_COMPILER} MPCXX=${CMAKE_CXX_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-fortran #--enable-cxx #--enable-parallel
		#COMMAND ${CMAKE_COMMAND} -E env FC=${MPI_Fortran_COMPILER} CXX=${MPI_CXX_COMPILER} CPP=${COMMANDER3_CPP_COMPILER} CC=${MPI_C_COMPILER} ./configure --prefix=<INSTALL_DIR> --enable-fortran #--enable-parallel #--enable-cxx #--enable-parallel
		)
	# adding hdf5_fortran and hdf5 into a variable (to be consistent with cmake docs)
	set(HDF5_Fortran_LIBRARIES 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhdf5_fortran.a" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhdf5.a" 
		"${ZLIB_LIBRARIES}")
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
	message(STATUS "FOUND HDF5 Fortran!")
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	message(STATUS "HDF5 Fortran INCLUDE DIRS are: ${HDF5_Fortran_INCLUDE_DIRS}")
	message(STATUS "HDF5 Fortran DEFINITIONS are: ${HDF5_Fortran_DEFINITIONS}")
	message(STATUS "HDF5 Fortran LIBRARIES are: ${HDF5_Fortran_LIBRARIES}")
endif()
