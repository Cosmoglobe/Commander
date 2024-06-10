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
# Description: This script determines the location of HDF5 on the host system.
# If it fails to do so, it will download, compile and install HDF5 from source.
# The HDF5 group provides sources for zlib and szip, which we are going to use
# here as well.
#================================================================================

#message(STATUS "---------------------------------------------------------------")

# TODO: make another variable for shared/static linking
# also ensure that if hdf5 wasn't compiled with autotools
# it still be working as before.
#if(NOT (HDF5_FORCE_COMPILE OR ALL_FORCE_COMPILE))
#if(USE_SYSTEM_HDF5 AND USE_SYSTEM_LIBS)
#	# Using static linking instead of dynamic
#	set(HDF5_USE_STATIC_LIBRARIES FALSE)#TRUE)
#	# Using parallel build instead of serial
#	set(HDF5_PREFER_PARALLEL TRUE)
#	#find_package(HDF5 1.10.0 COMPONENTS Fortran) # Fortran_HL)
#	find_package(HDF5 1.12.0 COMPONENTS Fortran Fortran_HL)
#endif()

if(COMPILE_HDF5)

	# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP
  #if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
	#	cmake_policy(SET CMP0135 NEW)
	#endif()

	if(NOT HDF5_Fortran_FOUND)
		message(STATUS "Missing component -- Fortran -- will be compiled from source.")	
	endif()
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting HDF5 from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${HDF5_SOURCE_DIR}/CMakeLists.txt")
		message(STATUS "No HDF5 sources were found; thus, will download it from source:\n${hdf5_url}")
		ExternalProject_Add(
			hdf5_src
			DEPENDS						required_libraries 
												zlib 
												libaec
			URL								"${hdf5_url}"
			URL_MD5						"${hdf5_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${HDF5_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing HDF5 sources inside:\n${HDF5_SOURCE_DIR}")
		add_custom_target(hdf5_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and Installing Static and Shared HDF5
	#------------------------------------------------------------------------------
	ExternalProject_Add(
		hdf5
		DEPENDS						required_libraries 
											zlib 
											libaec
											hdf5_src
		PREFIX						"${LIBS_BUILD_DIR}"
		SOURCE_DIR				"${HDF5_SOURCE_DIR}"
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
		LOG_CONFIGURE			ON 
		LOG_BUILD					ON 
		LOG_INSTALL				ON 
		# commands how to build the project
		DOWNLOAD_COMMAND	""
		CMAKE_ARGS
      #-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			#-DCMAKE_Fortran_MODULE_DIRECTORY:PATH=${CMAKE_Fortran_MODULE_DIRECTORY}
			# Specifying compilers
			-DCMAKE_Fortran_COMPILER=${MPI_Fortran_COMPILER}
			-DCMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
			-DCMAKE_C_COMPILER=${MPI_C_COMPILER}
			# Building both static and shared libraries
			# Note: On MacOS shared libraries are not supported!
			-DBUILD_SHARED_LIBS:BOOL=ON
			-DBUILD_STATIC_LIBS:BOOL=ON
			# Instructions for building language specific libraries
			-DHDF5_BUILD_FORTRAN:BOOL=ON
			-DHDF5_ENABLE_PARALLEL:BOOL=ON
			-DHDF5_BUILD_JAVA:BOOL=OFF
			-DHDF5_BUILD_CPP_LIB:BOOL=OFF
			-DHDF5_ENABLE_THREADSAFE:BOOL=OFF
			-DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=ON
			# Enabling ZLIB support
			# Note: for some reason this configuration works, even if we get
			# -- Could NOT find ZLIB (missing: ZLIB_DIR)
			# -- Found ZLIB: /usr/lib64/libz.so (found version "1.2.7")
			# -- Filter ZLIB is ON
			# during hdf5 configuration phase. Works while building only static libs
			-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
			-DZLIB_USE_EXTERNAL:BOOL=OFF
			-DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIRS}
			-DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARIES}
			# Enable SZip support
			-DHDF5_ENABLE_SZIP_SUPPORT:BOOL=ON
			-DHDF5_ENABLE_SZIP_ENCODING:BOOL=ON
			#-DSZIP_LIBRARY:FILEPATH=${LIBAEC_LIBRARIES}
			-DSZIP_LIBRARY:FILEPATH=${SZIP_LIBRARIES}
			-DSZIP_INCLUDE_DIR:PATH=${LIBAEC_INCLUDE_DIRS}
		)
	# adding hdf5_fortran and hdf5 into a variable (to be consistent with cmake docs)
	set(HDF5_Fortran_LIBRARIES 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_fortran${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_f90cstub${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_fortran${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_f90cstub${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		#"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5${CMAKE_STATIC_LIBRARY_SUFFIX}" 
		"${ZLIB_LIBRARIES}"
		"${SZIP_LIBRARIES}"
		)
	#hdf5_hl_fortran hdf5_hl_f90cstub hdf5_fortran hdf5_f90cstub hdf5_hl hdf5
	#[=[
		When compiled with CMake, HDF5 creates additional folders inside include
		to store .mod/.o files. We need to add them as well to compile sucessfully.	
	#]=]
	set(HDF5_Fortran_INCLUDE_DIRS 
			"${CMAKE_INSTALL_PREFIX}/include"
			#"${CMAKE_Fortran_MODULE_DIRECTORY}/static"
			"${CMAKE_INSTALL_PREFIX}/include/static"
			)
	include_directories(${HDF5_Fortran_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "HDF5 Fortran LIBRARIES will be: ${HDF5_Fortran_LIBRARIES}")
	#message(STATUS "HDF5 Fortran INCLUDE DIRS will be: ${HDF5_Fortran_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	add_custom_target(hdf5
		ALL ""
		DEPENDS required_libraries 
						zlib
						libaec
		)
	# When compiling with CMake, HDF5_Fortran_INCLUDE_DIRS appears
	# to be empty, but HDF5_INCLUDE_DIRS is not. However, we still need
	# to specify "static"/"shared"directory, which contains .mod/.o files.
	set(HDF5_Fortran_INCLUDE_DIRS ${HDF5_Fortran_INCLUDE_DIRS}
		${HDF5_INCLUDE_DIRS}
		${HDF5_INCLUDE_DIRS}/shared #static
		)
	include_directories(${HDF5_Fortran_INCLUDE_DIRS} 
		#${HDF5_INCLUDE_DIRS}
		#${HDF5_INCLUDE_DIRS}/shared #static
		)
	#message(STATUS "HDF5 Fortran INCLUDE DIRS are: ${HDF5_Fortran_INCLUDE_DIRS}")
	#message(STATUS "HDF5 Fortran LIBRARIES are: ${HDF5_Fortran_LIBRARIES}")
	#message(STATUS ${HDF5_Fortran_DEFINITIONS})
	#message(STATUS ${HDF5_Fortran_LIBRARY})
	#message(STATUS ${HDF5_LIBRARY})
	#message(STATUS ${HDF5_LIBRARIES})
	#message(STATUS ${HDF5_Fortran_LIBRARY_NAMES})
	#message(STATUS ${HDF5_INCLUDE_DIR})
	#message(STATUS ${HDF5_INCLUDE_DIRS})
endif()
