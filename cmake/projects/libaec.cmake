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
# Description: This script determines the location of SZip on the host system.
# If it fails to do so, it will download, compile and install SZip from source.
# LibAEC is (not strictly) required by HDF5.
#================================================================================

#message(STATUS "---------------------------------------------------------------")
#
#if(USE_SYSTEM_LIBAEC AND USE_SYSTEM_LIBS AND NOT USE_SYSTEM_HDF5)
#	# TODO: Add maybe SZip with find_package or something like that.
#	# Need to have a proper find package or something like that for SZip/LibAEC
#	#set(zlib_minimal_accepted_version "1.2.11")
#	#find_package(ZLIB 1.2.11)
#	# Require ZLib to be of the most recent version
#	#if(ZLIB_VERSION_STRING VERSION_LESS_EQUAL zlib_minimal_accepted_version)
#	#	message(STATUS "Required version -- ${zlib_minimal_accepted_version} -- will be compiled from source.")
#	#endif()
#endif()

if(COMPILE_LIBAEC)
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting LIBAEC from source.
	#------------------------------------------------------------------------------
	if(NOT EXISTS "${LIBAEC_SOURCE_DIR}/CMakeLists.txt")
		message(STATUS "No LIBAEC sources were found; thus, will download it from source:\n${libaec_git_url}")
		ExternalProject_Add(
			libaec_src
			GIT_REPOSITORY		"${libaec_git_url}"
			GIT_TAG						"${libaec_git_tag}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${LIBAEC_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
		message(STATUS "Found an existing LIBAEC sources inside:\n${LIBAEC_SOURCE_DIR}")
		add_custom_target(libaec_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and Installing Static and Shared LibAEC
	#------------------------------------------------------------------------------
	list(APPEND _AEC_LIB_TYPES_ static shared)
	list(APPEND _AEC_LIB_BOOL_VALS_ -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON)
	foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _AEC_LIB_TYPES_ _AEC_LIB_BOOL_VALS_)
		ExternalProject_Add(
			libaec_${_lib_type_}
			DEPENDS						required_libraries 
												libaec_src
			PREFIX						"${LIBS_BUILD_DIR}"
			SOURCE_DIR				"${LIBAEC_SOURCE_DIR}"
			INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_CONFIGURE			ON
			LOG_BUILD					ON
			LOG_INSTALL				ON
			# Disabling download
			DOWNLOAD_COMMAND	""
			# commands how to build the project
			CMAKE_ARGS
				-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
				# Specifying installations paths for binaries and libraries
				-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
				-DCMAKE_INSTALL_LIBDIR:PATH=lib #${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
				-DCMAKE_INSTALL_INCLUDEDIR:PATH=include
				# Specifying compilers
				-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
				# Building only shared libraries
				${_bool_val_}
			)
	endforeach()
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(libaec 
		ALL ""
		DEPENDS libaec_static
						libaec_shared
		)
	#------------------------------------------------------------------------------
	# Defining the variable which will show the path to the compiled libraries
	set(LIBAEC_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include"
		)
	include_directories(${LIBAEC_INCLUDE_DIRS}) 
	set(SZIP_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}sz${CMAKE_SHARED_LIBRARY_SUFFIX}" 
		)
	set(LIBAEC_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}aec${CMAKE_SHARED_LIBRARY_SUFFIX}" 
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}sz${CMAKE_SHARED_LIBRARY_SUFFIX}" 
		)
	#------------------------------------------------------------------------------
	#message(STATUS "LIBAEC LIBRARIES will be: ${LIBAEC_LIBRARIES}")
	#message(STATUS "LIBAEC INCLUDE DIRS will be: ${LIBAEC_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	add_custom_target(${project} ALL "")
	#------------------------------------------------------------------------------
	#message(STATUS "LIBAEC LIBRARIES are: ${LIBAEC_LIBRARIES}")
	#message(STATUS "LIBAEC INCLUDE DIRS are: ${LIBAEC_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
