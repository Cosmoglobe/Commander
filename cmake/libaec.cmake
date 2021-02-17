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
# Description: This script determines the location of SZip on the host system.
# If it fails to do so, it will download, compile and install SZip from source.
# LibAEC is (not strictly) required by HDF5.
#================================================================================

message(STATUS "---------------------------------------------------------------")

#if(NOT (LIBAEC_FORCE_COMPILE OR HDF5_FORCE_COMPILE OR ALL_FORCE_COMPILE))
if(USE_SYSTEM_LIBAEC AND USE_SYSTEM_LIBS AND NOT USE_SYSTEM_HDF5)
	# TODO: Add maybe SZip with find_package or something like that.
	# Need to have a proper find package or something like that for SZip/LibAEC
	#set(zlib_minimal_accepted_version "1.2.11")
	#find_package(ZLIB 1.2.11)
	# Require ZLib to be of the most recent version
	#if(ZLIB_VERSION_STRING VERSION_LESS_EQUAL zlib_minimal_accepted_version)
	#	message(STATUS "Required version -- ${zlib_minimal_accepted_version} -- will be compiled from source.")
	#endif()
endif()

if(NOT (HDF5_FOUND AND LIBAEC_FOUND)) #OR (ZLIB_VERSION_STRING VERSION_LESS_EQUAL zlib_minimal_accepted_version)) 
	#message(STATUS "Required version -- 1.2.11 -- will be compiled from source.")
	#------------------------------------------------------------------------------
	# Getting LIBAEC from source.
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}_src
		DEPENDS required_libraries 
		GIT_REPOSITORY "${${project}_git_url}"
		GIT_TAG "${${project}_git_tag}"
		#URL "${${project}_url}"
		#URL_MD5 "${${project}_md5}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		#SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		#INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_DOWNLOAD ON
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		)
	#------------------------------------------------------------------------------
	# Building Static LibAEC
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}_static
		DEPENDS required_libraries 
						${project}_src
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		# commands how to build the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			# some systems can create lib64 instead of lib and we need to ensure 
			# we get only lib, otherwise everything breaks.
			-DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
			-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
			# Specifying compilers
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
			# Building only shared libraries
			-DBUILD_SHARED_LIBS:BOOL=OFF
		)
	#------------------------------------------------------------------------------
	# Building Shared LibAEC
	#------------------------------------------------------------------------------
	ExternalProject_Add(${project}_shared
		DEPENDS required_libraries 
						${project}_src
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		SOURCE_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}_src"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		LOG_DIR "${CMAKE_LOG_DIR}"
		LOG_CONFIGURE ON
		LOG_BUILD ON
		LOG_INSTALL ON
		# Disabling download
		DOWNLOAD_COMMAND ""
		# commands how to build the project
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=Release
			# Specifying installations paths for binaries and libraries
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			# some systems can create lib64 instead of lib and we need to ensure 
			# we get only lib, otherwise everything breaks.
			-DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
			-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
			# Specifying compilers
			-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
			# Building only shared libraries
			-DBUILD_SHARED_LIBS:BOOL=ON
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
	# Adding new target -- libaec -- to ensure that only after all libraries built
	# the project can use this target.
	add_custom_target(${project} 
		ALL ""
		DEPENDS ${project}_static
						${project}_shared
		)
	#------------------------------------------------------------------------------
	message(STATUS "LIBAEC LIBRARIES will be: ${LIBAEC_LIBRARIES}")
	message(STATUS "LIBAEC INCLUDE DIRS will be: ${LIBAEC_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
else()
	add_custom_target(${project} ALL "")
	#------------------------------------------------------------------------------
	message(STATUS "LIBAEC LIBRARIES are: ${LIBAEC_LIBRARIES}")
	message(STATUS "LIBAEC INCLUDE DIRS are: ${LIBAEC_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
