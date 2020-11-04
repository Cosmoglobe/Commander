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
# Description: This script determines the location of Doxygen on the host system.
# If it fails to do so, it will download, compile and install Doxygen from source.
# In addition, it tries to determine the location of doxygen dependencies - Flex
# and Bison. If it fails to do so, it will install them from source as well.
# Note: Doxygen is not required dependency, but it is used to create out-of-source
# documentation. 
#================================================================================

if(DOXYGEN_BUILD_DOCS)
	message(STATUS "---------------------------------------------------------------")
	set(ENV{PATH}
		 #$ENV{PATH}:${out_install_dir}/include/:${out_lib_dir}/:${out_bin_dir}/curl
		 $ENV{PATH}:${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
		 )

	if(NOT FLEX_FORCE_COMPILE)
		find_package(FLEX)
	endif()
	if(NOT FLEX_FOUND)
		message(STATUS "Will download flex from source.")
		# Creating configure command for Flex
		set(flex_configure_command 
			"${CMAKE_COMMAND}" "-E" "env" 
			"FC=${COMMANDER3_Fortran_COMPILER}" 
			"CXX=${COMMANDER3_CXX_COMPILER}" 
			"CPP=${COMMANDER3_CPP_COMPILER}" 
			"CC=${COMMANDER3_C_COMPILER}" 
			"./configure" 
			"--prefix=<INSTALL_DIR>"
			)
		#------------------------------------------------------------------------------
		# Getting Flex from source
		ExternalProject_Add(flex
			URL ${flex_url}
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/flex"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/flex/src/flex"
			INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
			LOG_DIR "${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD ON
			LOG_CONFIGURE ON
			LOG_BUILD ON
			LOG_INSTALL ON
			CONFIGURE_COMMAND "${flex_configure_command}"
			)
		set(FLEX_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/flex)
		message(STATUS "${FLEX_EXECUTABLE}")
	else()
		message(STATUS "Found Flex")
		add_custom_target(flex ALL "")
	endif()

	if(NOT BISON_FORCE_COMPILE)
		find_package(BISON)
	endif()
	if(NOT BISON_FOUND)
		message(STATUS "Will download Bison from source")
		# Creating configure command for Bison
		set(bison_configure_command 
			"${CMAKE_COMMAND}" "-E" "env" 
			"FC=${COMMANDER3_Fortran_COMPILER}" 
			"CXX=${COMMANDER3_CXX_COMPILER}" 
			"CPP=${COMMANDER3_CPP_COMPILER}" 
			"CC=${COMMANDER3_C_COMPILER}" 
			"./configure" 
			"--prefix=<INSTALL_DIR>"
			)
		#------------------------------------------------------------------------------
		# Getting Bison from source
		ExternalProject_Add(bison
			URL ${bison_url}
			PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/bison"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/bison/src/bison"
			INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
			LOG_DIR "${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD ON
			LOG_CONFIGURE ON
			LOG_BUILD ON
			LOG_INSTALL ON
			CONFIGURE_COMMAND "${bison_configure_command}"
			)
		set(BISON_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bison)	
	else()
		message(STATUS "Found Bison")
		add_custom_target(bison ALL "")
	endif()

	if(NOT DOXYGEN_FORCE_COMPILE)
		find_package(Doxygen)
	endif()
	if(NOT DOXYGEN_FOUND)
		message(STATUS "Will download Doxygen from source.")
		#------------------------------------------------------------------------------
		# Getting Doxygen from source
		ExternalProject_Add(${project}
			DEPENDS flex bison
			URL "${${project}_url}"
			#GIT_REPOSITORY "${${project}_url}"
			DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}" #"${download_dir}"
			BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
			INSTALL_DIR ${CMAKE_INSTALL_PREFIX} #"${out_install_dir}"
			LOG_DIR "${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD ON
			LOG_CONFIGURE ON
			LOG_BUILD ON
			LOG_INSTALL ON
			CMAKE_ARGS
			-DFLEX_EXECUTABLE:PATH=${FLEX_EXECUTABLE}
			-DBISON_EXECUTABLE:PATH=${BISON_EXECUTABLE}
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			-DCMAKE_BUILD_TYPE=Release
			)
	else()
		add_custom_target(${project} ALL "")
		# refer to cmake docs about doxygen variables and options
		message(STATUS "Building Commander3 documentation from source code with Doxygen.")
		set(DOXYGEN_PROJECT_NAME "Commander3")
		set(DOXYGEN_PROJECT_NUMBER "${CMAKE_PROJECT_VERSION}")
		set(DOXYGEN_PROJECT_BRIEF "${CMAKE_PROJECT_DESCRIPTION}")
		set(DOXYGEN_GENERATE_HTML YES)
		set(DOXYGEN_GENERATE_LATEX YES)
		set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/docs)
		#
		set(DOXYGEN_FULL_PATH_NAMES NO)
		# Set the OPTIMIZE_FOR_FORTRAN tag to YES if your project consists of Fortran
		# sources. Doxygen will then generate output that is tailored for Fortran.
		# The default value is: NO.
		set(DOXYGEN_OPTIMIZE_FOR_FORTRAN YES)
		# If the EXTRACT_ALL tag is set to YES, doxygen will assume all entities in
		# documentation are documented, even if no documentation was available. Private
		# class members and static file members will be hidden unless the
		# EXTRACT_PRIVATE respectively EXTRACT_STATIC tags are set to YES.
		# Note: This will also disable the warnings about undocumented members that are
		# normally produced when WARNINGS is set to YES.
		# The default value is: NO.
		set(DOXYGEN_EXTRACT_ALL YES)
		
		doxygen_add_docs(${project}_docs
			ALL
			${COMMANDER3_SOURCE_DIR}
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
			)
	endif()

else()
	# making an empty target so the project will compile regardless of doxygen
	add_custom_target(${project} ALL "")
endif()

#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
