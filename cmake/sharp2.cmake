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
# Description: This script determines the location of Libsharp2 on the host system.
# If it fails to do so, it will download, compile and install Libsharp2 from source.
# Note: Starting from version 3.60, HEALPix comes bundled together with Libsharp2.
# Therefore, this file is depricated and will be removed in future.
#================================================================================

message(STATUS "---------------------------------------------------------------")

find_package(SHARP2)
if(NOT SHARP2_FOUND)
	# Creating configure command for SHARP2 
	set(sharp2_C_FLAGS "-DUSE_MPI -std=c99 -O3 -ffast-math")
	set(sharp2_configure_command 
		"autoreconf" "-i" "&&" 
		"${CMAKE_COMMAND}" "-E" "env" 
		"FC=${COMMANDER3_Fortran_COMPILER}" 
		"CXX=${COMMANDER3_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${COMMANDER3_C_COMPILER}" 
		"CFLAGS=${sharp2_C_FLAGS}" 
		"./configure"
		)
	#------------------------------------------------------------------------------
	# Getting Sharp2 from source
	ExternalProject_Add(${project}
		URL "${${project}_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${CMAKE_INSTALL_PREFIX}"
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		INSTALL_COMMAND ""
		# LibSharp doesn't have an install command, so we need to manually move files
		# into installation folder and hope they will work :)
		COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}" "${CMAKE_INSTALL_PREFIX}/${project}" #"${out_install_dir}/lib" 
		)

	# defining the variable which will show the path to the compiled libraries
	set(SHARP2_LIBRARIES ${CMAKE_INSTALL_PREFIX}/${project}/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
	message(STATUS "SHARP2 LIBRARIES will be: ${SHARP2_LIBRARIES}")
else()
	add_custom_target(${project} ALL "")
	message(STATUS "SHARP2 LIBRARIES are: ${SHARP2_LIBRARIES}")
	message(STATUS "SHARP2 INCLUDE DIRS are: ${SHARP2_INCLUDE_DIR}")
endif()
