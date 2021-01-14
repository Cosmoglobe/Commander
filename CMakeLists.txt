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
# Ensure out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(FATAL_ERROR 
			"\n-------------------------------------------------------------------"
			"\n		IN-SOURCE BUILDS ARE NOT ALLOWED!" 
			"\nTo sucessfully run CMake configuration, please, do the following:"
			"\n1. Delete the CMakeCache.txt file and the CMakeFiles directory"
			" created inside current (root) directory. For instance, on Linux do: "
			"\n		$ rm -r CMakeCache.txt CMakeFiles/"
			"\n2. Create new directory (any name will do, but usually it is called build)"
			" and re-run CMake configuration from it. For example, on Linux do:"
			"\n		$ mkdir build && cd build"
			"\n		$ cmake [CMake variable 1] [...] [CMake variable n] .."
			"\n3. Run compilation command to install Commander3, e.g.: "
			"\n		$ cmake --build . --target install -j n"
			"\n where n is the number of processors to utilize."
			"\nFor more information, please refer to official documentation: "
			"\n		https://docs.beyondplanck.science"
			"\nIn case you have any questions, feel free to post them on our forum:"
			"\n		https://forums.beyondplanck.science/"
			"\n-------------------------------------------------------------------"
			)
endif()
# Just a short message to signify start of configuration script
message(STATUS "===============================================================")
message(STATUS "Welcome to automatic Commander3 installer!")
message(STATUS "Issued under: GNU GPLv3 license")
message(STATUS "Official repo: https://github.com/Cosmoglobe/Commander")
message(STATUS "Documentation: https://docs.beyondplanck.science")
message(STATUS "Forum: https://forums.beyondplanck.science/")
message(STATUS "===============================================================")
message(STATUS "Scanning your system. Please, wait a moment...")
# The necessary minimum to set-up cmake-based project
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(Commander3
	DESCRIPTION "Complete installer for Commander 3"
	VERSION 1.2.0	
	LANGUAGES Fortran CXX C
	)
# Looking for python
find_package(Python REQUIRED COMPONENTS Interpreter)

# Including native cmake commands 
include(CMakePrintHelpers)
include(CMakeDependentOption)
# Enables content population at configure time
include(FetchContent)
# Including external project routines to enable superbuild pattern
include(ExternalProject)
# Verify that Fortran and CXX compilers can properly talk to each other
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
# Check if given fortran source compiles
include(CheckFortranSourceCompiles)
# The -fPIC flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# The -flto

# Adding project modules dir (with .cmake extension)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# Including custom cmake modules
include(project_instructions)
# Adding commander3 source files
# Note: I have added commander.cmake as an "independent" project 
# so there is no need for other cmake lists. Advantage of this 
# approach is that everything located in the same folder (cmake)
#add_subdirectory(commander3)
