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

# just a short message to signify start of configuration script
message(STATUS "===============================================================")
message(STATUS "Welcome to automatic Commander3 installer!")
message(STATUS "Issued under: GNU GPLv3 license")
message(STATUS "Official repo: https://github.com/Cosmoglobe/Commander")
message(STATUS "Documentation: https://docs.beyondplanck.science")
message(STATUS "===============================================================")
message(STATUS "Scanning your system. Please, wait a moment...")
#message(STATUS "---------------------------------------------------------------")
# the necessary minimum to set-up cmake
cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(CommanderSuperbuild
	DESCRIPTION "Complete installer for Commander 3"
	VERSION 1.0.0	
	LANGUAGES Fortran CXX C #Python
	)

# ensuring that we do out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

# including native cmake commands 
include(CMakePrintHelpers)
include(CMakeDependentOption)
# enables content population at configure time
include(FetchContent)
# including external project routines to enable superbuild pattern
include(ExternalProject)
# verify that Fortran and CXX compilers can properly talk to each other
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
# check if given fortran source compiles
include(CheckFortranSourceCompiles)
# the -fPIC flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# the -flto

# adding project modules dir (with .cmake extension)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# including cmake modules
include(project_instructions)
# adding commander3 source files
# Note: I have added commander.cmake as an "independent" project 
# so there is no need for other cmake lists. Advantage of this 
# approach is that everything located in the same folder (cmake)
#add_subdirectory(commander3)
