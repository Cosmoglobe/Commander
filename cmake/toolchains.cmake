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
# This file contains general instructions how to
# fetch and build the Commander dependencies
#==============================================================================


#==============================================================================
# COMPILE AND LINK OPTIONS 
#==============================================================================
# Setting project configuration -- variable is empty by default
# - "Debug" builds library/executable w/o optimization and w/ debug symbols; 
# - "Release" builds library/executable w/ optimization and w/o debug symbols;
# - "RelWithDebInfo" builds library/executable w/ less aggressive optimizations and w/ debug symbols;
# - "MinSizeRel" builds library/executable w/ optimizations that do not increase object code size. 
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo
		CACHE STRING
		"Specifies the Build type. Available options are: Release, Debug, RelWithDebInfo, MinSizeRel. Default: Release." FORCE)
endif()
#------------------------------------------------------------------------------
# Currently supported Compilers (CMake v3.21):
# https://cmake.org/cmake/help/v3.21/variable/CMAKE_LANG_COMPILER_ID.html
#
# - Absoft = Absoft Fortran (absoft.com)
# - ADSP = Analog VisualDSP++ (analog.com)
# - AppleClang = Apple Clang (apple.com)
# - ARMCC = ARM Compiler (arm.com)
# - ARMClang = ARM Compiler based on Clang (arm.com)
# - Bruce = Bruce C Compiler
# - CCur = Concurrent Fortran (ccur.com)
# - Clang = LLVM Clang (clang.llvm.org)
# - Cray = Cray Compiler (cray.com)
# - Embarcadero, Borland = Embarcadero (embarcadero.com)
# - Flang = Flang LLVM Fortran Compiler
# - Fujitsu = Fujitsu HPC compiler (Trad mode)
# - FujitsuClang = Fujitsu HPC compiler (Clang mode)
# - G95 = G95 Fortran (g95.org)
# - GNU = GNU Compiler Collection (gcc.gnu.org)
# - GHS = Green Hills Software (www.ghs.com)
# - HP = Hewlett-Packard Compiler (hp.com)
# - IAR = IAR Systems (iar.com)
# - Intel = Intel Compiler (intel.com)
# - IntelLLVM = Intel LLVM-Based Compiler (intel.com)
# - MSVC = Microsoft Visual Studio (microsoft.com)
# - NVHPC = NVIDIA HPC SDK Compiler (nvidia.com)
# - NVIDIA = NVIDIA CUDA Compiler (nvidia.com)
# - OpenWatcom = Open Watcom (openwatcom.org)
# - PGI = The Portland Group (pgroup.com)
# - PathScale = PathScale (pathscale.com)
# - ROCMClang = ROCm Toolkit Clang-based Compiler (rocmdocs.amd.com)
# - SDCC = Small Device C Compiler (sdcc.sourceforge.net)
# - SunPro = Oracle Solaris Studio (oracle.com)
# - TI = Texas Instruments (ti.com)
# - TinyCC = Tiny C Compiler (tinycc.org)
# - XL, VisualAge, zOS = IBM XL (ibm.com)
# - XLClang = IBM Clang-based XL (ibm.com)
#------------------------------------------------------------------------------
# Setting custom compile/link flags for each release type. These are
# additional flags, which user can define if he/she is not satisfied
# with the existing ones.
#------------------------------------------------------------------------------
set(COMMANDER3_Fortran_COMPILER_FLAGS "" #"${CMAKE_Fortran_FLAGS}"
	CACHE STRING
	"List of all additional flags user wants to add to configuration."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS ""
	CACHE STRING
	"List of additional linker flags user wants to add to configuration."
	)
# setting default compiler flags, but user will be able to overwrite them
# (although it is not really recommended to do so).
set(COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE ""
	CACHE STRING
	"List of compiler flags for Release version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG ""
	CACHE STRING
	"List of compiler flags for Debug version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO ""
	CACHE STRING
	"List of compiler flags for RelWithDebInfo version."
	)
set(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL ""
	CACHE STRING
	"List of compiler flags for MinSizeRel version."
	)

# The same as the above but linker flags
set(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE ""
	CACHE STRING
	"List of linker flags for Release version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG ""
	CACHE STRING
	"List of linker flags for Debug version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO ""
	CACHE STRING
	"List of linker flags for RelWithDebInfo version."
	)
set(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL ""
	CACHE STRING
	"List of linker flags for MinSizeRel version."
	)
# Commander also has one cpp file, so I add one flag
# similarly to the ones specified in the config files
set(COMMANDER3_CXX_COMPILER_FLAGS "-O3")

#------------------------------------------------------------------------------
# Specific flags for each subproject
# According to installation guidelines, we need to 
# specify slightly different flags for different compilers
# For easier debugging, I have created compiler variables
#------------------------------------------------------------------------------
set(COMMANDER3_C_COMPILER "${MPI_C_COMPILER}")
set(COMMANDER3_CXX_COMPILER "${MPI_CXX_COMPILER}")
set(COMMANDER3_Fortran_COMPILER "${MPI_Fortran_COMPILER}")

#------------------------------------------------------------------------------
# To get rid of "C preprocessor fails sanity check" error
# we need to link the CPP as follows for each subproject
# Different OS requires different linking
#------------------------------------------------------------------------------
if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
endif()

#------------------------------------------------------------------------------
# Specifying flags per Fortran compiler
# Intel
# Note: from https://www.osc.edu/documentation/knowledge_base/compilation_guide
# With the Intel compilers, use -xHost and -O2 or higher. 
# With the gnu compilers, use -march=native and -O3. 
# The PGI compilers by default use the highest available instruction set, so no additional flags are necessary.
#------------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	include(intel)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	include(gnu)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
	include(nvhpc)
#------------------------------------------------------------------------------
# AOCC Flang
#------------------------------------------------------------------------------
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
	if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
			"-O0" 
			"-g" 
			#"-fno-strict-aliasing"
			#"-fopenmp" 
			#"-fbacktrace" 
			#"-fexternal-blas"
			#"-C" 
			#"-Wall" 
			#"-Wextra" 
			#"-Warray-temporaries"
			#"-Wconversion-extra" 
			#"-pedantic" 
			#"-fcheck=all" 
			#"-ffpe-trap=invalid,zero,overflow,underflow" 
			#"-ffunction-sections" 
			#"-pipe"
			#"-ffpe-trap=zero"
			#"-fPIC"
			)
	endif()
#------------------------------------------------------------------------------
# PGI	
# TODO: For some reason commander dependencies crashes on this one
# so figure out how to make it work
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	# Compiler flags
	if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
			"-O4"# -fast -mp=all -traceback -Mconcur -fPIC" 
			"-DNDEBUG"
			"-fast" 
			"-mp=all"
			"-traceback" 
			"-Mconcur"
			"-fPIC"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
			"-O0" # -mp=all -gopt -fast -traceback -Minfo -Mconcur -C -fPIC" 
			"-mp=all"
			"-gopt" 
			"-fast" 
			"-traceback" 
			"-Minfo" 
			"-Mconcur"
			"-C"
			"-fPIC"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
			"-O2"# -mp=all -gopt -fast -traceback -Minfo -Mconcur -C -fPIC" 
			"-DNDEBUG"
			"-mp=all"
			"-gopt" 
			"-fast" 
			"-traceback" 
			"-Minfo" 
			"-Mconcur"
			"-C"
			"-fPIC"
			)
	endif()
	if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
			"-O0"# -mp=all -fast -traceback -Mconcur -fPIC" 
			"-mp=all"
			"-fast" 
			"-traceback" 
			"-Mconcur"
			"-C"
			"-fPIC"
			)
	endif()

	# Linker flags
	# the same logic as with compiler flags
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE 
			"-mp=all"# -gopt -Mconcur"
			"-gopt" 
			"-Mconcur"
			)
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG 
			"-mp=all"# -gopt -Mconcur"
			"-gopt" 
			"-Mconcur"
			)
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO 
			"-mp=all"# -Mconcur"
			"-Mconcur"
			)
	endif()
	if(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL MATCHES "")
		list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL 
			"-mp=all"# -Mconcur"
			"-Mconcur"
			)
	endif()
#------------------------------------------------------------------------------
# Flang
# TODO: need to figure out why healpix doesn't compile with flang
# and then add support for flang
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
#	# Compiler flags
#	#list(APPEND COMMANDER3_COMPILER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL "")
#	# Linker flags
#	#list(APPEND COMMANDER3_LINKER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
#------------------------------------------------------------------------------
# Flang
# TODO: need to figure out why healpix doesn't compile with flang
# and then add support for flang
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES Flang)
#	# Compiler flags
#	#list(APPEND COMMANDER3_COMPILER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO "")
#	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL "")
#	# Linker flags
#	#list(APPEND COMMANDER3_LINKER_FLAGS "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
#	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
endif()
#------------------------------------------------------------------------------
# Making a summary of compiler location and compile flags
#------------------------------------------------------------------------------
message(STATUS "---------------------------------------------------------------")
message(STATUS "SUMMARY ON COMPILERS:")
message(STATUS "---------------------------------------------------------------")
message(STATUS "Your system is: ${CMAKE_SYSTEM_NAME}, ${CMAKE_SYSTEM}")
message(STATUS "Fortran Compiler is: ${CMAKE_Fortran_COMPILER}")
message(STATUS "C Compiler is: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler is: ${CMAKE_CXX_COMPILER}")
message(STATUS "Commander3 configuration is: ${CMAKE_BUILD_TYPE}. Compiler flags to be applied:")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG} ${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_Fortran_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE} ${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO} ${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
	message(STATUS "${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL} ${COMMANDER3_Fortran_COMPILER_FLAGS};")#${COMMANDER3_COMPILER_FLAGS_ADDITIONAL}")
endif()
#------------------------------------------------------------------------------
# Making a summary of Host System 
#------------------------------------------------------------------------------
# CMake reference:
# https://cmake.org/cmake/help/v3.17/module/ProcessorCount.html 
# https://cmake.org/cmake/help/v3.17/command/cmake_host_system_information.html
include(ProcessorCount)
ProcessorCount(N_CORES)

cmake_host_system_information(RESULT N_LOGICAL_CORES  QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT N_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
# Processor
cmake_host_system_information(RESULT PROC_NAME QUERY PROCESSOR_NAME)
cmake_host_system_information(RESULT PROC_DESCRIPTION QUERY PROCESSOR_DESCRIPTION)
# OS information
cmake_host_system_information(RESULT HOST_OS_NAME QUERY OS_NAME)
cmake_host_system_information(RESULT HOST_OS_RELEASE QUERY OS_RELEASE)
cmake_host_system_information(RESULT HOST_OS_VERSION QUERY OS_VERSION)
cmake_host_system_information(RESULT HOST_OS_PLATFORM QUERY OS_PLATFORM)

message(STATUS ${HOST_NAME})
message(STATUS "${HOST_OS_NAME} ${HOST_OS_PLATFORM} ${HOST_OS_RELEASE} ${HOST_OS_VERSION}")
message(STATUS "${PROC_NAME} | ${PROC_DESCRIPTION}")

