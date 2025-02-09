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
		"Specifies the Build type. Available options are: Release, Debug, RelWithDebInfo, MinSizeRel. Default: RelWithDebInfo." FORCE)
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
set(COMMANDER3_CXX_COMPILER_FLAGS "-O0")

set(COMMANDER3_C_COMPILER "${MPI_C_COMPILER}")
set(COMMANDER3_CXX_COMPILER "${MPI_CXX_COMPILER}")
set(COMMANDER3_Fortran_COMPILER "${MPI_Fortran_COMPILER}")

if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
	set(COMMANDER3_CPP_COMPILER "${MPI_CXX_COMPILER} -E")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
	include(intel)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
	include(gnu)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
	include(nvhpc)
endif()
