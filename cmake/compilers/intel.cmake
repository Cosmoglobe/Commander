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
# Description: This file sets the default flags for Intel compilers. The compiler
# names are: ifort, icc, icpc, mpiifort, mpiicc and mpiicpc. 
# 
# Intel reference guide:
# https://software.intel.com/content/dam/develop/public/us/en/documents/quick-reference-guide-intel-compilers-v19-1-final-.pdf
# https://software.intel.com/content/www/us/en/develop/download/quick-reference-guide-to-optimization-with-intel-compilers-v19.html
# Useful resources:
# https://www.bu.edu/tech/support/research/software-and-programming/programming/compilers/intel-compiler-flags/
#
# Type in the command line:
# $ man ifort
# To see the Fortran reference guide

# the `man ifort` for Intel 2020 gives us:
#
# Options that Improve Run-Time Performance
#		The  following  command  line  options can be used to increase the run-time performance 
#   of code generated by the Intel(R) Fortran Compiler:
#
#		On systems  using  IA-32  architecture  and  systems  using  Intel(R)  64  architecture:
#		-ax<processor>, -ftz, -ip, -ipo, -march=<processor>, -mtune=<processor>, -O[n], -openmp,
#		-parallel, -prof-gen, -prof-use, -x<processor>.
#
# -Ofast -- sets compiler options -O3, -no-prec-div, and -fp-model fast=2.
# -xHost -- 
# -ip    -- Single file interprocedural optimizations, including selective inlining, within the current source file. 
# -fast  -- maximizes speed across the entire program. Sets the following options:
#		On macOS* systems: -ipo, -mdynamic-no-pic,-O3, -no-prec-div,-fp-model fast=2, and -xHost
#		On Windows* systems: /O3, /Qipo, /Qprec-div-, /fp:fast=2, and /QxHost
#		On Linux* systems: -ipo, -O3, -no-prec-div,-static, -fp-model fast=2, and -xHost
# -q[no-]opt-matmul -- Enables or disables a compiler-generated Matrix Multiply (matmul) library call. This option is enabled by  default  at  setting  O2  and  above.
# --fp-model=keyword:
# Other options are here:
# https://software.intel.com/content/www/us/en/develop/documentation/oneapi-dpcpp-cpp-compiler-dev-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/floating-point-options/fp-model-fp.html
#------------------------------------------------------------------------------
# Specifying flags per Fortran compiler
# Intel
# Note: from https://www.osc.edu/documentation/knowledge_base/compilation_guide
# With the Intel compilers, use -xHost and -O2 or higher. 
# With the GNU compilers, use -march=native and -O3. 
# The PGI compilers by default use the highest available instruction set, so no additional flags are necessary.
#------------------------------------------------------------------------------
# Compiler flags
# If user has not specified compilation flag, we use default configuration
if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
		"-O3"
		"-xHost" 
		"-fpe0"
		"-fPIC"
		#"-fp-model=strict"
		"-traceback" 
		"-qopenmp" #<= we are not using it at all, it is redundant 
		"-assume" "byterecl" # for I/O operations 
		#"-qopt-matmul" #<= increases linking time but doesn't increase performance 
		#"-DNDEBUG"
		#"-ipo" #  
		#"-parallel" 
		"-heap-arrays" "16384"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
		"-O0"  
		"-g" 
		"-xHost" 
		"-debug" "all"
		"-check" "all,noarg_temp_created"
		#"-warn" "all"
		"-fp-stack-check"
		"-fstack-protector-all"
		"-traceback" 
		#"-parallel" 
		"-qopenmp"
		"-C" 
		"-assume" "byterecl" 
		"-heap-arrays" "16384"
		"-fpe0"
		"-fPIC"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
		"-O2"
		"-xHost" 
		"-g" 
		"-fpe0"
		"-fPIC"
		#"-fp-model=strict"
		"-qopenmp" 
		"-assume" "byterecl" # for I/O operations 
		#"-qopt-matmul" #<= increases linking time but doesn't increase performance 
		"-traceback" 
		"-heap-arrays" "16384"
		#
		#"-O2"  
		#"-g" 
		#"-traceback" 
		#"-DNDEBUG" 
		#"-parallel" 
		#"-qopenmp"
		#"-qopt-matmul"
		#"-C"
		#"-assume" "byterecl" 
		#"-fpe0"
		#"-fPIC"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
		"-Os"
		"-traceback" 
		"-DNDEBUG" 
		"-parallel" 
		"-qopenmp" 
		"-C"
		"-assume" "byterecl" 
		"-heap-arrays" "16384"
		"-fpe0"
		"-fPIC"
		)
endif()
# Linker flags
# the same logic as with compiler flags
if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "-qopt-matmul")
endif()
if(COMMANDER3_Fortran_LINKER_FLAGS_DEBUG MATCHES "")
	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_DEBUG "")
endif()
if(COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO MATCHES "")
	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO "")
endif()
if(COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL MATCHES "")
	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL "")
endif()