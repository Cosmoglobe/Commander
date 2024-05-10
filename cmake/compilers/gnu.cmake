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
# Description: This file sets the default flags for GNU compilers. The compiler
# names are: gfortran, gcc, g++, mpifort (or mpif90), mpicc, mpicxx.
#------------------------------------------------------------------------------
# GNU - 9.3 - 10.x needs different flags
# setting different flags for different version
#------------------------------------------------------------------------------
if(COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
		"-O3"  
		"-DNDEBUG"
		"-fno-strict-aliasing"
		"-march=native" 
		"-flto" 
		"-fopenmp"
		"-fbacktrace" 
		"-fexternal-blas"
		"-ffpe-trap=zero"
		"-fPIC"
		#"-C"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
		"-O0"
		"-g" 
		"-Wuninitialized" # catching uninitialized variables
		"-fopenmp" 
		"-fbacktrace" 
		"-fexternal-blas"
		"-fPIC"
    "-fallow-argument-mismatch"
    #"-C" 
    #"-fno-strict-aliasing"
    #"-Wall" 
    #"-Wextra" 
    #"-Warray-temporaries"
    #"-Wconversion-extra" 
    #"-pedantic" 
    #"-Wmaybe-uninitialized"
    #"-fcheck=all" 
    #"-ffpe-trap=invalid,zero,overflow,underflow" 
    #"-ffunction-sections" 
    #"-pipe"
    #"-ffpe-trap=zero"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
		"-O2" 
		"-g" 
		"-Wuninitialized" # catching uninitialized variables
		"-fopenmp" 
		"-fbacktrace" 
		"-fexternal-blas"
		"-fPIC"
		#"-fno-strict-aliasing"
		#"-DNDEBUG" 
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
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
		"-Os"# -fno-strict-aliasing -DNDEBUG -fopenmp -fbacktrace -C -fexternal-blas -ffpe-trap=zero -fPIC" 
		"-fno-strict-aliasing"
		"-DNDEBUG" 
		"-fopenmp" 
		"-fbacktrace" 
		"-C"
		"-fexternal-blas"
		"-ffpe-trap=zero"
		"-fPIC"
		)
endif()
# adding different flags depending on the compiler version
list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS 
	  #"-Wfatal-errors"
		"-ffree-line-length-none" 
		"-fno-range-check"
	)
if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL "10")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS 
		"-fallow-argument-mismatch"
		)
else()
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS 
		"-Wno-argument-mismatch"
		)
endif()

# Linker flags
# the same logic as with compiler flags
if(COMMANDER3_Fortran_LINKER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_LINKER_FLAGS_RELEASE "-flto")
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
