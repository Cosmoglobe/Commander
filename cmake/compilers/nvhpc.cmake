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
# Description: This file sets the default flags for NVIDIA compilers. The compiler
# names are: nvfortran, nvc, nvc++, mpifort, mpicc, mpicxx.
#------------------------------------------------------------------------------
if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
		"-O4"
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
		"-O0"
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
		"-O2"  
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
		"-O0"  
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
