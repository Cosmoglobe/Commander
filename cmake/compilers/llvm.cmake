
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
# If user has not specified compilation flag, we use default configuration
if (COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE 
		"-O3"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG 
		"-O0"  
		"-g" 
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO 
		"-O2"
		"-xHost" 
		"-fpe0"
		"-fPIC"
		"-fp-model" "strict"
		#"-qopenmp" 
		"-assume" "byterecl" # for I/O operations 
		#"-qopt-matmul" #<= increases linking time but doesn't increase performance 
		"-g" 
		"-traceback" 
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
		#"-heap-arrays" "16384"
		#"-fpe0"
		#"-fPIC"
		)
endif()
if(COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL MATCHES "")
	list(APPEND COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL 
		"-Os"# -traceback -DNDEBUG -parallel -qopenmp -assume byterecl -heap-arrays 16384 -fpe0 -fPIC" 
		"-traceback" 
		"-DNDEBUG" 
		#"-parallel" 
		#"-qopenmp" 
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
