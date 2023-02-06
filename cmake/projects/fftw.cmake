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
# Description: This script determines the location of FFTW on the host system.
# If it fails to do so, it will download, compile and install FFTW from source.
#================================================================================

# TODO: I am using FFTW3_* instead of FFTW_* variables (because MPI is treated 
# separetely), that is why in the ned I have FFTW3_LIBRARIES instead of 
# FFTW_LIBRARIES. Perhaps this is incorrect? Need to take a look into FindFFTW.cmake
# and rewrite it if necessary.
#message(STATUS "---------------------------------------------------------------")
#if(USE_SYSTEM_FFTW AND USE_SYSTEM_LIBS)
#	find_package(FFTW 
#		COMPONENTS 
#		DOUBLE 
#		DOUBLE_THREADS
#		FLOAT 
#		#FLOAT_MPI 
#		FLOAT_OPENMP
#		FLOAT_THREADS
#		)
#endif()
# Is TRUE if one of the components were missing
if(COMPILE_FFTW)
	#------------------------------------------------------------------------------
	# Splitting the project into 5 steps:
	# 1. To download the project
	# 2. To compile with double static 
	# 3. To compile with double shared 
	# 4. To compile with float static 
	# 5. To compile with float shared 
	#
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting FFTW from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
	if(NOT EXISTS "${FFTW_SOURCE_DIR}/CMakeLists.txt")
    #message(STATUS "No FFTW sources were found; thus, will download it from source:\n${fftw_url}")
		ExternalProject_Add(
			fftw_src
			URL								"${fftw_url}"
			URL_MD5						"${fftw_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
			SOURCE_DIR				"${FFTW_SOURCE_DIR}"
			BINARY_DIR				"${FFTW_SOURCE_DIR}" 
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# Ommiting Configuration, build and install steps
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
    #message(STATUS "Found an existing FFTW sources inside:\n${FFTW_SOURCE_DIR}")
		add_custom_target(fftw_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
	# Compiling and Installing Static and Shared FFTW
	#------------------------------------------------------------------------------
	# Looping over libraries we need to compile
	list(APPEND _FFTW_NAMES_ double float)
	list(APPEND _FFTW_ARGS_ -DENABLE_FLOAT:BOOL=OFF -DENABLE_FLOAT:BOOL=ON)
	list(APPEND _FFTW_LIB_TYPE_ shared static)
	list(APPEND _FFTW_LIB_BOOL_VAL_ -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF)
	foreach(_fftw_component_ _fftw_arg_ IN ZIP_LISTS _FFTW_NAMES_ _FFTW_ARGS_)
		foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _FFTW_LIB_TYPE_ _FFTW_LIB_BOOL_VAL_)
			ExternalProject_Add(
				fftw_${_fftw_component_}_${_lib_type_}
				DEPENDS           fftw_src
				PREFIX            "${LIBS_BUILD_DIR}"
				SOURCE_DIR        "${FFTW_SOURCE_DIR}"
				#BINARY_DIR        "${FFTW_SOURCE_DIR}"
				INSTALL_DIR       "${CMAKE_INSTALL_PREFIX}"
				LOG_DIR           "${CMAKE_LOG_DIR}"
        LOG_CONFIGURE     ON 
        LOG_BUILD         ON 
        LOG_INSTALL       ON 
				# Disabling download
				DOWNLOAD_COMMAND  ""
				CMAKE_ARGS
          -DCMAKE_BUILD_TYPE:STRING=Release
					# Specifying installations paths for binaries and libraries
					-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
					# Specifying compilers
          #-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
					-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
					-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
					# Building both static and shared libraries
					${_bool_val_}
					# Which libraries to produce
					-DENABLE_OPENMP:BOOL=ON
					-DENABLE_THREADS:BOOL=ON
          -DENABLE_SSE:BOOL=${FFTW_ENABLE_SSE}
          -DENABLE_SSE2:BOOL=${FFTW_ENABLE_SSE2}
					-DENABLE_AVX:BOOL=${FFTW_ENABLE_AVX}
					-DENABLE_AVX2:BOOL=${FFTW_ENABLE_AVX2}
					${_fftw_arg_}
					# ensuring it will be installed inside `lib` and not `lib64`
					-DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
				)
		endforeach()
	endforeach()
	#------------------------------------------------------------------------------
	# Adding fftw3, fftw3_threads, and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
	set(FFTW_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"	
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		)
	set(FFTW_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include"	
		)
	include_directories(${FFTW_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "FFTW LIBRARIES will be: ${FFTW_LIBRARIES}")
	#message(STATUS "FFTW INCLUDE DIRS will be: ${FFTW_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(fftw 
		ALL ""
		DEPENDS fftw_float_static
						fftw_double_static
		        fftw_float_shared
						fftw_double_shared
		)
	#------------------------------------------------------------------------------
elseif(COMPILE_AMDFFTW)
  # TODO: Add compilation of FFTW from AMD
	#------------------------------------------------------------------------------
  # Getting AMD FFTW from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
  if(NOT EXISTS "${AMDFFTW_SOURCE_DIR}/CMakeLists.txt")
    #message(STATUS "No FFTW sources were found; thus, will download it from source:\n${fftw_url}")
		ExternalProject_Add(
			amdfftw_src
      GIT_REPOSITORY		"${amdfftw_git_url}"
      GIT_TAG						"${amdfftw_git_tag}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
      SOURCE_DIR				"${AMDFFTW_SOURCE_DIR}"
      BINARY_DIR				"${AMDFFTW_SOURCE_DIR}" 
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# Ommiting Configuration, build and install steps
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
    #message(STATUS "Found an existing FFTW sources inside:\n${FFTW_SOURCE_DIR}")
		add_custom_target(amdfftw_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
  # Compiling and Installing Static and Shared AMD FFTW
	#------------------------------------------------------------------------------
  # Note: There is some problem with enabling SSE support for AMD FFTW. Even the 
  # official configure command omits it:
  # ./configure --enable-sse2 --enable-avx --enable-avx2 --enable-mpi 
  # --enable-openmp --enable-shared --enable-amd-opt --enable-amd-mpifft 
  # --prefix=<your-install-dir>
  # So it will be ignored here as well.
  # 
  # In addition, there is an error when compiling tests (via CMake) with OpenMP 
  # support (cannot find symbols for omp_num_threads etc.), so the tests are 
  # disabled. The issue is resolved by adding `-fopenmp` flag while compiling and 
  # linking, which is not the case for the tests.
  #
  # Also, when compiling with configure, it requires to have both non-omp and omp 
  # versions of the libraries, whereas CMake version doe snto complaine for some 
  # reason.
	#------------------------------------------------------------------------------
  list(APPEND amdfftw_double_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
    #"FC=${MPI_Fortran_COMPILER}" 
    "CXX=${CMAKE_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
    "CC=${CMAKE_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
    "${AMDFFTW_SOURCE_DIR}/configure" 
    "--enable-shared"
    "--enable-openmp"
    "--enable-amd-opt"
    "--enable-sse2" 
    "--enable-avx" 
    "--enable-avx2"
    "--enable-amd-mpifft"
		"--prefix=<INSTALL_DIR>"
    )
	ExternalProject_Add(
		amdfftw_double
		DEPENDS				    amdfftw_src
    #amdfftw_float
		PREFIX						"${LIBS_BUILD_DIR}"
    SOURCE_DIR				"${AMDFFTW_SOURCE_DIR}"
    BINARY_DIR				"${AMDFFTW_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
    LOG_CONFIGURE			ON 
		LOG_BUILD					ON 
		LOG_INSTALL				ON
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${amdfftw_double_configure_command}"
		)
	list(APPEND amdfftw_float_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
    #"FC=${MPI_Fortran_COMPILER}" 
    "CXX=${CMAKE_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
    "CC=${CMAKE_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
    "${AMDFFTW_SOURCE_DIR}/configure" 
    "--enable-shared"
    "--enable-openmp"
    "--enable-float"
    "--enable-amd-opt"
    "--enable-sse2" 
    "--enable-avx" 
    "--enable-avx2"
    "--enable-amd-mpifft"
		"--prefix=<INSTALL_DIR>"
    )
	ExternalProject_Add(
		amdfftw_float 
		DEPENDS				    amdfftw_src
                      amdfftw_double
		PREFIX						"${LIBS_BUILD_DIR}"
    SOURCE_DIR				"${AMDFFTW_SOURCE_DIR}"
    BINARY_DIR				"${AMDFFTW_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
    LOG_CONFIGURE			ON 
		LOG_BUILD					ON 
		LOG_INSTALL				ON
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${amdfftw_float_configure_command}"
		)
  
	# Looping over libraries we need to compile
  # Precision (double, float)
  #list(APPEND _AMDFFTW_NAMES_ double float)
  #list(APPEND _AMDFFTW_PRECISION_ARGS_ -DENABLE_FLOAT:BOOL=OFF -DENABLE_FLOAT:BOOL=ON)
  ## with OpenMP or not
  #list(APPEND _AMDFFTW_NAMES_OMP_ omp threads noomp)
  #list(APPEND _AMDFFTW_OMP_ARGS_ -DENABLE_OPENMP:BOOL=ON -DENABLE_THREADS:BOOL=ON -DENABLE_OPENMP:BOOL=OFF)
  ##list(APPEND _AMDFFTW_OMP_ARGS_ -DENABLE_THREADS:BOOL=ON -DENABLE_THREADS:BOOL=OFF)
  ## static or shared
  #list(APPEND _AMDFFTW_LIB_TYPE_ shared static)
  #list(APPEND _AMDFFTW_LIB_BOOL_VAL_ -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF)
  #foreach(_component_ _precision_arg_ IN ZIP_LISTS _AMDFFTW_NAMES_ _AMDFFTW_PRECISION_ARGS_)
  #  foreach(_omp_ _omp_arg_ IN ZIP_LISTS _AMDFFTW_NAMES_OMP_ _AMDFFTW_OMP_ARGS_)
  #    foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _AMDFFTW_LIB_TYPE_ _AMDFFTW_LIB_BOOL_VAL_)
  #      ExternalProject_Add(
  #        amdfftw_${_component_}_${_omp_}_${_lib_type_}
  #        DEPENDS           amdfftw_src
  #        PREFIX            "${LIBS_BUILD_DIR}"
  #        SOURCE_DIR        "${AMDFFTW_SOURCE_DIR}"
  #        #BINARY_DIR        "${FFTW_SOURCE_DIR}"
  #        INSTALL_DIR       "${CMAKE_INSTALL_PREFIX}"
  #        LOG_DIR           "${CMAKE_LOG_DIR}"
  #        LOG_CONFIGURE     ON 
  #        LOG_BUILD         ON 
  #        LOG_INSTALL       ON 
  #        # Disabling download
  #        DOWNLOAD_COMMAND  ""
  #        CMAKE_ARGS
  #          -DCMAKE_BUILD_TYPE:STRING=Release
  #          # Specifying installations paths for binaries and libraries
  #          -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  #          # Specifying compilers
  #          #-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  #          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  #          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  #          # Building both static and shared libraries
  #          ${_bool_val_}
  #          # Which libraries to produce (only libfftw3_omp)
  #          #-DENABLE_THREADS:BOOL=ON
  #          #-DENABLE_OPENMP:BOOL=ON
  #          ${_omp_arg_}
  #          #-DENABLE_SSE:BOOL=${FFTW_ENABLE_SSE}
  #          -DENABLE_SSE2:BOOL=${FFTW_ENABLE_SSE2}
  #          -DENABLE_AVX:BOOL=${FFTW_ENABLE_AVX}
  #          -DENABLE_AVX2:BOOL=${FFTW_ENABLE_AVX2}
  #          # AMD optimizations
  #          -DENABLE_AMD_OPT:BOOL=ON
  #          -DENABLE_AMD_APP_OPT:BOOL=ON      # HPC optimizations, supported: float, double 
  #          -DENABLE_AMD_FAST_PLANNER:BOOL=ON # supported in float and double
  #          -DENABLE_AMD_TRANS:BOOL=ON
  #          -DENABLE_MPI:BOOL=OFF
  #          -DENABLE_AMD_MPIFFT:BOOL=ON 
  #          ${_precision_arg_}
  #          # ensuring it will be installed inside `lib` and not `lib64`
  #          -DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  #          # This gives an error, described above
  #          -DBUILD_TESTS:BOOL=OFF
  #        )
  #    endforeach()
	#	endforeach()
	#endforeach()
	#------------------------------------------------------------------------------
	# Adding fftw3, fftw3_threads, and fftws3_omp into a library variable
	# Defining this variable just to not to overwrite FFTW_LIBRARIES created by FindFFTW
  # Note: if we compile with configure, we need to include both non-omp and omp versions 
  # of the libraries. With CMake, however, we should use only _omp one.
	set(FFTW_LIBRARIES
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}"		
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX}"		
    #"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"	
    #"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_threads${CMAKE_STATIC_LIBRARY_SUFFIX}"
    #
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX}"		
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f_omp${CMAKE_STATIC_LIBRARY_SUFFIX}"		
		)
	set(FFTW_INCLUDE_DIRS
		"${CMAKE_INSTALL_PREFIX}/include"	
		)
	include_directories(${FFTW_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
  #add_custom_target(amdfftw_float
  #  ALL ""
  #  DEPENDS amdfftw_float_omp_static
  #          amdfftw_float_threads_static
  #          amdfftw_float_noomp_static
  #          amdfftw_float_omp_shared
  #          amdfftw_float_threads_shared
  #          amdfftw_float_noomp_shared
  #  )
  #add_custom_target(amdfftw_double
  #  ALL ""
  #  DEPENDS amdfftw_double_omp_static
  #          amdfftw_double_threads_static
  #          amdfftw_double_noomp_static
  #          amdfftw_double_omp_shared
  #          amdfftw_double_threads_shared
  #          amdfftw_double_noomp_shared
  #  )
  add_custom_target(fftw 
		ALL ""
    DEPENDS amdfftw_float       
		        amdfftw_double      
	#	DEPENDS amdfftw_float_static
	#					amdfftw_double_static
	#	        amdfftw_float_shared
	#					amdfftw_double_shared
		)
	#------------------------------------------------------------------------------
  
else()
	# adding empty targets in case FFTW was found on the system
	add_custom_target(fftw ALL "")
	set(FFTW_LIBRARIES
		${FFTW_DOUBLE_LIB}
		${FFTW_DOUBLE_THREADS_LIB}
		${FFTW_FLOAT_LIB}
		${FFTW_FLOAT_OPENMP_LIB}
		${FFTW_FLOAT_THREADS_LIB}
		#${FFTW_FLOAT_MPI_LIB}
		)
	include_directories(${FFTW_INCLUDE_DIRS})
	#------------------------------------------------------------------------------
	#message(STATUS "FFTW LIBRARIES will be: ${FFTW_LIBRARIES}")
	#message(STATUS "FFTW INCLUDE DIRS will be: ${FFTW_INCLUDE_DIRS}")
	#------------------------------------------------------------------------------
endif()
