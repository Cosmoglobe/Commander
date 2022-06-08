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
# Description: This script defines how to install OpenBLAS or AMD BLIS/FLAME from
# source on the host system. As a result, the `blas` target is created to be con-
# sumed by Commander3. If BLAS/LAPACK implementation was found by find_package(),
# then `blas` target is empty.
#================================================================================
if(COMPILE_OPENBLAS)
	#------------------------------------------------------------------------------
	# Note: the explicit splitting for download and install step is done on purpose
	# to avoid errors when you want to recompile libraries for different owls etc.
	# In addition, this will allow us to download sources only once and then just 
	# reuse it whenever possible.
	#------------------------------------------------------------------------------
	# Getting OpenBLAS from source
	#------------------------------------------------------------------------------
	# Checking whether we have source directory and this directory is not empty.
  #
  if(NOT EXISTS "${OPENBLAS_SOURCE_DIR}/CMakeLists.txt")
    #message(STATUS "No BLAS sources were found; thus, will download it from source:\n${blas_url}")
		ExternalProject_Add(
			openblas_src
			DEPENDS						required_libraries
			URL								"${openblas_url}"
			URL_MD5						"${openblas_md5}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
      SOURCE_DIR				"${OPENBLAS_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
			LOG_DOWNLOAD			ON
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
			)
	else()
    #message(STATUS "Found an existing BLAS sources inside:\n${OPENBLAS_SOURCE_DIR}")
		add_custom_target(openblas_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
  # Compiling and Installing Static and Shared OpenBLAS
	#------------------------------------------------------------------------------
  list(APPEND _OPENBLAS_LIB_TYPES_ static shared)
  list(APPEND _OPENBLAS_LIB_BOOL_VALS_ -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON)
  foreach(_lib_type_ _bool_val_ IN ZIP_LISTS _OPENBLAS_LIB_TYPES_ _OPENBLAS_LIB_BOOL_VALS_)
    ExternalProject_Add(
			openblas_${_lib_type_}
      DEPENDS           required_libraries 
                        openblas_src
      PREFIX            "${LIBS_BUILD_DIR}"
      SOURCE_DIR        "${OPENBLAS_SOURCE_DIR}"
      INSTALL_DIR       "${CMAKE_INSTALL_PREFIX}"
      LOG_DIR	          "${CMAKE_LOG_DIR}"
      LOG_CONFIGURE     ON 
      LOG_BUILD         ON 
			LOG_INSTALL       ON 
      # commands how to build the project
      DOWNLOAD_COMMAND  ""
      CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_Fortran_COMPILER=${MPI_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${MPI_C_COMPILER}
        -DCMAKE_INSTALL_LIBDIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        -DUSE_OPENMP=OFF #<= HKE said we are not using OMP version of BLAS/LAPACK
        # To fix error Nils had with SMP_SERVER, this variable should be set to `false`
        # For me, however, it does absolutely nothing for some reason
        -DUSE_THREAD=OFF
        # Building both static and shared libraries to be consistent 
        # with Make and Autotools
				${_bool_val_}
      )
  endforeach()
	#------------------------------------------------------------------------------
	# In case of static linking, we do not need to specify linker flags.
	# NOTE: If OpenBLAS 0.3.12 is installed with CMake, there is a special
	# variable, which controls the installation of libraries and it is 
	# CMAKE_INSTALL_LIBDIR. We need to set it manually as above to "lib"
	# otherwise on some platforms the dir will be named "lib" and on others
	# "lib64".
	#------------------------------------------------------------------------------
	set(BLAS_LINKER_FLAGS "")
	set(BLAS_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	set(LAPACK_LINKER_FLAGS "")
	set(LAPACK_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	#------------------------------------------------------------------------------
  # Creating ALIASes
	#------------------------------------------------------------------------------
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
    )
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_OPTIONS "${BLAS_LINKER_FLAGS}"
    )
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
    )
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_OPTIONS "${LAPACK_LINKER_FLAGS}"
    )
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(blas   
    ALL ""
		DEPENDS openblas_static
						openblas_shared
    )
	#------------------------------------------------------------------------------
	#message(STATUS "BLAS LINKER FLAGS will be:   ${BLAS_LINKER_FLAGS}")
	#message(STATUS "BLAS LIBRARIES will be:      ${BLAS_LIBRARIES}")
	#message(STATUS "LAPACK LINKER FLAGS will be: ${LAPACK_LINKER_FLAGS}")
	#message(STATUS "LAPACK LIBRARIES will be:    ${LAPACK_LIBRARIES}")
	#------------------------------------------------------------------------------
elseif(COMPILE_FLAME)
	#------------------------------------------------------------------------------
  # Note: BLIS and FLAME have CMake implementation only for Windows.
	#------------------------------------------------------------------------------
  # Getting BLIS from source
	#------------------------------------------------------------------------------
  # Note: For BLIS configure to properly work it needs to find the 
  # blis.pc.in file located inside the root directory. Therefore, 
  # either the file should be copied to `${LIBS_BUILD_DIR}/src/blis-build` 
  # or the `BINARY_DIR` should have the same value as `SOURCE_DIR`. The 
  # latter, of course, undermines the purpose of `subbuilds` directory;
  # thus, the `copy` command is added to the `configure`. For future
  # references, the respective error is the following:
  # `gmake[4]: *** No rule to make target 'blis.pc.in', needed by`
	#------------------------------------------------------------------------------
  if(NOT EXISTS "${BLIS_SOURCE_DIR}/CMakeLists.txt")
    #message(STATUS "No BLAS sources were found; thus, will download it from source:\n${blas_url}")
		ExternalProject_Add(
			blis_src
			DEPENDS						required_libraries
      GIT_REPOSITORY		"${blis_git_url}"
      GIT_TAG						"${blis_git_tag}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
      SOURCE_DIR				"${BLIS_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
      LOG_DOWNLOAD			ON 
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
      COMMAND ${CMAKE_COMMAND} -E copy "${BLIS_SOURCE_DIR}/blis.pc.in" "${LIBS_BUILD_DIR}/src/blis-build"
			)
	else()
    #message(STATUS "Found an existing BLAS sources inside:\n${BLIS_SOURCE_DIR}")
		add_custom_target(blis_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
  # Getting FLAME from source
	#------------------------------------------------------------------------------
  # Note: LAPACK sources inside `src` directory is copied over to `subbuilds` , 
  # because, otherwise, it complains about missing files.
	#------------------------------------------------------------------------------
  if(NOT EXISTS "${FLAME_SOURCE_DIR}/CMakeLists.txt")
    #message(STATUS "No BLAS sources were found; thus, will download it from source:\n${blas_url}")
		ExternalProject_Add(
			flame_src
			DEPENDS						required_libraries
      GIT_REPOSITORY		"${flame_git_url}"
      GIT_TAG						"${flame_git_tag}"
			PREFIX						"${LIBS_BUILD_DIR}"
			DOWNLOAD_DIR			"${CMAKE_DOWNLOAD_DIRECTORY}"
      SOURCE_DIR				"${FLAME_SOURCE_DIR}"
			LOG_DIR						"${CMAKE_LOG_DIR}"
      LOG_DOWNLOAD			ON 
			# commands how to build the project
			CONFIGURE_COMMAND ""
			BUILD_COMMAND			""
			INSTALL_COMMAND		""
      COMMAND ${CMAKE_COMMAND} -E copy_directory "${FLAME_SOURCE_DIR}/src" "${LIBS_BUILD_DIR}/src/flame-build/src"
			)
	else()
    #message(STATUS "Found an existing BLAS sources inside:\n${BLIS_SOURCE_DIR}")
		add_custom_target(flame_src
			ALL ""
			)
	endif()
	#------------------------------------------------------------------------------
  # Compiling and Installing Static and Shared BLIS
	#------------------------------------------------------------------------------
	list(APPEND 
    blis_configure_command 
		"${CMAKE_COMMAND}" "-E" "env" 
    "FC=${MPI_Fortran_COMPILER}" 
		"CXX=${MPI_CXX_COMPILER}" 
		"CPP=${COMMANDER3_CPP_COMPILER}" 
		"CC=${MPI_C_COMPILER}" 
		"MPICC=${MPI_C_COMPILER}" 
		#"./configure" 
    "${BLIS_SOURCE_DIR}/configure" 
    "--prefix=<INSTALL_DIR>"
    "auto"
    )
	ExternalProject_Add(
		blis      
		DEPENDS						blis_src
		PREFIX						"${LIBS_BUILD_DIR}"
    SOURCE_DIR				"${BLIS_SOURCE_DIR}"
    #BINARY_DIR				"${BLIS_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
    LOG_CONFIGURE			ON 
    LOG_BUILD					ON 
    LOG_INSTALL				ON 
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${blis_configure_command}"
		)
	#------------------------------------------------------------------------------
  # Compiling and Installing Static and Shared FLAME
	#------------------------------------------------------------------------------
  # Note: Usage of MPi compilers result in the errors:
  # `Unrecognised command line option ...`
  # So the usual GCC/GFortran is used instead. 
	#------------------------------------------------------------------------------
	list(APPEND 
    flame_configure_command 
    #"rm" "Makefile" "&&"
    #"autoreconf" "-i" "&&" 
		"${CMAKE_COMMAND}" "-E" "env" 
    "FC=${CMAKE_Fortran_COMPILER}"#${MPI_Fortran_COMPILER}" 
    "CXX=${CMAKE_CXX_COMPILER}"#${MPI_CXX_COMPILER}" 
    "CPP=${COMMANDER3_CPP_COMPILER}" 
    "CC=${CMAKE_C_COMPILER}"#${MPI_C_COMPILER}" 
    #"MPICC=${MPI_C_COMPILER}" 
    #"./configure"
    "${FLAME_SOURCE_DIR}/configure" 
    "--prefix=<INSTALL_DIR>"
    "--enable-static-build" 
    "--enable-dynamic-build" 
    "--enable-blas-ext-gemmt" 
    "--enable-optimizations" 
    "--enable-warnings" 
    "--enable-max-arg-list-hack" 
    "--enable-lapack2flame"
    )
	ExternalProject_Add(
		flame
		DEPENDS						flame_src
                      blis 
		PREFIX						"${LIBS_BUILD_DIR}"
    SOURCE_DIR				"${FLAME_SOURCE_DIR}"
    #BINARY_DIR				"${BLIS_SOURCE_DIR}" 
		INSTALL_DIR				"${CMAKE_INSTALL_PREFIX}"
		LOG_DIR						"${CMAKE_LOG_DIR}"
    LOG_CONFIGURE			ON 
    LOG_BUILD					ON 
    LOG_INSTALL				ON 
		# Disabling download
		DOWNLOAD_COMMAND	""
		# Commands to configure, build and install the project
		CONFIGURE_COMMAND "${flame_configure_command}"
    # Note: For some reason CMake default install command results in errors,
    # but if `make install` is explicitly stated, everything works fine.
    INSTALL_COMMAND   "make" "install"
		)
	set(BLAS_LINKER_FLAGS "")
	set(BLAS_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}blis${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	set(LAPACK_LINKER_FLAGS "")
	set(LAPACK_LIBRARIES
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}flame${CMAKE_STATIC_LIBRARY_SUFFIX}"
		"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}aocldtl${CMAKE_STATIC_LIBRARY_SUFFIX}"
		)
	#------------------------------------------------------------------------------
  # Creating ALIASes
	#------------------------------------------------------------------------------
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
    )
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_OPTIONS "${BLAS_LINKER_FLAGS}"
    )
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
    )
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_OPTIONS "${LAPACK_LINKER_FLAGS}"
    )
	#------------------------------------------------------------------------------
	# Creating Unified Target
	#------------------------------------------------------------------------------
	add_custom_target(blas 
		ALL ""
		DEPENDS blis
            flame
		)
else()
	# to avoid cmake errors we create an empty target
	add_custom_target(blas 
		ALL ""
		DEPENDS required_libraries
		)
#	set(blas_lib 
#    ${BLAS_LINKER_FLAGS} 
#    ${BLAS_LIBRARIES} 
#    ${LAPACK_LINKER_FLAGS} 
#    ${LAPACK_LIBRARIES}
#    )
  #------------------------------------------------------------------------------
	#message(STATUS "BLAS LINKER FLAGS:   ${BLAS_LINKER_FLAGS}")
	#message(STATUS "BLAS LIBRARIES:      ${BLAS_LIBRARIES}")
	#message(STATUS "LAPACK LINKER FLAGS: ${LAPACK_LINKER_FLAGS}")
	#message(STATUS "LAPACK LIBRARIES:    ${LAPACK_LIBRARIES}")
	#------------------------------------------------------------------------------
endif()
