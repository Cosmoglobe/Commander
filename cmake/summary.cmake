#================================================================================
# Author: Maksym Brilenkov
#================================================================================
# Description: This file contains Summary Information Output to the user's screen
# regarding System and CMake configuration which will be used to compile the code
#================================================================================
# To make pretty output of libraries in terminal window
function(change_list_output input_list output_list)
	string(REPLACE ";" "\n   |   - " dummy_list "${input_list}")
	set(${output_list} "${dummy_list}" PARENT_SCOPE)
endfunction()
#------------------------------------------------------------------------------
# Making a summary of Host System 
#------------------------------------------------------------------------------
# CMake reference:
# https://cmake.org/cmake/help/v3.21/module/ProcessorCount.html 
# https://cmake.org/cmake/help/v3.21/command/cmake_host_system_information.html
include(ProcessorCount)
ProcessorCount(N_CORES)

# CPU 
cmake_host_system_information(RESULT CPU_NAME              QUERY PROCESSOR_NAME)
cmake_host_system_information(RESULT CPU_DESCRIPTION       QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT N_LOGICAL_CORES       QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT N_PHYSICAL_CORES      QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT CPU_HAS_SERIAL_NUMBER QUERY HAS_SERIAL_NUMBER)
cmake_host_system_information(RESULT CPU_SERIAL_NUMBER     QUERY PROCESSOR_SERIAL_NUMBER)
cmake_host_system_information(RESULT CPU_IS_64BIT          QUERY IS_64BIT)
cmake_host_system_information(RESULT CPU_HAS_FPU           QUERY HAS_FPU)
cmake_host_system_information(RESULT CPU_HAS_MMX           QUERY HAS_MMX)
# One if processor supports Ext. MMX instructions
cmake_host_system_information(RESULT CPU_HAS_MMX_PLUS      QUERY HAS_MMX_PLUS) 
cmake_host_system_information(RESULT CPU_HAS_SSE           QUERY HAS_SSE)
cmake_host_system_information(RESULT CPU_HAS_SSE2          QUERY HAS_SSE2)
cmake_host_system_information(RESULT CPU_HAS_SSE_FP        QUERY HAS_SSE_FP)
# One if processor supports SSE MMX instructions
cmake_host_system_information(RESULT CPU_HAS_SSE_MMX       QUERY HAS_SSE_MMX)
# OS information
#cmake_host_system_information(RESULT HOST_OS_NAME QUERY OS_NAME)
#cmake_host_system_information(RESULT HOST_NAME             QUERY HOSTNAME)
cmake_host_system_information(RESULT HOST_OS_RELEASE       QUERY OS_RELEASE)
cmake_host_system_information(RESULT HOST_OS_VERSION       QUERY OS_VERSION)
cmake_host_system_information(RESULT HOST_OS_PLATFORM      QUERY OS_PLATFORM)
# RAM information
cmake_host_system_information(RESULT TOT_VIRTUAL_MEMORY    QUERY TOTAL_VIRTUAL_MEMORY)
cmake_host_system_information(RESULT AVAIL_VIRTUAL_MEMORY  QUERY AVAILABLE_VIRTUAL_MEMORY)
cmake_host_system_information(RESULT TOT_PHYSICAL_MEMORY   QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT AVAIL_PHYSICAL_MEMORY QUERY AVAILABLE_PHYSICAL_MEMORY)

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(_comm3_flags_ "${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG} ${COMMANDER3_Fortran_COMPILER_FLAGS};")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
	set(_comm3_flags_ "${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE} ${COMMANDER3_Fortran_COMPILER_FLAGS};")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	set(_comm3_flags_ "${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO} ${COMMANDER3_Fortran_COMPILER_FLAGS};")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
	set(_comm3_flags_ "${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL} ${COMMANDER3_Fortran_COMPILER_FLAGS};")
endif()

message(STATUS "==============================================================="
				"\n   |                    CONFIGURATION SUMMARY                    |"
				"\n   ===============================================================")
message(STATUS "|                          HARDWARE                           |"
				"\n   ---------------------------------------------------------------")
message(STATUS "| Operating System:"
        "\n   |------------------"
        "\n   | * Name            -- ${CMAKE_HOST_SYSTEM_NAME}"
        "\n   | * Platform        -- ${HOST_OS_PLATFORM}"
				"\n   | * Release         -- ${HOST_OS_RELEASE}"
				"\n   | * Version         -- ${HOST_OS_VERSION}"
        "\n   |------------------")
message(STATUS "| CPU             :"
        "\n   |------------------"
				"\n   | * Name            -- ${CPU_NAME}"
				"\n   | * Description     -- ${CPU_DESCRIPTION}"
				"\n   | * Physical Cores  -- ${N_PHYSICAL_CORES}"
				"\n   | * Logical Cores   -- ${N_LOGICAL_CORES}"
				"\n   | * Is 64Bit        -- ${CPU_IS_64BIT}"
				"\n   | * Has FPU         -- ${CPU_HAS_FPU}"
				"\n   | * Has SSE         -- ${CPU_HAS_SSE}"
				"\n   | * Has SSE2        -- ${CPU_HAS_SSE2}"
				"\n   | * Has SSE FP      -- ${CPU_HAS_SSE_FP}"
				"\n   | * Has SSE MMX     -- ${CPU_HAS_SSE_MMX}"
				"\n   | * Has MMX         -- ${CPU_HAS_MMX}"
				"\n   | * Has Ext. MMX    -- ${CPU_HAS_MMX_PLUS}"
        "\n   |------------------")
message(STATUS "| Memory          :"
        "\n   |------------------"
				"\n   | * Virtual [total] -- ${TOT_VIRTUAL_MEMORY} [MiB]"
				"\n   | * Virtual [avail] -- ${AVAIL_VIRTUAL_MEMORY} [MiB]"
				"\n   | * Phys    [total] -- ${TOT_PHYSICAL_MEMORY} [MiB]"
				"\n   | * Phys    [avail] -- ${AVAIL_PHYSICAL_MEMORY} [MiB]"
				"\n   ---------------------------------------------------------------")
message(STATUS "|                         COMPILERS                           |"
				"\n   ---------------------------------------------------------------")
message(STATUS "| Fortran         :"
        "\n   |------------------"
				"\n   | * Name            -- ${CMAKE_Fortran_COMPILER}"
				"\n   | * ID              -- ${CMAKE_Fortran_COMPILER_ID}"
				"\n   | * Version         -- ${CMAKE_Fortran_COMPILER_VERSION}"
				"\n   | * MPI Name        -- ${MPI_Fortran_COMPILER}"
				"\n   | * MPI Exec.       -- ${MPIEXEC_EXECUTABLE}"
        "\n   |------------------")
message(STATUS "| C++             :"
        "\n   |------------------"
				"\n   | * Name            -- ${CMAKE_CXX_COMPILER}"
				"\n   | * ID              -- ${CMAKE_CXX_COMPILER_ID}"
				"\n   | * Version         -- ${CMAKE_CXX_COMPILER_VERSION}"
				"\n   | * MPI Name        -- ${MPI_CXX_COMPILER}"
				"\n   | * MPI Exec.       -- ${MPIEXEC_EXECUTABLE}"
        "\n   |------------------")
message(STATUS "| C               :"
        "\n   |------------------"
				"\n   | * Name            -- ${CMAKE_C_COMPILER}"
				"\n   | * ID              -- ${CMAKE_C_COMPILER_ID}"
				"\n   | * Version         -- ${CMAKE_C_COMPILER_VERSION}"
				"\n   | * MPI Name        -- ${MPI_C_COMPILER}"
				"\n   | * MPI Exec.       -- ${MPIEXEC_EXECUTABLE}"
				"\n   ---------------------------------------------------------------")
message(STATUS "|                        INSTALLATION                         |"
				"\n   ---------------------------------------------------------------"
				"\n   | * CMake Exec.     -- ${CMAKE_COMMAND}"
				"\n   | * CMake Ver.      -- ${CMAKE_VERSION}"
				"\n   | * Prefix          -- ${CMAKE_INSTALL_PREFIX}"
				"\n   | * Downloads Dir   -- ${CMAKE_DOWNLOAD_DIRECTORY}"
				"\n   | * Build Type      -- ${CMAKE_BUILD_TYPE}"
				"\n   | * Compiler Flags  -- ${_comm3_flags_}"
				"\n   ---------------------------------------------------------------")
message(STATUS "|                         LIBRARIES                           |"
				"\n   ---------------------------------------------------------------")
change_list_output("${LIBM_LIBRARY}" _changed_list_)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
message(STATUS "| Linux Math      :"
			"\n   |------------------"
			"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
			"\n   |------------------")
endif()
change_list_output("${BLAS_LIBRARIES}" _changed_list_)
message(STATUS "| BLAS            :"
        "\n   |------------------"
				"\n   | * Found           -- ${BLAS_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
change_list_output("${LAPACK_LIBRARIES}" _changed_list_)
message(STATUS "| LAPACK          :"
        "\n   |------------------"
				"\n   | * Found           -- ${LAPACK_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
change_list_output("${ZLIB_LIBRARIES}" _changed_list_)
if(NOT COMPILE_ZLIB)
message(STATUS "| ZLIB            :"
        "\n   |------------------"
				"\n   | * Found           -- ${ZLIB_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
else()
	message(STATUS "| ZLIB            :"
        "\n   |------------------"
				"\n   | * Found           -- ${ZLIB_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${zlib_url}"
        "\n   |------------------")
endif()
change_list_output("${HDF5_Fortran_LIBRARIES}" _changed_list_)
if(NOT COMPILE_HDF5)
message(STATUS "| HDF5            :"
        "\n   |------------------"
				"\n   | * Found           -- ${HDF5_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
else()
message(STATUS "| HDF5            :"
        "\n   |------------------"
				"\n   | * Found           -- ${HDF5_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${hdf5_url}"
				"\n   | * Depends       :\n   |   - ZLIB\n   |   - LIBAEC"
        "\n   |------------------")
endif()
change_list_output("${FFTW_LIBRARIES}" _changed_list_)
if(NOT COMPILE_FFTW)
message(STATUS "| FFTW            :"
        "\n   |------------------"
				"\n   | * Found           -- ${FFTW_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
else()
message(STATUS "| FFTW            :"
        "\n   |------------------"
				"\n   | * Found           -- ${FFTW_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${fftw_url}"
        "\n   |------------------")
endif()
change_list_output("${CFITSIO_LIBRARIES}" _changed_list_)
if(NOT COMPILE_CFITSIO)
message(STATUS "| CFITSIO         :"
        "\n   |------------------"
				"\n   | * Found           -- ${CFITSIO_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
        "\n   |------------------")
elseif(NOT (COMPILE_CFITSIO AND CFITSIO_USE_CURL))
message(STATUS "| CFITSIO         :"
        "\n   |------------------"
				"\n   | * Found           -- ${CFITSIO_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${cfitsio_url}"
				"\n   | * Depends       :\n   |   - ZLIB"
        "\n   |------------------")
elseif((NOT COMPILE_CFITSIO) AND CFITSIO_USE_CURL)
message(STATUS "| CFITSIO         :"
        "\n   |------------------"
				"\n   | * Found           -- ${CFITSIO_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${cfitsio_url}"
				"\n   | * Depends       :\n   |   - ZLIB\n   |   - CURL"
        "\n   |------------------")
endif()
change_list_output("${HEALPIX_LIBRARIES}" _changed_list_)
if(NOT COMPILE_HEALPIX)
message(STATUS "| HEALPIX         :"
        "\n   |------------------"
				"\n   | * Found           -- ${HEALPIX_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}")
else()
message(STATUS "| HEALPIX         :"
        "\n   |------------------"
				"\n   | * Found           -- ${HEALPIX_FOUND}"
				"\n   | * Name(s)       :\n   |   - ${_changed_list_}"
				"\n   | * Source URL    :\n   |   - ${healpix_url}"
				"\n   | * Depends       :\n   |   - CFITSIO")
endif()
message(STATUS "===============================================================")
