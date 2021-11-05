#================================================================================
# Author: Maksym Brilenkov
#================================================================================
# Description: This file contains Summary Information Output to the user's screen
# regarding System and CMake configuration which will be used to compile the code
#================================================================================
# TODO:
# Put these variables inside variables.cmake so you can decide which library to 
# compile how (with MMX, SSE support etc.)
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

#			if(CPU_IS_64BIT)
#				"\n   | * Is 64Bit        -- yes" 
#			else()
#				"\n   | * Is 64Bit        -- no" 
#			endif()

message(STATUS ${HOST_NAME})
message(STATUS "${HOST_OS_NAME} ${HOST_OS_PLATFORM} ${HOST_OS_RELEASE} ${HOST_OS_VERSION}")
message(STATUS "${PROC_NAME} | ${PROC_DESCRIPTION}")

message(STATUS "===============================================================")
message(STATUS "|                    CONFIGURATION SUMMARY                    |")
message(STATUS "===============================================================")
message(STATUS "|                          GENERAL                            |")
message(STATUS "---------------------------------------------------------------")
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
				"\n   | * Logical Cores   -- ${N_LOGICAL_CORES}"
        "\n   |------------------")
message(STATUS "| Memory          :"
        "\n   |------------------"
				"\n   | * Virtual [total] -- ${TOT_VIRTUAL_MEMORY} [MiB]"
				"\n   | * Virtual [avail] -- ${AVAIL_VIRTUAL_MEMORY} [MiB]"
				"\n   | * Phys    [total] -- ${TOT_PHYSICAL_MEMORY} [MiB]"
				"\n   | * Phys    [avail] -- ${AVAIL_PHYSICAL_MEMORY} [MiB]")
message(STATUS "---------------------------------------------------------------")
message(STATUS "|                         COMPILERS                           |")
message(STATUS "---------------------------------------------------------------")
message(STATUS "")
message(STATUS "| Name                   : ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "| OS Platform              : ${HOST_OS_PLATFORM}")
message(STATUS "| OS Release               : ${HOST_OS_RELEASE}")
message(STATUS "| OS Version                : ${HOST_OS_VERSION}")
message(STATUS "---------------------------------------------------------------")
message(STATUS "| CMake v.         : ${CMAKE_VERSION}")
message(STATUS "| CMake command         : ${CMAKE_COMMAND}")
message(STATUS "---------------------------------------------------------------")
message(STATUS "|                         COMPILERS                           |")
message(STATUS "---------------------------------------------------------------")
message(STATUS "| Fortran         : ${CMAKE_Fortran_COMPILER}")
message(STATUS "| Compiler ID     : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "| Compiler v.: ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "| C compiler              : ${CMAKE_C_COMPILER}")
message(STATUS "| C compiler id           : ${CMAKE_C_COMPILER_ID}")
message(STATUS "| C compiler version      : ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "| C++ compiler            : ${CMAKE_CXX_COMPILER}")
message(STATUS "| C++ compiler id         : ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "| C++ compiler version    : ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Using ccache if found : ${USE_CCACHE}")
message(STATUS "  BLAS                : ${BLAS_INFO}")
message(STATUS "  BLAS_HAS_SBGEMM     : ${BLAS_HAS_SBGEMM}")
message(STATUS "  LAPACK              : ${LAPACK_INFO}")
message(STATUS "Hardware")
message(STATUS "Compilers:")
message(STATUS "")
message(STATUS "---------------------------------------------------------------")
message(STATUS "")
message(STATUS "")
message(STATUS "---------------------------------------------------------------")
message(STATUS "")
message(STATUS "")
message(STATUS "===============================================================")
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
