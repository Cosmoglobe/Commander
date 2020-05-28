# Project: OpenMP 
# File determines location of OpenMP on the system  
# Author: Maksym Brilenkov

find_package(OpenMP REQUIRED)

# printing out status of the search
#message(STATUS "
#MPI Libs Path: ${MPI_Fortran_LIBRARIES}
#MPI Include Path: ${MPI_Fortran_INCLUDE_PATH}
#MPI Compile Options are: ${MPI_Fortran_COMPILE_OPTIONS} 
#MPI Link options are: ${MPI_Fortran_LINK_FLAGS}
#THREADS Libs Path: ${CMAKE_THREAD_LIBS_INIT}"#MPIexec: ${MPIEXEC_EXECUTABLE}"
#)

# to avoid cmake errors we create and empty target
add_custom_target(${project} ALL "")
set(${project}_lib OpenMP::OpenMP_Fortran)
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})

