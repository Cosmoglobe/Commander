# Project: MPI 
# File determines location of MPI on the system  
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
find_package(Threads)
# printing out status of the search
message(STATUS "MPI Libs Path: ${MPI_Fortran_LIBRARIES}")
message(STATUS "MPI Include Path: ${MPI_Fortran_INCLUDE_PATH}")
message(STATUS "MPI Compile Options are: ${MPI_Fortran_COMPILE_OPTIONS}") 
message(STATUS "MPI Link options are: ${MPI_Fortran_LINK_FLAGS}")
message(STATUS "THREADS Libs Path: ${CMAKE_THREAD_LIBS_INIT}")#MPIexec: ${MPIEXEC_EXECUTABLE}"

# to avoid cmake errors we create and empty target
add_custom_target(${project} ALL "")
set(${project}_lib MPI::MPI_Fortran)
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${MPI_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES} Threads::Threads)

