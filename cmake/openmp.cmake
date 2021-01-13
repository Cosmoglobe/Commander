# Project: OpenMP 
# File determines location of OpenMP on the system  
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
find_package(OpenMP REQUIRED)

# to avoid cmake errors we create and empty target
add_custom_target(${project} ALL "")
set(${project}_lib OpenMP::OpenMP_Fortran)
# setting compilation and linking flags
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${OpenMP_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${OpenMP_Fortran_LIBRARIES})

message(STATUS "OPENMP Fortran LIBRARIES are: ${OpenMP_Fortran_LIBRARIES}")
