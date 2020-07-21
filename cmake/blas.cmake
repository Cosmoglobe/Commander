# Project: BLAS 
# File determines location of BLAS on the system  
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")
# require BLAS and LAPACK
# From docs: Note C, CXX or Fortran must be enabled
# to detect a BLAS/LAPACK library. C or CXX must be
# enabled to use Intel Math Kernel Library (MKL).
#set(BLA_VENDOR Intel10_64lp)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# to avoid cmake errors we create and empty target
add_custom_target(${project} ALL "")
set(${project}_lib ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
message(STATUS "BLAS LINKER FLAGS:   ${BLAS_LINKER_FLAGS}")
message(STATUS "BLAS LIBRARIES:      ${BLAS_LIBRARIES}")
message(STATUS "LAPACK LINKER FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK LIBRARIES:    ${LAPACK_LIBRARIES}")
