# Project: LibSharp
# File which contains setup for the project 
# Author: Maksym Brilenkov

# define this variable here for easier reference in the future

#message("MY CMAKE FLAGS ARE "${CMAKE_C_FLAGS})#"${${project}_configure_command}")

ExternalProject_Add(${project}
	URL "${${project}_url}"
	PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
	DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
	BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
	INSTALL_DIR "${CMAKE_INSTALL_OUTPUT_DIRECTORY}"
	# commands how to build the project
	CONFIGURE_COMMAND "${${project}_configure_command}"
	#COMMAND CFLAGS="-DUSE_MPI" ${download_dir}/${project}/src/${project}/configure --prefix=<INSTALL_DIR>
	COMMAND CFLAGS=-DUSE_MPI ./configure #--prefix=<INSTALL_DIR>
	#BUILD_IN_SOURCE 1	
	INSTALL_COMMAND ""
	#DEPENDS mpi
	# LibSharp doesn't have an install command, so we need to manually move files
	# into installation folder and hope they will work :)
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}/.libs" "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" #"${out_install_dir}/lib" 
	)

# defining the variable which will show the path to the compiled libraries
set(SHARP2_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
# creating a library statis library:
# addign standard cmake suffixes 
# (cmake will figure stuff out depending on the platform)
#add_library(sharp2_lib STATIC ${out_lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(sharp2_lib PROPERTIES LINKER_LANGUAGE Fortran)
#message(STATUS "SHARP2 LIBS ARE: ${sharp2_lib}")

#target_compile_options(${project} PUBLIC -DUSE_MPI)

# linking libraries
#add_library(${project}_lib STATIC IMPORTED)
#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")

#message("${${project}_lib}")
#add_library(${project}_lib STATIC ${out_install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})

#target_link_libraries(${${project}_lib} MPI::MPI_Fortran)
#message("THIS IS LIBSHARP TEST MUHAHAHAHAHAH ${project}_lib")
#message(${out_install_dir}/lib/${${${project}_lib}_name})
#########file(COPY "${download_dir}/${project}/src/${project}/auto"  DESTINATION "${out_install_dir}")

# linking libraries to LAPACK etc.
# Linking BLAS library
#target_link_libraries(${${project}_lib} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
# Linking LAPACK library
#target_link_libraries(${${project}_lib} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
# Linking MPI
#target_link_libraries(${${project}_lib} ${MPI_Fortran_LINK_FLAGS} ${MPI_fortran_LINKER_FLAGS} ${MPI_Fortran_COMPILE_FLAGS} ${MPI_Fortran_LIBRARIES})
