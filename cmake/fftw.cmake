# Project: FFTW
# File which contains setup for current project 
# Author: Maksym Brilenkov

# looking for FFTW3 library on the system
find_package(FFTW3)
if(NOT FFTW3_FOUND)
	ExternalProject_Add(${project}
		URL "${${project}_url}"
		PREFIX "${CMAKE_DOWNLOAD_DIRECTORY}/${project}"
		DOWNLOAD_DIR "${CMAKE_DOWNLOAD_DIRECTORY}"
		BINARY_DIR "${CMAKE_DOWNLOAD_DIRECTORY}/${project}/src/${project}"
		INSTALL_DIR "${out_install_dir}"
		# commands how to build the project
		CONFIGURE_COMMAND "${${project}_configure_command}"
		COMMAND ./configure --prefix=<INSTALL_DIR> 
		#BUILD_IN_SOURCE 1	
		)
	set(FFTW3_LIBRARIES ${out_lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
	add_custom_target(${project} ALL "")
	message(STATUS "FFTW3 LFLAGS are: ${FFTW3_LFLAGS}")
	# this is an include dirs provided by custom FindFFTW3.cmake
	include_directories(${FFTW3_INCLUDES})
endif()
#message(STATUS
#"
#FFTW VERSION: ${FFTW3_VERSION}
#FFTW INCLUDE DIRS: ${FFTW3_INCLUDE_DIRS}
#FFTW LIBRARIES: ${FFTW3_LIBRARIES}
#")
#message("${${project}_configure_command}")




# adding the library to the project
#add_library(${project}_lib STATIC IMPORTED)
#set(${${project}_lib}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}3${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${${project}_lib} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${${${project}_lib}_name}")
#message("The ${${${project}_lib}_name} Path is " ${out_install_dir}/lib/${${${project}_lib}_name})
#include_directories(${out_install_dir}/include)

#add_library(${project} STATIC IMPORTED)
#set(lib_${project}_name ${CMAKE_STATIC_LIBRARY_PREFIX}${project}${CMAKE_STATIC_LIBRARY_SUFFIX})
#set_target_properties(${project} PROPERTIES IMPORTED_LOCATION "${out_install_dir}/lib/${lib_${project}_name}") 
