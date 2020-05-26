# This file contains general instructions how to
# fetch and build the commander dependencies
# Author: Maksym Brilenkov

#=================================================
# Setting up default values for
# the main project variables
set(CMAKE_DOWNLOAD_DIRECTORY "${CMAKE_SOURCE_DIR}/build/downloads"
	CACHE STRING
	"Directory where to download commander dependencies' source files")
set(CMAKE_INSTALL_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/install"
	CACHE STRING
	"Directory where to install commander dependencies")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_OUTPUT_DIRECTORY}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_OUTPUT_DIRECTORY}/bin")
set(DOXYGEN_BUILD_DOCS TRUE 
	CACHE BOOL
	"Determine whether to use doxygen or not")

# Commander source dir
set(COMMANDER3_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/commander")

# tempita source dir
set(tempita_dir ${CMAKE_SOURCE_DIR}/src/python)
# setting dir to download .tar.gz files 
#set(download_dir ${CMAKE_SOURCE_DIR}/build/downloads)
# setting an install dir
#set(out_install_dir ${CMAKE_SOURCE_DIR}/build/install)
# setting dir for output binaries
#set(out_bin_dir ${CMAKE_SOURCE_DIR}/build/install/bin)
# setting dir for output libraries
#set(out_lib_dir ${CMAKE_SOURCE_DIR}/build/install/lib)

# adding custom FindSomething.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


unset(projects)

# project names
list(APPEND projects 
	tempita
	blas # blas-lapack module 
	mpi
	openmp
	curl
	zlib
	sharp2
	fftw
	cfitsio
	hdf5
	doxygen
	healpix
	)

# projects urls
set(curl_url "https://github.com/curl/curl/releases/download/curl-7_69_0/curl-7.69.0.zip")#"https://github.com/curl/curl/releases/download/curl-7_69_1/curl-7.69.1.tar.gz")
set(fftw_url "http://fftw.org/fftw-3.3.8.tar.gz")
set(hdf5_url "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz")#"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz")
set(sharp2_url "https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz")#"https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/archive/master/libsharp-master.tar.gz") #"https://github.com/Libsharp/libsharp/archive/v1.0.0.tar.gz")
set(cfitsio_url "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz")
set(healpix_url "https://sourceforge.net/projects/healpix/files/Healpix_3.50/Healpix_3.50_2018Dec10.zip/download")#"https://sourceforge.net/projects/healpix/files/Healpix_3.60/Healpix_3.60_2019Dec18.zip/download")#"https://sourceforge.net/projects/healpix/files/latest/download")
# there is some weird errors appearing for doxygen v1.8.17 and above, so will stick with this one
set(doxygen_url "https://github.com/doxygen/doxygen/archive/Release_1_8_16.tar.gz")#"https://github.com/doxygen/doxygen/archive/Release_1_8_18.tar.gz")#"https://github.com/doxygen/doxygen.git")
set(flex_url "http://sourceforge.net/projects/flex/files/flex-2.5.39.tar.gz/download")#"https://sourceforge.net/projects/flex/files/flex-2.6.0.tar.xz/download")
set(bison_url "http://ftp.gnu.org/gnu/bison/bison-3.6.tar.gz")#"http://ftp.gnu.org/gnu/bison/bison-3.6.2.tar.gz")

# projects configure commands (by default it is assumes to be cmake)
set(curl_configure_command "")
set(fftw_configure_command "")
set(hdf5_configure_command "")
set(sharp2_configure_command "autoreconf" "-i")# <= you need to specify command separately, othewise it won't work#"autoconf")
set(cfitsio_configure_command "")
set(healpix_configure_command "")
set(tempita_configure_command "")

# projects install commands
#set(fftw_install_command "make install")
##set(hdf5_install_command )
##set(sharp2_install_command "")
#set(cfitsio_install_command )

#message(${sharp_install_command})



# include all project configuration files
foreach(project ${projects})
	include("${project}")
	#target_link_libraries(${${project}_lib} ${MPI_Fortran_LIBRARIES})
endforeach()
# 
include_directories(${CMAKE_INSTALL_OUTPUT_DIRECTORY}/include)
#include_directories(${out_install_dir}/lib)
#include_directories(${out_install_dir}/mod)
# this one here is not advised to use, but desperate times need desperate measures
#link_directories(${out_install_dir}/mod)
