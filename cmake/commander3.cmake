# Project: Commander3 
# Link: https://gitlab.mpcdf.mpg.de/mtr/libsharp/-/tree/master
# This file contains instructions on how to install commander3
# on your machine.
# Author: Maksym Brilenkov

message(STATUS "---------------------------------------------------------------")

# To resolve conflict with the flags, we need to compile *.cpp
# into a stand alone library and then just link it to rest of 
# the files. However, with CMake it is possible to compile the
# file into object, which will not create *.a files. We go with
# this approach.
add_library(comm_system_backend 
	STATIC 
	${COMMANDER3_SOURCE_DIR}/comm_system_backend.cpp
	)
target_compile_options(comm_system_backend
	PUBLIC
	# setting flags depending on configuration
	#"$<$<CONFIG:Release>:${COMMANDER3_CXX_COMPILER_FLAGS_RELEASE}>"
	#"$<$<CONFIG:Debug>:${COMMANDER3_CXX_COMPILER_FLAGS_DEBUG}>"
	#"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_CXX_COMPILER_FLAGS_RELWITHDEBINFO}>"
	#"$<$<CONFIG:MinSizeRel>:${COMMANDER3_CXX_COMPILER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_CXX_COMPILER_FLAGS}
	)
#target_link_options(${comm_system_backend}
#	PUBLIC
#	"$<$<CONFIG:Release>:${COMMANDER3_CXX_LINKER_FLAGS_RELEASE}>"
#	"$<$<CONFIG:Debug>:${COMMANDER3_CXX_LINKER_FLAGS_DEBUG}>"
#	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_CXX_LINKER_FLAGS_RELWITHDEBINFO}>"
#	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_CXX_LINKER_FLAGS_MINSIZEREL}>"
#	# setting other compiler dependent flags
#	${COMMANDER3_CXX_LINKER_FLAGS}
#	)

# setting the directory where to output all .mod and .o files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/mod)
# TODO: add all sources manually instead of this command, as 
# there seems to be a problem with tempita language
#file(GLOB_RECURSE sources *.f90 *.cpp *.f)
set(sources
	${COMMANDER3_SOURCE_DIR}/commander.f90
	${COMMANDER3_SOURCE_DIR}/ARS_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_fft_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_QUcov_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_4D_map_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_0D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_rms_mod.f90
	# TOD processing/simulations modules
	${COMMANDER3_SOURCE_DIR}/comm_tod_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_mapmaking_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_LFI_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_gain_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_noise_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_orbdipole_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_pointing_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_WMAP_mod.f90
	#${COMMANDER3_SOURCE_DIR}/comm_tod_simulations_mod.f90
	#
	${COMMANDER3_SOURCE_DIR}/comm_F_int_1D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_output_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_B_bl_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_2D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_param_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_beam_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_int_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_physdust_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_B_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_line_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_powlaw_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_bp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_F_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_ptsrc_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_zodi_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_bp_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_freefree_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_shared_arr_mod.f90
	${COMMANDER3_SOURCE_DIR}/d1mach.f
	${COMMANDER3_SOURCE_DIR}/comm_chisq_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_gain_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_shared_output_mod.f90
	${COMMANDER3_SOURCE_DIR}/drc3jj.f
	${COMMANDER3_SOURCE_DIR}/comm_Cl_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_hdf_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_signal_mod.f90
	${COMMANDER3_SOURCE_DIR}/hashtbl_4Dmap.f90
	${COMMANDER3_SOURCE_DIR}/comm_cmb_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_spindust2_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/hashtbl.f90
	${COMMANDER3_SOURCE_DIR}/comm_cmb_relquad_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_huffman_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_spindust_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/InvSamp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_line_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_status_mod.f90
	${COMMANDER3_SOURCE_DIR}/locate_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_conviqt_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_map_mod.f90
	#comm_system_backend.cpp
	${COMMANDER3_SOURCE_DIR}/math_tools.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_MBB_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_system_mod.f90
	${COMMANDER3_SOURCE_DIR}/powell_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_precond_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_md_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_task_mod.f90
	${COMMANDER3_SOURCE_DIR}/sharp.f90
	${COMMANDER3_SOURCE_DIR}/comm_cr_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_mpi_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_template_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/sort_utils.f90
	${COMMANDER3_SOURCE_DIR}/comm_data_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_N_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_tod_bandpass_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_defs.f90
	${COMMANDER3_SOURCE_DIR}/comm_noise_mod.f90
	${COMMANDER3_SOURCE_DIR}/spline_1D_mod.f90
	${COMMANDER3_SOURCE_DIR}/spline_2D_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_diffuse_comp_mod.f90
	${COMMANDER3_SOURCE_DIR}/comm_nonlin_mod.f90
	#progressbar_mod.f90
	)

# manually setting Fortran compiler flags for different compilers
#if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	# To prevent error with length and turn type mismatch into a warnings
	# -fallow-argument-mismatch <= is for 10x compilers
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -ffree-line-length-none -fallow-argument-mismatch -std=legacy")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl -ffree-line-length-none -fallow-argument-mismatch")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fallow-argument-mismatch")
	# -Wno-argument-mismatch <= is th flag for 9x and older compilers
	#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-argument-mismatch")
	#	message(STATUS "CMAKE_Fortran_FLAGS are: ${CMAKE_Fortran_FLAGS}")
	#elseif("${CMAKE_Fortran_COMPILER}" MATCHES "Intel")
	#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -qopt-matmul")
	#endif()
# setting executable name
set(commander3 commander3)
add_executable(${commander3} "")
# make sure that commander executable will be built last
add_dependencies(${commander3} ${projects}) #fftw_ fftw_f)
target_sources(${commander3}
	PUBLIC	
	${sources}
	)
set_property(TARGET ${commander3} PROPERTY ENABLE_EXPORTS TRUE)
set_property(TARGET ${commander3} PROPERTY LINKER_LANGUAGE Fortran)
# adding compiler flags to commander3 target
target_compile_options(${commander3}
	PUBLIC
	# setting flags depending on configuration
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_COMPILER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_COMPILER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_COMPILER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_COMPILER_FLAGS}
	)
# adding linker flags to commander3 target
target_link_options(${commander3}
	PUBLIC
	"$<$<CONFIG:Release>:${COMMANDER3_Fortran_LINKER_FLAGS_RELEASE}>"
	"$<$<CONFIG:Debug>:${COMMANDER3_Fortran_LINKER_FLAGS_DEBUG}>"
	"$<$<CONFIG:RelWithDebInfo>:${COMMANDER3_Fortran_LINKER_FLAGS_RELWITHDEBINFO}>"
	"$<$<CONFIG:MinSizeRel>:${COMMANDER3_Fortran_LINKER_FLAGS_MINSIZEREL}>"
	# setting other compiler dependent flags
	${COMMANDER3_Fortran_LINKER_FLAGS}
	)

#message("cmake dl libs are ${CMAKE_DL_LIBS}")
# LINKING ORDER IN LIBRARIES IS IMPORTANT!
target_link_libraries(${commander3} 
	PUBLIC	
	# linking MPI
	MPI::MPI_Fortran
	# linking OpenMP
	OpenMP::OpenMP_Fortran
	# including MKL
	#-qopt-matmul
	${BLAS_LINKER_FLAGS} 
	${BLAS_LIBRARIES}
	${LAPACK_LINKER_FLAGS} 
	${LAPACK_LIBRARIES}
	#-ffree-line-length-none
	#-fno-strict-overflow
	#hdf5_lib
	#${HDF5_Fortran_LIBRARIES}
	#${HDF5_Fortran_LIBRARIES}
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhdf5_fortran.a" #<= getting other errors when linking this one as well
	# including sharp2
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libsharp2.a"
	#"${out_lib_dir}/libsharp2.a"
	#${SHARP2_LIBRARIES}
	# including healpix
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhealpix.a"
	#"${out_lib_dir}/libhealpix.a"
	${HEALPIX_LIBRARIES}
	# including cfitsio
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libcfitsio.a"
	#"${out_lib_dir}/libcfitsio.a"
	${CFITSIO_LIBRARIES}
	# including hdf5 - first fortran and then general
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhdf5_fortran.a" #<= getting other errors when linking this one as well
	#"${out_lib_dir}/libhdf5_fortran.a"
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhdf5.a"
	#"${out_lib_dir}/libhdf5.a"
	# to avoid error error dlclose@@GLIBC_2.2.5', so 
	# we need to link math library
	-lm
	# and -ldl (dl library)
	${CMAKE_DL_LIBS}
	${HDF5_Fortran_LIBRARIES}
	#${HDF5_LIBRARIES}
	#${hdf5_fortran}
	# hdf5 requires zlib (?), otherwise will get some stupid error
	#"-lz"
	#"/usr/lib64/libz.so"
	#ZLIB::ZLIB
	#zlib_lib
	#${ZLIB_LIBRARIES}
	# adding curl
	#${CURL_LIBRARIES}
	# these are sort of curl dependencies
	#-lcrypto 
	#-lssl
	#CURL::libcurl
	# now, I am not even sure we need to include this one :)
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libfftw3.a"
	#"${out_lib_dir}/libfftw3.a"
	${FFTW3_LIBRARIES}
	# Linking commander *.cpp file(s)
	comm_system_backend
	)

# installing commander into appropriate folder
install(TARGETS ${commander3} RUNTIME DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
