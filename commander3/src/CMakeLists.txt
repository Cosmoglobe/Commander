# setting the directory where to output all .mod and .o files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_OUTPUT_DIRECTORY}/mod)
# TODO: add all sources manually instead of this command, as 
# there seems to be a problem with tempita language
#file(GLOB_RECURSE sources *.f90 *.cpp *.f)
set(sources
	commander.f90
	ARS_mod.f90
	comm_fft_mod.f90
	comm_N_QUcov_mod.f90
	comm_tod_mapmaking_mod.f90
	comm_4D_map_mod.f90
	comm_F_int_0D_mod.f90
	comm_N_rms_mod.f90
	comm_tod_mod.f90
	comm_tod_LFI_mod.f90
	comm_tod_gain_mod.f90
	comm_F_int_1D_mod.f90
	comm_output_mod.f90
	comm_tod_noise_mod.f90
	comm_B_bl_mod.f90
	comm_F_int_2D_mod.f90
	comm_param_mod.f90
	comm_tod_orbdipole_mod.f90
	comm_beam_mod.f90
	comm_F_int_mod.f90
	comm_physdust_comp_mod.f90
	comm_tod_pointing_mod.f90
	comm_B_mod.f90
	comm_F_line_mod.f90
	comm_powlaw_comp_mod.f90
	comm_utils.f90
	comm_bp_mod.f90
	comm_F_mod.f90
	comm_ptsrc_comp_mod.f90
	comm_zodi_mod.f90
	comm_bp_utils.f90
	comm_freefree_comp_mod.f90
	comm_shared_arr_mod.f90
	d1mach.f
	comm_chisq_mod.f90
	comm_gain_mod.f90
	comm_shared_output_mod.f90
	drc3jj.f
	comm_Cl_mod.f90
	comm_hdf_mod.f90
	comm_signal_mod.f90
	hashtbl_4Dmap.f90
	comm_cmb_comp_mod.f90
	#comm_hdf_mod.f90.in
	comm_spindust2_comp_mod.f90
	hashtbl.f90
	comm_cmb_relquad_comp_mod.f90
	comm_huffman_mod.f90
	comm_spindust_comp_mod.f90
	InvSamp_mod.f90
	comm_comp_mod.f90
	comm_line_comp_mod.f90
	comm_status_mod.f90
	locate_mod.f90
	comm_conviqt_mod.f90
	comm_map_mod.f90
	comm_system_backend.cpp
	math_tools.f90
	comm_cr_mod.f90
	comm_MBB_comp_mod.f90
	comm_system_mod.f90
	powell_mod.f90
	comm_cr_precond_mod.f90
	comm_md_comp_mod.f90
	comm_task_mod.f90
	sharp.f90
	comm_cr_utils.f90
	comm_mpi_mod.f90
	comm_template_comp_mod.f90
	sort_utils.f90
	comm_data_mod.f90
	comm_N_mod.f90
	comm_tod_bandpass_mod.f90
	comm_defs.f90
	comm_noise_mod.f90
	spline_1D_mod.f90
	spline_2D_mod.f90
	comm_diffuse_comp_mod.f90
	comm_nonlin_mod.f90
	comm_tod_WMAP_mod.f90
	)

# manually setting Fortran compiler flags for different compilers
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	# To prevent error with length and turn type mismatch into a warnings
	# -fallow-argument-mismatch <= is for 10x compilers
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -ffree-line-length-none -fallow-argument-mismatch -std=legacy")
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl -ffree-line-length-none -fallow-argument-mismatch")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fallow-argument-mismatch")
	# -Wno-argument-mismatch <= is th flag for 9x and older compilers
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-argument-mismatch")
	message(STATUS "CMAKE_Fortran_FLAGS are: ${CMAKE_Fortran_FLAGS}")
endif()
#
set(executable commander3)
add_executable(${executable} "")
# make sure that commander executable will be built last
add_dependencies(${executable} ${projects})
target_sources(${executable}
	PUBLIC	
	${sources}
	)
set_property(TARGET ${executable} PROPERTY ENABLE_EXPORTS TRUE)
set_property(TARGET ${executable} PROPERTY LINKER_LANGUAGE Fortran)

message("cmake dl libs are ${CMAKE_DL_LIBS}")
# the linking order in libraries is important
target_link_libraries(${executable} 
	PUBLIC	
	# linking MPI
	MPI::MPI_Fortran
	# linking OpenMP
	OpenMP::OpenMP_Fortran
	#-ffree-line-length-none
	#-fno-strict-overflow
	#hdf5_lib
	#${HDF5_Fortran_LIBRARIES}
	#${HDF5_Fortran_LIBRARIES}
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhdf5_fortran.a" #<= getting other errors when linking this one as well
	# including sharp2
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libsharp2.a"
	#"${out_lib_dir}/libsharp2.a"
	${SHARP2_LIBRARIES}
	# including healpix
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libhealpix.a"
	#"${out_lib_dir}/libhealpix.a"
	${HEALPIX_LIBRARIES}
	# including cfitsio
	#"/mn/stornext/u3/maksymb/cmake_tests/CommanderSuperbuild/build/install/lib/libcfitsio.a"
	#"${out_lib_dir}/libcfitsio.a"
	${CFITSIO_LIBRARIES}
	# including MKL
	${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
	${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
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
	)

# installing commander into appropriate folder
install(TARGETS ${executable} RUNTIME DESTINATION ${CMAKE_INSTALL_OUTPUT_DIRECTORY}/bin)
