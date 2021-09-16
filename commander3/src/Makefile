#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
# 
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

.SUFFIXES:

FSOURCES=d1mach.f drc3jj.f
F90SOURCES=hashtbl.f90 \
           hashtbl_4dmap.f90 \
           powell_mod.f90 \
           comm_hdf_mod.f90 \
           sharp.f90 \
           sort_utils.f90 \
           ars_mod.f90 \
           math_tools.f90 \
           locate_mod.f90 \
           spline_1d_mod.f90 \
           spline_2d_mod.f90 \
           invsamp_mod.f90 \
           comm_mpi_mod.f90 \
           comm_system_backend.f90 \
           comm_system_mod.f90 \
           comm_map_mod.f90 \
           comm_4d_map_mod.f90 \
           comm_defs.f90 \
           comm_conviqt_mod.f90 \
	   comm_utils.f90 \
           comm_shared_arr_mod.f90 \
           comm_huffman_mod.f90 \
           comm_cr_utils.f90 \
           comm_cr_precond_mod.f90 \
           comm_zodi_mod.f90 \
           comm_shared_output_mod.f90 \
           comm_status_mod.f90 \
	   comm_bp_utils.f90 \
           comm_param_mod.f90 \
           comm_comp_mod.f90 \
           comm_diffuse_comp_mod.f90 comm_diffuse_comp_smod.f90 \
	   comm_cmb_comp_mod.f90 \
           comm_cmb_relquad_comp_mod.f90 \
           comm_powlaw_comp_mod.f90 \
           comm_physdust_comp_mod.f90 \
           comm_mbb_comp_mod.f90 \
           comm_freefree_comp_mod.f90 \
           comm_line_comp_mod.f90 \
           comm_spindust_comp_mod.f90 \
           comm_spindust2_comp_mod.f90 \
           comm_ame_lognormal_mod.f90 \
           comm_md_comp_mod.f90 \
           comm_template_comp_mod.f90 \
           comm_ptsrc_comp_mod.f90 \
	   comm_gain_mod.f90 \
           comm_nonlin_mod.f90 \
           comm_n_mod.f90 \
           comm_n_rms_mod.f90 \
           comm_n_qucov_mod.f90 \
           comm_b_mod.f90 \
           comm_b_bl_mod.f90 \
           comm_beam_mod.f90 \
           comm_f_int_mod.f90 \
           comm_f_int_1d_mod.f90 \
           comm_f_int_0d_mod.f90 \
           comm_f_int_2d_mod.f90 \
           comm_fft_mod.f90 \
           comm_tod_noise_psd_mod.f90 \
           comm_tod_mod.f90 \
           comm_tod_mapmaking_mod.f90 \
           comm_tod_bandpass_mod.f90 \
           comm_tod_gain_mod.f90 \
           comm_tod_pointing_mod.f90 \
           comm_tod_lfi_mod.f90 \
           comm_tod_jump_mod.f90 \
           comm_tod_spider_mod.f90 \
           comm_tod_orbdipole_mod.f90 \
           comm_tod_noise_mod.f90 \
           comm_tod_driver_mod.f90 \
           comm_f_line_mod.f90 \
           comm_data_mod.f90 \
           comm_bp_mod.f90 \
           comm_cl_mod.f90 \
           comm_chisq_mod.f90 \
           comm_cr_mod.f90 \
           comm_signal_mod.f90 \
           comm_output_mod.f90 \
           comm_tod_wmap_mod.f90 \
	   comm_tod_simulations_mod.f90 \
           comm_tod_lb_mod.f90 \
           comm_camb_mod.f90 \
           comm_tod_adc_mod.f90 \
           commander.f90

commander: $(subst .f,.o,$(FSOURCES)) $(subst .f90,.o,$(F90SOURCES))
	$(MPF90) -o commander $+ $(LINK) $(MPFCLIBS)

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod commander

%.o %.mod %.smod: %.f90 
	$(MPF90) $(F90COMP) -c $<
	@touch $@

%.o %.mod %.smod: %.f
	$(MPF90) $(F90COMP) -c $<
	@touch $@

%.o %.mod %.smod: %.cpp
	$(MPCXX) $(CXXCOMP) -c $<
	@touch $@


locate_mod.o :             math_tools.mod
spline_1d_mod.o :          locate_mod.mod math_tools.mod
spline_2d_mod.o :          comm_utils.mod
invsamp_mod.o :            spline_1d_mod.mod
comm_map_mod.o :           sharp.mod comm_hdf_mod.mod comm_param_mod.mod
comm_conviqt_mod.o :       sharp.mod comm_map_mod.mod comm_shared_arr_mod.mod
comm_hdf_mod.o :           comm_utils.mod
comm_fft_mod.o :           locate_mod.mod comm_utils.mod comm_param_mod.mod
comm_tod_noise_psd_mod.o : comm_utils.mod
comm_tod_mod.o :           comm_tod_noise_psd_mod.mod comm_fft_mod.mod comm_map_mod.mod comm_hdf_mod.mod comm_param_mod.mod \
			   comm_huffman_mod.mod comm_zodi_mod.mod comm_tod_orbdipole_mod.mod
comm_tod_mapmaking_mod.o:  comm_tod_mod.mod
comm_tod_bandpass_mod.o:   comm_tod_mod.mod
comm_tod_noise_mod.o:      comm_tod_mod.mod comm_tod_jump_mod.mod invsamp_mod.mod
comm_tod_gain_mod.o:       comm_tod_mod.mod comm_tod_noise_mod.mod
comm_tod_orbdipole_mod.o:  comm_map_mod.mod comm_utils.mod
comm_tod_pointing_mod.o:   comm_tod_mod.mod
comm_utils.o :             comm_defs.mod spline_1d_mod.mod
comm_shared_arr_mod.o :    comm_utils.mod
comm_huffman_mod.o :       comm_utils.mod
comm_system_mod.o :        comm_system_backend.o
comm_shared_output_mod.o : comm_system_mod.mod
comm_status_mod.o :        comm_shared_output_mod.mod
comm_utils.o :             sort_utils.mod comm_system_mod.mod sharp.mod comm_mpi_mod.mod
ars_mod.o :                sort_utils.mod
comm_bp_utils.o :          comm_utils.mod
comm_param_mod.o :         comm_utils.mod comm_status_mod.mod hashtbl.mod
comm_cr_utils.o :          comm_utils.mod
comm_cr_precond_mod.o :    comm_utils.mod
comm_comp_mod.o :          comm_bp_mod.mod comm_param_mod.mod comm_cr_utils.mod comm_data_mod.mod
comm_cl_mod.o :            comm_map_mod.mod
comm_camb_mod.o :          comm_utils.mod 
comm_diffuse_comp_mod.o :  comm_comp_mod.mod comm_data_mod.mod comm_f_int_mod.mod comm_cl_mod.mod \
                           comm_cr_precond_mod.mod comm_beam_mod.mod
comm_diffuse_comp_smod.o : comm_diffuse_comp_mod.mod 
comm_diffuse_comp_smod.o : comm_diffuse_comp_mod.smod
comm_template_comp_mod.o : comm_comp_mod.mod
comm_ptsrc_comp_mod.o :    comm_template_comp_mod.mod comm_f_int_0d_mod.mod comm_f_int_2d_mod.mod \
                           comm_hdf_mod.mod comm_cr_precond_mod.mod powell_mod.mod
comm_cmb_comp_mod.o :      comm_diffuse_comp_mod.mod comm_f_int_0d_mod.mod
comm_cmb_relquad_comp_mod.o : comm_template_comp_mod.mod comm_f_int_0d_mod.mod
comm_powlaw_comp_mod.o :   comm_diffuse_comp_mod.mod comm_f_int_1d_mod.mod
comm_physdust_comp_mod.o : comm_diffuse_comp_mod.mod comm_f_int_1d_mod.mod
comm_spindust_comp_mod.o : comm_diffuse_comp_mod.mod comm_f_int_1d_mod.mod
comm_ame_lognormal_mod.o : comm_diffuse_comp_mod.mod comm_f_int_1d_mod.mod
comm_spindust2_comp_mod.o : comm_diffuse_comp_mod.mod comm_f_int_1d_mod.mod
comm_mbb_comp_mod.o :      comm_diffuse_comp_mod.mod comm_f_int_2d_mod.mod
comm_freefree_comp_mod.o : comm_diffuse_comp_mod.mod comm_f_int_2d_mod.mod
comm_line_comp_mod.o :     comm_diffuse_comp_mod.mod comm_f_line_mod.mod
comm_md_comp_mod.o :       comm_diffuse_comp_mod.mod comm_f_line_mod.mod
comm_signal_mod.o :        comm_powlaw_comp_mod.mod comm_cmb_comp_mod.mod comm_cmb_relquad_comp_mod.mod \
                           comm_spindust_comp_mod.mod comm_spindust2_comp_mod.mod comm_mbb_comp_mod.mod \
                           comm_freefree_comp_mod.mod comm_line_comp_mod.mod \
                           comm_md_comp_mod.mod comm_template_comp_mod.mod \
                           comm_cr_mod.mod comm_physdust_comp_mod.mod comm_ame_lognormal_mod.o
comm_nonlin_mod.o :        comm_data_mod.mod comm_comp_mod.mod comm_chisq_mod.mod \
                           comm_gain_mod.mod comm_output_mod.mod comm_signal_mod.mod \
                           comm_utils.mod
comm_gain_mod.o :          comm_data_mod.mod comm_comp_mod.mod comm_chisq_mod.mod
comm_n_mod.o :             comm_param_mod.mod
comm_n_rms_mod.o :         comm_n_mod.mod
comm_noise_mod.o :         comm_n_rms_mod.mod comm_n_qucov_mod.mod
comm_bp_mod.o :            comm_bp_utils.mod comm_param_mod.mod
comm_f_int_mod.o :         comm_utils.mod
comm_f_line_mod.o :        comm_f_int_mod.mod 
comm_f_int_0D_mod.o :      comm_f_int_mod.mod 
comm_f_int_1D_mod.o :      comm_f_int_mod.mod spline_1d_mod.mod
comm_f_int_2D_mod.o :      comm_f_int_mod.mod spline_2d_mod.mod comm_utils.mod
comm_b_mod.o :             comm_map_mod.mod
comm_b_bl_mod.o :          comm_b_mod.mod 
comm_beam_mod.o :          comm_b_bl_mod.mod
comm_data_mod.o :          comm_map_mod.mod comm_noise_mod.mod comm_bp_mod.mod comm_beam_mod.mod comm_tod_mod.mod \
                           comm_tod_lfi_mod.mod comm_tod_spider_mod.mod comm_tod_wmap_mod.mod comm_tod_lb_mod.mod
comm_cr_mod.o :            comm_data_mod.mod comm_comp_mod.mod \
                           comm_template_comp_mod.mod comm_diffuse_comp_mod.mod comm_output_mod.mod
comm_chisq_mod.o :         comm_data_mod.mod comm_comp_mod.mod
comm_output_mod.o :        comm_comp_mod.mod comm_hdf_mod.mod comm_cl_mod.mod
comm_4d_map_mod.o :        hashtbl_4dmap.mod comm_hdf_mod.mod comm_utils.mod
comm_zodi_mod.o :          comm_utils.mod comm_param_mod.mod comm_bp_utils.mod
comm_tod_driver_mod.o :    comm_tod_mod.mod comm_zodi_mod.mod \
                           comm_tod_mapmaking_mod.mod comm_tod_bandpass_mod.mod \
                           comm_tod_pointing_mod.mod comm_tod_orbdipole_mod.mod \
                           comm_tod_gain_mod.mod comm_tod_noise_mod.mod \
                           comm_tod_simulations_mod.mod comm_shared_arr_mod.mod 
comm_tod_lfi_mod.o :       comm_tod_driver_mod.mod comm_shared_arr_mod.mod \
                           comm_conviqt_mod.mod comm_4d_map_mod.mod comm_tod_adc_mod.mod 
comm_tod_spider_mod.o :    comm_tod_mod.mod comm_shared_arr_mod.mod comm_conviqt_mod.mod comm_4d_map_mod.mod comm_zodi_mod.mod \
                           comm_tod_mapmaking_mod.mod comm_tod_bandpass_mod.mod comm_tod_pointing_mod.mod comm_tod_orbdipole_mod.mod comm_tod_gain_mod.mod comm_tod_noise_mod.mod \
                           comm_tod_jump_mod.mod 
comm_tod_jump_mod.o :      comm_utils.mod sort_utils.mod
comm_tod_wmap_mod.o :      comm_tod_mod.mod comm_shared_arr_mod.mod comm_conviqt_mod.mod comm_4d_map_mod.mod comm_zodi_mod.mod \
                           comm_tod_mapmaking_mod.mod comm_tod_bandpass_mod.mod comm_tod_pointing_mod.mod \
			   comm_tod_orbdipole_mod.mod comm_tod_gain_mod.mod comm_tod_noise_mod.mod
comm_tod_simulations_mod.o: comm_hdf_mod.mod comm_fft_mod.mod comm_shared_arr_mod.mod spline_1d_mod.mod comm_param_mod.mod comm_utils.mod
comm_tod_adc_mod.o: spline_1d_mod.mod comm_tod_mod.mod


commander.o :              comm_signal_mod.mod comm_data_mod.mod comm_cr_mod.mod comm_chisq_mod.mod \
                           comm_output_mod.mod comm_nonlin_mod.mod comm_tod_simulations_mod.mod

