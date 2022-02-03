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
           comm_n_lcut_mod.f90 \
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
           comm_tod_gain_mod.f90 comm_tod_gain_smod.f90 \
           comm_tod_pointing_mod.f90 \
           comm_tod_lfi_mod.f90 comm_tod_lfi_smod.f90 \
           comm_tod_quiet_mod.f90 comm_tod_quiet_smod.f90 \
           comm_tod_jump_mod.f90 \
           comm_tod_spider_mod.f90 comm_tod_spider_smod.f90 \
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
           comm_tod_adc_mod.f90 comm_tod_adc_smod.f90 \
	   comm_timing_mod.f90 \
	   commander.f90

all: commander

#MODDIR := .mod
#ifneq ($(MODDIR),)
#  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
#  F90COMP+= -module $(MODDIR)
#endif

commander: $(subst .f,.o,$(FSOURCES)) $(subst .f90,.o,$(F90SOURCES))
	$(MPF90) -o commander $+ $(LINK) $(MPFCLIBS)

.PHONY: clean
clean:
	-rm -f *.o commander *.{mod,smod}

%.o %.mod %.smod: %.f90 
	$(MPF90) $(F90COMP) -c $<
	@touch $@

%.o %.mod %.smod: %.f
	$(MPF90) $(F90COMP) -c $<
	@touch $@

%.o %.mod %.smod: %.cpp
	$(MPCXX) $(CXXCOMP) -c $<
	@touch $@

ars_mod.o                   : sort_utils.o
comm_4d_map_mod.o           : comm_utils.o hashtbl_4dmap.o comm_hdf_mod.o
comm_ame_lognormal_mod.o    : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_2d_mod.o comm_data_mod.o
commander.o                 : comm_param_mod.o comm_data_mod.o comm_signal_mod.o comm_cr_mod.o comm_chisq_mod.o comm_output_mod.o comm_comp_mod.o comm_nonlin_mod.o comm_tod_simulations_mod.o comm_tod_gain_mod.o
comm_b_bl_mod.o             : comm_param_mod.o comm_map_mod.o comm_b_mod.o comm_utils.o spline_1d_mod.o
comm_beam_mod.o             : comm_b_mod.o comm_b_bl_mod.o
comm_b_mod.o                : comm_param_mod.o comm_map_mod.o spline_1d_mod.o
comm_bp_mod.o               : comm_param_mod.o comm_bp_utils.o
comm_bp_utils.o             : sort_utils.o comm_utils.o comm_hdf_mod.o
comm_chisq_mod.o            : comm_data_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_ptsrc_comp_mod.o comm_template_comp_mod.o
comm_cl_mod.o               : comm_param_mod.o comm_map_mod.o math_tools.o comm_hdf_mod.o comm_bp_utils.o invsamp_mod.o
comm_cmb_comp_mod.o         : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_f_int_0d_mod.o comm_data_mod.o comm_bp_utils.o
comm_cmb_relquad_comp_mod.o : comm_param_mod.o comm_comp_mod.o comm_template_comp_mod.o comm_f_int_mod.o comm_f_int_0d_mod.o comm_data_mod.o comm_bp_utils.o
comm_comp_mod.o             : comm_param_mod.o comm_bp_utils.o comm_bp_mod.o comm_map_mod.o comm_cr_utils.o comm_data_mod.o comm_hdf_mod.o
comm_conviqt_mod.o          : sharp.o comm_map_mod.o comm_utils.o comm_shared_arr_mod.o comm_param_mod.o
comm_cr_mod.o               : comm_comp_mod.o comm_data_mod.o comm_param_mod.o comm_diffuse_comp_mod.o comm_ptsrc_comp_mod.o comm_template_comp_mod.o comm_cr_utils.o comm_cr_precond_mod.o math_tools.o comm_output_mod.o
comm_cr_precond_mod.o       : comm_utils.o
comm_cr_utils.o             : comm_utils.o
comm_data_mod.o             : comm_param_mod.o comm_bp_mod.o comm_noise_mod.o comm_beam_mod.o comm_map_mod.o comm_tod_mod.o comm_tod_lfi_mod.o comm_tod_spider_mod.o comm_tod_wmap_mod.o comm_tod_lb_mod.o comm_tod_quiet_mod.o locate_mod.o
comm_diffuse_comp_mod.o     : comm_param_mod.o comm_comp_mod.o comm_map_mod.o comm_data_mod.o comm_f_int_mod.o comm_cl_mod.o math_tools.o comm_cr_utils.o comm_cr_precond_mod.o comm_hdf_mod.o invsamp_mod.o powell_mod.o comm_beam_mod.o
comm_diffuse_comp_smod.o    : comm_diffuse_comp_mod.o
comm_fft_mod.o              : comm_utils.o comm_param_mod.o locate_mod.o
comm_f_int_0d_mod.o         : comm_f_int_mod.o comm_comp_mod.o comm_bp_mod.o comm_utils.o
comm_f_int_1d_mod.o         : comm_f_int_mod.o comm_comp_mod.o comm_bp_mod.o spline_1d_mod.o comm_utils.o
comm_f_int_2d_mod.o         : comm_f_int_mod.o comm_comp_mod.o comm_bp_mod.o spline_2d_mod.o comm_utils.o
comm_f_int_mod.o            : comm_utils.o comm_bp_mod.o
comm_f_line_mod.o           : comm_f_int_mod.o comm_comp_mod.o comm_bp_mod.o comm_utils.o
comm_f_mod.o                : comm_comp_mod.o comm_param_mod.o comm_map_mod.o comm_bp_mod.o comm_f_int_mod.o comm_f_int_1d_mod.o
comm_freefree_comp_mod.o    : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_1d_mod.o comm_data_mod.o
comm_gain_mod.o             : comm_param_mod.o comm_data_mod.o comm_comp_mod.o
comm_hdf_mod.o              : comm_utils.o
comm_huffman_mod.o          : comm_utils.o
comm_line_comp_mod.o        : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_f_line_mod.o comm_data_mod.o comm_bp_utils.o
comm_map_mod.o              : sharp.o comm_hdf_mod.o comm_param_mod.o
comm_mbb_comp_mod.o         : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_2d_mod.o comm_data_mod.o
comm_md_comp_mod.o          : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_f_line_mod.o comm_data_mod.o comm_bp_utils.o
comm_n_mod.o                : comm_param_mod.o comm_map_mod.o
comm_noise_mod.o            : comm_n_mod.o comm_n_rms_mod.o comm_n_qucov_mod.o comm_n_lcut_mod.o 
comm_nonlin_mod.o           : comm_param_mod.o comm_data_mod.o comm_comp_mod.o comm_chisq_mod.o comm_gain_mod.o comm_line_comp_mod.o comm_diffuse_comp_mod.o comm_signal_mod.o comm_utils.o
comm_n_qucov_mod.o          : comm_n_mod.o comm_param_mod.o comm_map_mod.o
comm_n_lcut_mod.o           : comm_n_mod.o comm_param_mod.o comm_map_mod.o
comm_n_rms_mod.o            : comm_n_mod.o comm_param_mod.o comm_map_mod.o
comm_output_mod.o           : comm_param_mod.o comm_comp_mod.o comm_chisq_mod.o comm_hdf_mod.o
comm_param_mod.o            : comm_utils.o comm_status_mod.o hashtbl.o
comm_physdust_comp_mod.o    : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_1d_mod.o comm_data_mod.o spline_2d_mod.o
comm_powlaw_comp_mod.o      : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_1d_mod.o comm_data_mod.o
comm_ptsrc_comp_mod.o       : math_tools.o comm_param_mod.o comm_comp_mod.o comm_f_int_mod.o comm_f_int_0d_mod.o comm_f_int_2d_mod.o comm_data_mod.o comm_hdf_mod.o comm_cr_utils.o comm_cr_precond_mod.o locate_mod.o spline_1d_mod.o invsamp_mod.o powell_mod.o
comm_shared_arr_mod.o       : comm_utils.o
comm_shared_output_mod.o    : comm_utils.o comm_system_mod.o
comm_signal_mod.o           : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_cmb_comp_mod.o comm_cmb_relquad_comp_mod.o comm_powlaw_comp_mod.o comm_physdust_comp_mod.o comm_spindust_comp_mod.o comm_spindust2_comp_mod.o comm_ame_lognormal_mod.o comm_mbb_comp_mod.o comm_freefree_comp_mod.o comm_line_comp_mod.o comm_md_comp_mod.o comm_template_comp_mod.o comm_ptsrc_comp_mod.o comm_cr_mod.o comm_cr_utils.o comm_hdf_mod.o comm_data_mod.o comm_chisq_mod.o
comm_spindust2_comp_mod.o   : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_2d_mod.o comm_data_mod.o spline_2d_mod.o
comm_spindust_comp_mod.o    : comm_param_mod.o comm_comp_mod.o comm_diffuse_comp_mod.o comm_map_mod.o comm_f_int_1d_mod.o comm_data_mod.o spline_1d_mod.o
comm_status_mod.o           : comm_utils.o comm_system_mod.o comm_shared_output_mod.o comm_timing_mod.o
comm_task_mod.o             : comm_utils.o
comm_template_comp_mod.o    : math_tools.o comm_param_mod.o comm_comp_mod.o comm_f_int_mod.o comm_data_mod.o comm_hdf_mod.o comm_cr_utils.o comm_cr_precond_mod.o locate_mod.o spline_1d_mod.o
comm_tod_adc_mod.o          : math_tools.o spline_1d_mod.o comm_tod_mod.o comm_map_mod.o comm_param_mod.o invsamp_mod.o comm_hdf_mod.o
comm_tod_bandpass_mod.o     : comm_tod_mod.o comm_utils.o
comm_tod_driver_mod.o       : comm_tod_mod.o comm_tod_gain_mod.o comm_tod_noise_mod.o comm_tod_pointing_mod.o comm_tod_bandpass_mod.o comm_tod_orbdipole_mod.o comm_tod_simulations_mod.o comm_tod_mapmaking_mod.o comm_zodi_mod.o comm_shared_arr_mod.o
comm_tod_gain_mod.o         : comm_tod_mod.o comm_tod_noise_mod.o comm_utils.o comm_fft_mod.o math_tools.o
comm_tod_gain_smod.o        : comm_tod_gain_mod.o
comm_tod_jump_mod.o         : comm_utils.o
comm_tod_lb_mod.o           : comm_tod_mod.o comm_param_mod.o comm_map_mod.o comm_conviqt_mod.o comm_huffman_mod.o comm_hdf_mod.o comm_fft_mod.o comm_shared_arr_mod.o spline_1d_mod.o comm_4d_map_mod.o comm_tod_driver_mod.o comm_utils.o
comm_tod_lfi_mod.o          : comm_tod_mod.o comm_param_mod.o comm_map_mod.o comm_conviqt_mod.o comm_huffman_mod.o comm_hdf_mod.o comm_fft_mod.o comm_shared_arr_mod.o spline_1d_mod.o comm_4d_map_mod.o comm_tod_driver_mod.o comm_utils.o comm_tod_adc_mod.o
comm_tod_lfi_smod.o         : comm_tod_lfi_mod.o
comm_tod_quiet_mod.o				: comm_tod_mod.o comm_param_mod.o comm_map_mod.o comm_conviqt_mod.o comm_huffman_mod.o comm_hdf_mod.o comm_fft_mod.o comm_shared_arr_mod.o spline_1d_mod.o comm_4d_map_mod.o comm_tod_driver_mod.o comm_utils.o
comm_tod_quiet_smod.o				: comm_tod_quiet_mod.o
comm_tod_mapmaking_mod.o    : comm_tod_mod.o comm_utils.o comm_shared_arr_mod.o comm_map_mod.o comm_param_mod.o
comm_tod_mod.o              : comm_utils.o comm_param_mod.o comm_hdf_mod.o comm_map_mod.o comm_fft_mod.o comm_huffman_mod.o comm_conviqt_mod.o comm_zodi_mod.o comm_tod_orbdipole_mod.o comm_tod_noise_psd_mod.o
comm_tod_noise_mod.o        : comm_tod_mod.o comm_utils.o comm_fft_mod.o invsamp_mod.o comm_tod_noise_psd_mod.o
comm_tod_noise_psd_mod.o    : comm_utils.o
comm_tod_orbdipole_mod.o    : comm_utils.o comm_defs.o comm_map_mod.o spline_1d_mod.o
comm_tod_pointing_mod.o     : comm_tod_mod.o comm_map_mod.o comm_utils.o
comm_tod_simulations_mod.o  : comm_hdf_mod.o comm_fft_mod.o comm_shared_arr_mod.o spline_1d_mod.o comm_param_mod.o comm_utils.o
comm_tod_spider_mod.o       : comm_tod_mod.o comm_param_mod.o comm_map_mod.o comm_conviqt_mod.o comm_huffman_mod.o comm_hdf_mod.o comm_fft_mod.o comm_shared_arr_mod.o spline_1d_mod.o comm_4d_map_mod.o comm_zodi_mod.o comm_tod_mapmaking_mod.o comm_tod_pointing_mod.o comm_tod_gain_mod.o comm_tod_bandpass_mod.o comm_tod_orbdipole_mod.o comm_tod_jump_mod.o
comm_tod_spider_smod.o      : comm_tod_spider_mod.o
comm_tod_wmap_mod.o         : comm_tod_mod.o comm_param_mod.o comm_map_mod.o comm_conviqt_mod.o comm_huffman_mod.o comm_hdf_mod.o comm_fft_mod.o spline_1d_mod.o comm_4d_map_mod.o comm_zodi_mod.o comm_tod_mapmaking_mod.o comm_tod_pointing_mod.o comm_tod_gain_mod.o comm_tod_bandpass_mod.o comm_tod_orbdipole_mod.o comm_tod_driver_mod.o comm_utils.o
comm_utils.o                : comm_system_mod.o comm_defs.o sort_utils.o spline_1d_mod.o comm_mpi_mod.o
comm_zodi_mod.o             : comm_utils.o comm_param_mod.o
invsamp_mod.o               : spline_1d_mod.o
sharp_test.o                : comm_map_mod.o comm_param_mod.o ars_mod.o
spline_1d_mod.o             : locate_mod.o math_tools.o
spline_2d_mod.o             : spline_1d_mod.o math_tools.o comm_utils.o
comm_hdf_mod.o              : comm_utils.o
comm_timing_mod.o           : comm_utils.o 
