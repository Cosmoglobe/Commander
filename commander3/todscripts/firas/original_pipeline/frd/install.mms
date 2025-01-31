! Installs facility FRD_REFERENCE_DATA.
!
! Author: R. Kummerer
!	  STX
!	  June 6, 1988
!
! Modified:
!
!         R. Kummerer, Oct 27, 1989. SER 3959. Remove FRD_HKP_STOL since
!		ENGPLOTS menu pages are hard coded.
!         R. Kummerer, Nov  9, 1989. SER 4974. Restore FRD_RPAD_APOD accidentally
!		deleted during last delivery of INSTALL file.
!         H. Wang, Mar 6, 1990. SER 3253.
!               Add FRD_CALRES for Fep_dwellplot reference data. 
!	  N. Gonzales, Apr 10, 1991. SER 7977.
!		Add FRD_AVE_CALRES for calibration resistance counts.
!	  N. Gonzales, Apr 11, 1991. SER 7980.
!		Add FRD_GFG for fakeit and gain reference data.
!	  S. Alexander, Apr 25, 1991. Add FRD_VABSAA for attitude quantities.
!	  S. Alexander, May 2, 1991. Add FRD_MTMSWEEP for sweep flyback times.
!         N. Gonzales, August 28, 1991. Add FRD_GTRANS for FDQ ref. datasets.
!         H. Wang, Dec. 9, 1991, Add FRD_MINCOADD for FEC, FSS, FIC ref. 
!	  datasets. SER 9338.
!         S. Alexander, Dec 20, 1991. Add FRD_BASIS; SER 7985.
!         S. Alexander, Dec 23, 1991. Add FRD_CMDGAIN; SER 7985.
!         S. Alexander, Dec 31, 1991. Add FRD_DGTL_TRANSIENT; SER 7985.
!         S. Alexander, Jul 13, 1992. Add FRD_NYQUIST; SER 9790.
!         S. Alexander, Aug  3, 1992. Add FRD_SAMPRATE; SER 9859.
!	  S. Alexander, Aug  6, 1992. Add FRD_REFTEMPS; SER 7985.
!         S. Alexander, Apr  5, 1993. Add FRD_APOD; SER 8836.
!         S. Alexander, Apr  5, 1993. Add FRD_EXTRACT_MODEL; SER 8178.
!         S. Alexander, May 24, 1993. Add FRD_VIBCORR; SER 8292.
!         S. Alexander, Dec 13, 1993. Add FRD_GLITCH_CORR; SER 11702.
!         S. Alexander, May  2, 1994. Remove FRD_BOLOMETER_PARMS,
!                                     FRD_GRT_SWITCH, FRD_MODEL_SWITCH,
!                                     and FRD_RPAD_APOD; SPR 11741
!         S. Brodd,     Sep  2, 1994. Add FRD_VARIANCES; SER 11409.
!         S. Brodd,     Sep 19, 1994. Add FRD_EXTRACT_PIXEL; SER 11894.
!         S. Brodd,     Sep 19, 1994. Add FRD_GALDISK; SER 11895.
!         S. Brodd,     Sep 26, 1994. Remove FRD_GALDISK; SPR 11919.
!         S. Brodd,     Nov  3, 1994. Add new IDL procedures; SER 11977.
!         S. Brodd,     Aug 28, 1995. Add FRD_ELEX_TRANSFCNL, FRD_NYQUISTL,
!                                     and FRD_APODL. SER 12244.
!         S. Brodd,     Dec 12, 1995. Add FRD_EXTRACT_MODELL; SPR 12282.
!         S. Brodd,     Jan  2, 1996. Add FRD_FLV; SPR 12286.
!         S. Brodd,     Jan  2, 1996. Add FRD_VIBCORRL; SPR 12287.
!         S. Brodd,     Aug  1, 1996. Add FRD_EXTRACT_PIXELL; SPR 12336.
!         S. Brodd,     Jun 17, 1997. Add FRD_DB_GRADIENT; SPR 12348.

.suffixes
.suffixes .olb .obj .tlb .cld .hlb .hlp .exe .pro

.cld.tlb
  $(libr) $(librflags) $(mms$target) $(mms$source)

FRD : csdr$cld:frd.cld, -
      csdr$library:csdrmsg.olb(frd_msg), -
      csdr$help:csdrhelp.hlb(frd.hlp), -
      csdr$system:frd_define_idx_params.exe, -
      csdr$system:frd_l_define_limits.exe, -
      csdr$system:frd_cc_thresholds.exe, -
      csdr$system:frd_mincoadd.exe, -
      csdr$system:frd_glitch_profile.exe, -
      csdr$system:frd_elex_transfcn.exe, -
      csdr$system:frd_elex_transfcnl.exe, -
      csdr$system:frd_dgtl_transfcn.exe, -
      csdr$system:frd_gtrans.exe, -
      csdr$system:frd_calres.exe, -
      csdr$system:frd_extrema.exe, -
      csdr$system:frd_ave_calres.exe, -
      csdr$system:frd_gfg.exe, -
      csdr$system:frd_vabsaa.exe, -
      csdr$system:frd_mtmsweep.exe, -
      csdr$system:frd_basis.exe, -
      csdr$system:frd_cmdgain.exe, -
      csdr$system:frd_dgtl_transient.exe, -
      csdr$system:frd_nyquist.exe, -
      csdr$system:frd_nyquistl.exe, -
      csdr$system:frd_samprate.exe, -
      csdr$system:frd_reftemps.exe, -
      csdr$system:frd_apod.exe, -
      csdr$system:frd_apodl.exe, -         
      csdr$system:frd_extract_model.exe, -
      csdr$system:frd_extract_modell.exe, -
      csdr$system:frd_vibcorr.exe, -
      csdr$system:frd_vibcorrl.exe, -
      csdr$system:frd_glitch_corr.exe, -
      csdr$system:frd_variances.exe, -
      csdr$system:frd_flv.exe, -
      csdr$system:frd_extract_pixel.exe, -
      csdr$system:frd_extract_pixell.exe, -
      csdr$idl:frd_cal_weights.pro, -
      csdr$idl:frd_compute_csp.pro, -
      csdr$idl:frd_compute_hr.pro, -
      csdr$idl:frd_compute_lr.pro, -
      csdr$idl:frd_fex_cvs.pro, -
      csdr$idl:frd_fex_mcs.pro, -
      csdr$idl:frd_merge_csp.pro, -
      csdr$idl:frd_db_gradient.pro
  ! Facility FRD has been installed.

csdr$cld:frd.cld : frd.cld
  COPY frd.cld csdr$cld:frd.cld;0
  ! FRD.CLD has been installed.

csdr$system:frd_define_idx_params.exe : frd_define_idx_params.exe
  COPY frd_define_idx_params.exe csdr$system:frd_define_idx_params.exe;0
  ! FRD_DEFINE_IDX_PARAMS.EXE has been installed.

csdr$system:frd_l_define_limits.exe : frd_l_define_limits.exe
  COPY frd_l_define_limits.exe csdr$system:frd_l_define_limits.exe;0
  ! FRD_L_DEFINE_LIMITS.EXE has been installed.

csdr$system:frd_cc_thresholds.exe : frd_cc_thresholds.exe
  COPY frd_cc_thresholds.exe csdr$system:frd_cc_thresholds.exe;0
  ! FRD_CC_THRESHOLDS.EXE has been installed.

csdr$system:frd_mincoadd.exe : frd_mincoadd.exe
  COPY frd_mincoadd.exe csdr$system:frd_mincoadd.exe;0
  ! FRD_MINCOADD.EXE has been installed.

csdr$system:frd_glitch_profile.exe : frd_glitch_profile.exe
  COPY frd_glitch_profile.exe csdr$system:frd_glitch_profile.exe;0
  ! FRD_GLITCH_PROFILE.EXE has been installed.

csdr$system:frd_elex_transfcn.exe : frd_elex_transfcn.exe
  COPY frd_elex_transfcn.exe csdr$system:frd_elex_transfcn.exe;0
  ! FRD_ELEX_TRANSFCN.EXE has been installed.

csdr$system:frd_elex_transfcnl.exe : frd_elex_transfcnl.exe
  COPY frd_elex_transfcnl.exe csdr$system:frd_elex_transfcnl.exe;0
  ! FRD_ELEX_TRANSFCNL.EXE has been installed.

csdr$system:frd_dgtl_transfcn.exe : frd_dgtl_transfcn.exe
  COPY frd_dgtl_transfcn.exe csdr$system:frd_dgtl_transfcn.exe;0
  ! FRD_DGTL_TRANSFCN.EXE has been installed.

csdr$system:frd_gtrans.exe : frd_gtrans.exe
  COPY frd_gtrans.exe csdr$system:frd_gtrans.exe;0
  ! FRD_GTRANS.EXE has been installed.

csdr$system:frd_calres.exe : frd_calres.exe
  COPY frd_calres.exe csdr$system:frd_calres.exe;0
  ! FRD_CALRES.EXE has been installed.

csdr$system:frd_extrema.exe : frd_extrema.exe
  COPY frd_extrema.exe csdr$system:frd_extrema.exe;0
  ! FRD_EXTREMA.EXE has been installed.

csdr$system:frd_ave_calres.exe : frd_ave_calres.exe
  COPY frd_ave_calres.exe csdr$system:frd_ave_calres.exe;0
  ! FRD_AVE_CALRES.EXE has been installed.

csdr$system:frd_gfg.exe : frd_gfg.exe
  COPY frd_gfg.exe csdr$system:frd_gfg.exe;0
  ! FRD_GFG.EXE has been installed.

csdr$system:frd_vabsaa.exe : frd_vabsaa.exe
  COPY frd_vabsaa.exe csdr$system:frd_vabsaa.exe;0
  ! FRD_VABSAA.EXE has been installed.

csdr$system:frd_mtmsweep.exe : frd_mtmsweep.exe
  COPY frd_mtmsweep.exe csdr$system:frd_mtmsweep.exe;0
  ! FRD_MTMSWEEP.EXE has been installed.            

csdr$system:frd_basis.exe : frd_basis.exe
  COPY frd_basis.exe csdr$system:frd_basis.exe;0
  ! FRD_BASIS.EXE has been installed.            

csdr$system:frd_cmdgain.exe : frd_cmdgain.exe
  COPY frd_cmdgain.exe csdr$system:frd_cmdgain.exe;0
  ! FRD_CMDGAIN.EXE has been installed.            

csdr$system:frd_dgtl_transient.exe : frd_dgtl_transient.exe
  COPY frd_dgtl_transient.exe csdr$system:frd_dgtl_transient.exe;0
  ! FRD_DGTL_TRANSIENT.EXE has been installed.            

csdr$system:frd_nyquist.exe : frd_nyquist.exe
  COPY frd_nyquist.exe csdr$system:frd_nyquist.exe;0
  ! FRD_NYQUIST.EXE has been installed.            

csdr$system:frd_nyquistl.exe : frd_nyquistl.exe
  COPY frd_nyquistl.exe csdr$system:frd_nyquistl.exe;0
  ! FRD_NYQUISTL.EXE has been installed.            

csdr$system:frd_samprate.exe : frd_samprate.exe
  COPY frd_samprate.exe csdr$system:frd_samprate.exe;0
  ! FRD_SAMPRATE.EXE has been installed.            

csdr$system:frd_reftemps.exe : frd_reftemps.exe
  COPY frd_reftemps.exe csdr$system:frd_reftemps.exe;0
  ! FRD_REFTEMPS.EXE has been installed.            

csdr$system:frd_apod.exe : frd_apod.exe
  COPY frd_apod.exe csdr$system:frd_apod.exe;0
  ! FRD_APOD.EXE has been installed.            

csdr$system:frd_apodl.exe : frd_apodl.exe
  COPY frd_apodl.exe csdr$system:frd_apodl.exe;0
  ! FRD_APODL.EXE has been installed.            

csdr$system:frd_extract_model.exe : frd_extract_model.exe
  COPY frd_extract_model.exe csdr$system:frd_extract_model.exe;0
  ! FRD_EXTRACT_MODEL.EXE has been installed.            

csdr$system:frd_extract_modell.exe : frd_extract_modell.exe
  COPY frd_extract_modell.exe csdr$system:frd_extract_modell.exe;0
  ! FRD_EXTRACT_MODELL.EXE has been installed.            

csdr$system:frd_vibcorr.exe : frd_vibcorr.exe
  COPY frd_vibcorr.exe csdr$system:frd_vibcorr.exe;0
  ! FRD_VIBCORR.EXE has been installed.            

csdr$system:frd_vibcorrl.exe : frd_vibcorrl.exe
  COPY frd_vibcorrl.exe csdr$system:frd_vibcorrl.exe;0
  ! FRD_VIBCORRL.EXE has been installed.            

csdr$system:frd_glitch_corr.exe : frd_glitch_corr.exe
  COPY frd_glitch_corr.exe csdr$system:frd_glitch_corr.exe;0
  ! FRD_GLITCH_CORR.EXE has been installed.

csdr$system:frd_variances.exe : frd_variances.exe
  COPY frd_variances.exe csdr$system:frd_variances.exe;0
  ! FRD_VARIANCES.EXE has been installed.            

csdr$system:frd_flv.exe : frd_flv.exe
  COPY frd_flv.exe csdr$system:frd_flv.exe;0
  ! FRD_FLV.EXE has been installed.            

csdr$system:frd_extract_pixel.exe : frd_extract_pixel.exe
  COPY frd_extract_pixel.exe csdr$system:frd_extract_pixel.exe;0
  ! FRD_EXTRACT_PIXEL.EXE has been installed.            

csdr$system:frd_extract_pixell.exe : frd_extract_pixell.exe
  COPY frd_extract_pixell.exe csdr$system:frd_extract_pixell.exe;0
  ! FRD_EXTRACT_PIXELL.EXE has been installed.            

csdr$idl:frd_cal_weights.pro : frd_cal_weights.pro
 copy frd_cal_weights.pro csdr$idl:frd_cal_weights.pro;0

csdr$idl:frd_compute_csp.pro : frd_compute_csp.pro
 copy frd_compute_csp.pro csdr$idl:frd_compute_csp.pro;0

csdr$idl:frd_compute_hr.pro : frd_compute_hr.pro
 copy frd_compute_hr.pro csdr$idl:frd_compute_hr.pro;0

csdr$idl:frd_compute_lr.pro : frd_compute_lr.pro
 copy frd_compute_lr.pro csdr$idl:frd_compute_lr.pro;0

csdr$idl:frd_fex_cvs.pro : frd_fex_cvs.pro
 copy frd_fex_cvs.pro csdr$idl:frd_fex_cvs.pro;0

csdr$idl:frd_fex_mcs.pro : frd_fex_mcs.pro
 copy frd_fex_mcs.pro csdr$idl:frd_fex_mcs.pro;0

csdr$idl:frd_merge_csp.pro : frd_merge_csp.pro
 copy frd_merge_csp.pro csdr$idl:frd_merge_csp.pro;0

csdr$idl:frd_db_gradient.pro : frd_db_gradient.pro
 copy frd_db_gradient.pro csdr$idl:frd_db_gradient.pro;0
