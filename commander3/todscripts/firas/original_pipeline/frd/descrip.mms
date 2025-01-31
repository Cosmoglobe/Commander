! Builds facility FRD_REFERENCE_DATA
!
! Author: R. Kummerer
!         STX
!         June 2, 1988
!
! Modified:
!         Shirley M. Read
!	  STX
!	  June 6, 1988
!	  Added FRD_Define_Idx_Params to the FRD Facility.
!         D. Bouler Feb 15, 1989. spr 2979, "configure galactic data files".
!		Added FRD_GALBEAM to .exe list. This generates a galactic
!		brightness file. 
!         R. Kummerer, April 20, 1989. SPR 3525, Add FRD_GALSPC.
!         R. Kummerer, July 21, 1989. SPR 4179. Move FRD_GALBRT and FRD_GALSPC
!		to FSF.
!         R. Kummerer, Oct 27, 1989. SER 3959. Remove FRD_HKP_STOL since
!		ENGPLOTS menu pages are hard coded.
!	  R. Kummerer, Feb 12, 1990. SPR 6018. VMS 5.2 options file changes.
!	  H. Wang, Mar 6, 1990. SPR 3253. VMS 5.8 
!                   fep_dwellplot need reference data in the archive
!         D. Bouler, Sept 28, 1990 SER 3028 
!               Add program FRD_EXTREMA to update extrema files.
!	  N. Gonzales, April 10, 1991. SER 7977. Add program FRD_AVE_CALRES.
!	  N. Gonzales, April 11, 1991. SER 7980. Add program FRD_GFG.
!	  S. Alexander, April 25, 1991. Add program FRD_VABSAA.
!	  S. Alexander, May 2, 1991. Add program FRD_MTMSWEEP.
!         N. Gonzales, August 28, 1991. Add program FRD_GTRANS
!         H. Wang, Dec. 9, 1991. SER 9338. Add program FRD_MINCOADD
!         S. Alexander, Dec 20, 1991. SER 7985. Add program FRD_BASIS
!         S. Alexander, Dec 23, 1991. SER 7985. Add program FRD_CMDGAIN
!         S. Alexander, Dec 30, 1991. SER 7985. Add IMSL$DIR to link statement 
!            for FRD_GLITCH_PROFILE.
!         S. Alexander, Dec 31, 1991. SER 7985. Add program FRD_DGTL_TRANSIENT
!         S. Alexander, Jul 13, 1992. Add program FRD_NYQUIST; SER 9790.
!	  S. Alexander, Aug  3, 1992. Add program FRD_SAMPRATE; SER 9859.
!	  S. Alexander, Aug  6, 1992. Add program FRD_REFTEMPS; SER 7985.
!	  S. Alexander, Apr  5, 1993. Add program FRD_APOD; SER 8836.
!         S. Alexander, Apr  5, 1993. Add program FRD_EXTRACT_MODEL; SER 8178.
!	  S. Alexander, May 24, 1993. Add program FRD_VIBCORR; SER 8292.
!	  S. Alexander, Jun 16, 1993. Add IMSL library to link statement of
!                                     FRD_GLITCH_PROFILE for IMSL version 2;
!                                     SER 10601.
!	  S. Alexander, Dec 13, 1993. Add program FRD_GLITCH_CORR. SER 11702.
!	  S. Alexander, May  2, 1994. Remove programs FRD_BOLOMETER_PARMS,
!                                     FRD_GRT_SWITCH, FRD_MODEL_SWITCH, and
!                                     FRD_RPAD_APOD. SPR 11741.
!	  S. Brodd,     Sep  2, 1994. Add program FRD_VARIANCES. SER 11409.
!	  S. Brodd,     Sep 19, 1994. Add program FRD_EXTRACT_PIXEL. SER 11894.
!	  S. Brodd,     Aug 28, 1995. Add programs FRD_ELEX_TRANSFCNL,
!                                     FRD_NYQUISTL, and FRD_APODL. SER 12244.
!	  S. Brodd,     Dec 19, 1995. Add program FRD_EXTRACT_MODELL. SPR 12282.
!	  S. Brodd,     Jan  2, 1996. Add program FRD_FLV SPR 12286.
!	  S. Brodd,     Jan  2, 1996. Add program FRD_VIBCORRL SPR 12287.
!	  S. Brodd,     Aug  1, 1996. Add program FRD_EXTRACT_PIXELL SPR 12336.

.suffixes
.suffixes .exe .olb .obj .for .cld .msg .tlb .txt

FRD : FRDBLD.TLB, FRD_L_DEFINE_LIMITS.EXE, -
      FRD_DEFINE_IDX_PARAMS.EXE, FRD_CC_THRESHOLDS.EXE, FRD_MINCOADD.EXE, -
      FRD_GLITCH_PROFILE.EXE, FRD_ELEX_TRANSFCN.EXE, FRD_DGTL_TRANSFCN.EXE, -
      FRD_CALRES.EXE, FRD_EXTREMA.EXE, FRD_AVE_CALRES.EXE, -
      FRD_GFG.EXE, FRD_VABSAA.EXE, FRD_MTMSWEEP.EXE, -
      FRD_GTRANS.EXE, FRD_BASIS.EXE, FRD_CMDGAIN.EXE, -
      FRD_DGTL_TRANSIENT.EXE, FRD_NYQUIST.EXE, -
      FRD_SAMPRATE.EXE, FRD_REFTEMPS.EXE, FRD_APOD.EXE, FRD_EXTRACT_MODEL.EXE, -
      FRD_VIBCORR.EXE, FRD_GLITCH_CORR.EXE, FRD_VARIANCES.EXE, -
      FRD_EXTRACT_PIXEL.EXE, FRD_ELEX_TRANSFCNL.EXE, FRD_NYQUISTL.EXE, -
      FRD_APODL.EXE, FRD_VIBCORRL.EXE, FRD_EXTRACT_MODELL.EXE, FRD_FLV.EXE, -
      FRD_EXTRACT_PIXELL.EXE
  !Facility FRD is up to date.

text_libs = FRDBLD.TLB/LIB + CSDR$LIBRARY:FUTLIB.TLB/LIB + -
            CSDR$LIBRARY:CSDRLIB.TLB/LIB
            
inc_files = FRD_MODEL_INVOC.TXT, -
            FRD_MODEL_SOLN.TXT, -
            FRD_MODEL_SOLNL.TXT, -
            FRD_MODEL_CONFIG.TXT, -
            FRD_MODEL_CONFIGL.TXT, -
            FRD_INVOC_VARIANCES.TXT

.for.obj
 $(fort) /EXTEND_SOURCE $(fflags) $(mms$source) + $(text_libs)

!Message file

frd_msg.obj : frd_msg.msg
  MESSAGE frd_msg

!Text library

frdbld.tlb : $(inc_files)
  LIBRARY/TEXT/CREATE frdbld $(inc_files)
  ! FRDBLD.TLB has been built


!Build FRD_L_DEFINE_LIMITS

frd_l_define_limits_obj = frd_l_define_limits.obj, -
			  frd_l_extract_flags.obj, -
			  frd_l_extract_gparms.obj, -
			  frd_l_extract_iparms.obj, -
			  frd_l_extract_rparms.obj, -
			  frd_l_process_eng.obj, -
			  frd_l_process_flags.obj, -
			  frd_l_process_sci.obj, -
			  frd_l_remtab.obj, -
			  frd_l_store_eng.obj, -
			  frd_l_store_flags.obj, -
			  frd_l_store_grt.obj, -
			  frd_l_store_sci.obj, -
			  frd_l_write_eng.obj, -
			  frd_l_write_flags.obj, -
			  frd_l_write_sci.obj, -
		          frd_msg.obj

frd_l_define_limits.exe :	frdbld.olb( $(frd_l_define_limits_obj) ), -
				csdr$library:futlib.olb, -
				csdr$library:futlib.tlb, -
				csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_l_define_limits), -
			csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:csdrlib/lib, -
			csdr$library:v5.opt/option
 ! FRD_L_DEFINE_LIMITS.EXE has been built

!Build FRD_DEFINE_IDX_PARAMS

frd_define_idx_params_obj = frd_define_idx_params.obj, -
			    frd_parse_idx_flags.obj, -
			    frd_parse_idx_tols.obj, -
			    frd_extract_idx_flagparms.obj, -
			    frd_extract_idx_tolparms.obj, -
			    frd_l_remtab.obj

frd_define_idx_params.exe :  frdbld.olb($(frd_define_idx_params_obj)), -
			     csdr$library:v5.opt
 $(link) $(linkflags)	frdbld.olb/lib/inc=(frd_define_idx_params), -
			csdr$library:v5.opt/option
 !FRD_DEFINE_IDX_PARAMS has been built.

!Build FRD_CC_THRESHOLDS

frd_cc_thresholds_obj = frd_cc_thresholds.obj,-
		          frd_msg.obj

frd_cc_thresholds.exe :	frdbld.olb( $(frd_cc_thresholds_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_cc_thresholds), -
			csdr$library:v5.opt/option
 ! FRD_CC_THRESHOLDS.EXE has been built

!Build FRD_MINCOADD

frd_mincoadd_obj = frd_mincoadd.obj,-
		          frd_msg.obj

frd_mincoadd.exe :	frdbld.olb( $(frd_mincoadd_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_mincoadd), -
			csdr$library:v5.opt/option
 ! FRD_MINCOADD.EXE has been built

!Build FRD_GLITCH_PROFILE

frd_glitch_profile_obj = frd_glitch_profile.obj,-
		          frd_msg.obj

frd_glitch_profile.exe :	frdbld.olb( $(frd_glitch_profile_obj) ), -
                                imsl$dir:imsl.olb, -
                                sys$library:vaxcrtl.olb, -
				csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_glitch_profile), -
                        imsl$dir:imsl.olb/lib, -
                        sys$library:vaxcrtl.olb/lib, -
			csdr$library:v5.opt/option
 ! FRD_GLITCH_PROFILE.EXE has been built

!Build FRD_ELEX_TRANSFCN

frd_elex_transfcn_obj = frd_elex_transfcn.obj,-
		          frd_msg.obj
frd_elex_transfcn.exe :	frdbld.olb( $(frd_elex_transfcn_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_elex_transfcn), -
			csdr$library:v5.opt/option
 ! FRD_ELEX_TRANSFCN.EXE has been built

!Build FRD_ELEX_TRANSFCNL

frd_elex_transfcnl_obj = frd_elex_transfcnl.obj,-
		          frd_msg.obj
frd_elex_transfcnl.exe :	frdbld.olb( $(frd_elex_transfcnl_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_elex_transfcnl), -
			csdr$library:v5.opt/option
 ! FRD_ELEX_TRANSFCNL.EXE has been built

!Build FRD_DGTL_TRANSFCN

frd_dgtl_transfcn_obj = frd_dgtl_transfcn.obj,-
		          frd_msg.obj
frd_dgtl_transfcn.exe :	frdbld.olb( $(frd_dgtl_transfcn_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_dgtl_transfcn), -
			csdr$library:v5.opt/option
 ! FRD_DGTL_TRANSFCN.EXE has been built

!Build FRD_GTRANS

frd_gtrans_obj = frd_gtrans.obj,-
		          frd_msg.obj
frd_gtrans.exe :	frdbld.olb( $(frd_gtrans_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_gtrans), -
			csdr$library:v5.opt/option
 ! FRD_GTRANS.EXE has been built

!Build FRD_CALRES

frd_calres_obj = frd_calres.obj,-
		          frd_msg.obj
frd_calres.exe :	frdbld.olb( $(frd_calres_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_calres), -
			csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_CALRES.EXE has been built

!Build FRD_EXTREMA

frd_extrema_obj =       frd_extrema.obj,-
		        frd_msg.obj
frd_extrema.exe :	frdbld.olb( $(frd_extrema_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_extrema), -
			csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_EXTREMA.EXE has been built

!Build FRD_AVE_CALRES

frd_ave_calres_obj = 	frd_ave_calres.obj,-
			frd_facr_get_options.obj,-
			frd_read_sum_avg.obj,-
			frd_msg.obj

frd_ave_calres.exe :	frdbld.olb( $(frd_ave_calres_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_ave_calres), -
			csdr$library:futlib.olb/lib/inc=(fut_error), -
                        csdr$library:csdrlib.olb/lib,-
			csdr$library:v5.opt/option
 ! FRD_AVE_CALRES has been built

!Build FRD_GFG

frd_gfg_obj = 	frd_gfg.obj,-
		frd_fakeit.obj,-
		frd_gain.obj,-
		frd_gfg_check.obj,-
		frd_gfg_close.obj,-
		frd_gfg_get_options.obj,-
		frd_gfg_open.obj,-
		frd_gfg_read2.obj,-
		frd_msg.obj

frd_gfg.exe :	frdbld.olb( $(frd_gfg_obj) ), -
		csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_gfg), -
			csdr$library:futlib.olb/lib/inc=(fut_error), -
                        csdr$library:csdrlib.olb/lib,-
			csdr$library:v5.opt/option
 ! FRD_GFG has been built

!Build FRD_VABSAA

frd_vabsaa_obj = frd_vabsaa.obj,-
		          frd_msg.obj

frd_vabsaa.exe :	frdbld.olb( $(frd_vabsaa_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_vabsaa), -
			csdr$library:v5.opt/option
 ! FRD_VABSAA.EXE has been built

!Build FRD_MTMSWEEP

frd_mtmsweep_obj = frd_mtmsweep.obj,-
		          frd_msg.obj

frd_mtmsweep.exe :	frdbld.olb( $(frd_mtmsweep_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_mtmsweep), -
			csdr$library:v5.opt/option
 ! FRD_MTMSWEEP.EXE has been built                                  

!Build FRD_BASIS

frd_basis_obj = frd_basis.obj,-
		          frd_msg.obj

frd_basis.exe	 :	frdbld.olb( $(frd_basis_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_basis), -
			csdr$library:v5.opt/option
 ! FRD_BASIS.EXE has been built                                  

!Build FRD_CMDGAIN

frd_cmdgain_obj = frd_cmdgain.obj,-
		          frd_msg.obj

frd_cmdgain.exe	 :	frdbld.olb( $(frd_cmdgain_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_cmdgain), -
			csdr$library:v5.opt/option
 ! FRD_CMDGAIN.EXE has been built                                  

!Build FRD_DGTL_TRANSIENT

frd_dgtl_transient_obj = frd_dgtl_transient.obj,-
		         frd_dgtl_transient_compress.obj,-
		         frd_dgtl_transient_fcn.obj,-
		         frd_msg.obj
frd_dgtl_transient.exe : frdbld.olb( $(frd_dgtl_transient_obj) ), -
			 csdr$library:futlib.olb, -
			 csdr$library:v5.opt
 $(link) $(linkflags)	 frdbld/lib/inc=(frd_dgtl_transient), -
			 csdr$library:futlib.olb/lib/inc=(fut_msg), -
			 csdr$library:v5.opt/option
 ! FRD_DGTL_TRANSIENT.EXE has been built

!Build FRD_NYQUIST

frd_nyquist_obj = frd_nyquist.obj,-
		  frd_msg.obj

frd_nyquist.exe	 :	frdbld.olb( $(frd_nyquist_obj) ), -
                        csdr$library:futlib.olb, -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_nyquist), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_NYQUIST.EXE has been built                                  

!Build FRD_NYQUISTL

frd_nyquistl_obj = frd_nyquistl.obj,-
		  frd_msg.obj

frd_nyquistl.exe	 :	frdbld.olb( $(frd_nyquistl_obj) ), -
                        csdr$library:futlib.olb, -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_nyquistl), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_NYQUISTL.EXE has been built                                  

!Build FRD_SAMPRATE

frd_samprate_obj = frd_samprate.obj,-
	           frd_msg.obj
frd_samprate.exe :	frdbld.olb( $(frd_samprate_obj) ), -
			csdr$library:futlib.olb, -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_samprate), -
			csdr$library:futlib.olb/lib, -
			csdr$library:v5.opt/option
 ! FRD_SAMPRATE.EXE has been built

!Build FRD_REFTEMPS

frd_reftemps_obj = frd_reftemps.obj, -
                   frd_msg.obj

frd_reftemps.exe : frdbld.olb( $(frd_reftemps_obj) ), -
                   csdr$library:v5.opt
 $(link) $(linkflags) frdbld/lib/inc=(frd_reftemps), -
                      csdr$library:v5.opt/option
 ! FRD_REFTEMPS has been built

!Build FRD_APOD

frd_apod_obj = frd_apod.obj,-
               frd_apod_fcn.obj,-
               frd_msg.obj

frd_apod.exe :	frdbld.olb( $(frd_apod_obj) ), -
                csdr$library:futlib.olb, -
	       	csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_apod), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_APOD.EXE has been built

!Build FRD_APODL

frd_apodl_obj = frd_apodl.obj,-
               frd_apod_fcnl.obj,-
               frd_msg.obj

frd_apodl.exe :	frdbld.olb( $(frd_apodl_obj) ), -
                csdr$library:futlib.olb, -
	       	csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_apodl), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_APODL.EXE has been built

!Build FRD_EXTRACT_MODEL

frd_extract_model_obj = frd_extract_model.obj,-
			frd_parse_model.obj,-
			frd_read_model.obj,-
			frd_write_model.obj,-
			frd_msg.obj

frd_extract_model.exe :	frdbld.tlb, -
                        frdbld.olb( $(frd_extract_model_obj) ), -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_extract_model), -
			csdr$library:futlib.olb/lib, -
                        csdr$library:csdrlib.olb/lib,-
			csdr$library:v5.opt/option
 ! FRD_EXTRACT_MODEL has been built

!Build FRD_EXTRACT_MODELL

frd_extract_modell_obj = frd_extract_modell.obj,-
 			 frd_parse_modell.obj,-
  			 frd_read_modell.obj,-
			 frd_write_modell.obj,-
			 frd_msg.obj

frd_extract_modell.exe : frdbld.tlb, -
                         frdbld.olb( $(frd_extract_modell_obj) ), -
			 csdr$library:v5.opt
 $(link) $(linkflags)	 frdbld/lib/inc=(frd_extract_modell), -
			 csdr$library:futlib.olb/lib, -
                         csdr$library:csdrlib.olb/lib,-
			 csdr$library:v5.opt/option
 ! FRD_EXTRACT_MODELL has been built

!Build FRD_VIBCORR

frd_vibcorr_obj = frd_vibcorr.obj,-
		  frd_msg.obj

frd_vibcorr.exe	 :	frdbld.olb( $(frd_vibcorr_obj) ), -
                        csdr$library:futlib.olb, -
			csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_vibcorr), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_VIBCORR.EXE has been built                                  

!Build FRD_VIBCORRL

frd_vibcorrl_obj = frd_vibcorrl.obj,-
 		   frd_msg.obj

frd_vibcorrl.exe : frdbld.olb( $(frd_vibcorrl_obj) ), -
                           csdr$library:futlib.olb, -
			   csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_vibcorrl), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_VIBCORRL.EXE has been built                                  

!Build FRD_VARIANCES

frd_variances_obj = frd_variances.obj,-
                    frd_parse_variances.obj,-
                    frd_read_variances.obj,-
                    frd_galactic_cut.obj,-
                    frd_msg.obj

frd_variances.exe : frdbld.olb($(frd_variances_obj)), -
                    csdr$library:futlib.olb, -
                    csdr$library:v5.opt
 $(link) $(linkflags) frdbld/lib/inc=(frd_variances), -
                      csdr$library:futlib.olb/lib/inc=(fut_msg), -
                      csdr$library:v5.opt/option
 ! FRD_VARIANCES.EXE has been built                                  

!Build FRD_FLV

frd_flv_obj = frd_flv.obj,-
              frd_parse_flv.obj,-
              frd_read_flv.obj,-
              frd_gal_cut_flv.obj,-
              frd_msg.obj

frd_flv.exe : frdbld.olb($(frd_flv_obj)), -
              csdr$library:futlib.olb, -
              csdr$library:v5.opt
 $(link) $(linkflags) frdbld/lib/inc=(frd_flv), -
                      csdr$library:futlib.olb/lib/inc=(fut_msg), -
                      csdr$library:v5.opt/option
 ! FRD_FLV.EXE has been built                                  

!Build FRD_GLITCH_CORR

frd_glitch_corr_obj = frd_glitch_corr.obj,-
		      frd_msg.obj

frd_glitch_corr.exe : frdbld.olb( $(frd_glitch_corr_obj) ), -
                      csdr$library:futlib.olb, -
		      csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_glitch_corr), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
			csdr$library:v5.opt/option
 ! FRD_GLITCH_CORR.EXE has been built

!Build FRD_EXTRACT_PIXEL

frd_extract_pixel_obj = frd_extract_pixel.obj,-
 		        frd_msg.obj

frd_extract_pixel.exe : frdbld.olb( $(frd_extract_pixel_obj) ), -
                        csdr$library:futlib.olb, -
    		        csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_extract_pixel), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
                        csdr$library:csdrlib.olb/lib,-
			csdr$library:v5.opt/option
 ! FRD_EXTRACT_PIXEL.EXE has been built

!Build FRD_EXTRACT_PIXELL

frd_extract_pixell_obj = frd_extract_pixell.obj,-
  		         frd_msg.obj

frd_extract_pixell.exe : frdbld.olb( $(frd_extract_pixell_obj) ), -
                        csdr$library:futlib.olb, -
    		        csdr$library:v5.opt
 $(link) $(linkflags)	frdbld/lib/inc=(frd_extract_pixell), -
                        csdr$library:futlib.olb/lib/inc=(fut_msg), -
                        csdr$library:csdrlib.olb/lib,-
			csdr$library:v5.opt/option
 ! FRD_EXTRACT_PIXELL.EXE has been built

!Dependencies

frd_l_define_limits.obj		: csdr$library:futlib.tlb
frd_l_extract_flags.obj		: csdr$library:futlib.tlb
frd_l_extract_gparms.obj	: csdr$library:futlib.tlb
frd_l_extract_iparms.obj	: csdr$library:futlib.tlb
frd_l_extract_rparms.obj	: csdr$library:futlib.tlb
frd_l_process_eng.obj		: csdr$library:futlib.tlb
frd_l_process_flags.obj		: csdr$library:futlib.tlb
frd_l_process_sci.obj		: csdr$library:futlib.tlb
frd_l_remtab.obj		: csdr$library:futlib.tlb
frd_l_store_eng.obj		: csdr$library:futlib.tlb
frd_l_store_flags.obj		: csdr$library:futlib.tlb
frd_l_store_sci.obj		: csdr$library:futlib.tlb
frd_l_write_eng.obj		: csdr$library:futlib.tlb
frd_l_write_flags.obj		: csdr$library:futlib.tlb
frd_l_write_sci.obj		: csdr$library:futlib.tlb
frd_define_idx_params.obj       : csdr$library:futlib.tlb
frd_parse_idx_flags.obj		: csdr$library:futlib.tlb, fex_idx_flag^
frd_parse_idx_tols.obj		: fex_idx_tols^
frd_extract_idx_flagparms.obj	: 
frd_extract_idx_tolparms.obj	: 
frd_cc_thresholds.obj		: fex_cth^
frd_mincoadd.obj		: 
frd_elex_transfcn.obj		: 
frd_elex_transfcnl.obj		: 
frd_dgtl_transfcn.obj		: 
frd_gtrans.obj                  : 
frd_calres.obj	        	: 
frd_ave_calres.obj		: csdr$library:futlib.tlb
frd_facr_get_options.obj	: csdr$library:csdrlib.tlb, -
				  csdr$library:futlib.tlb
frd_read_sum_avg.obj		: csdr$library:futlib.tlb
frd_gfg.obj			: csdr$library:futlib.tlb, -
				  csdr$library:csdrlib.tlb
frd_fakeit.obj 			: nfs_hkp^, fex_fakeit^
frd_gain.obj			: nfs_hkp^, fex_gain^
frd_gfg_check.obj		: csdr$library:futlib.tlb, -
				  csdr$library:csdrlib.tlb, -
				  nfs_hkp^, fex_fakeit^, fex_gain^
frd_gfg_close.obj               : csdr$library:futlib.tlb, -
				  csdr$library:csdrlib.tlb
frd_gfg_get_options.obj         : csdr$library:csdrlib.tlb
frd_gfg_open.obj                : csdr$library:futlib.tlb, -
				  csdr$library:csdrlib.tlb
frd_gfg_read2.obj		: csdr$library:csdrlib.tlb, -
				  nfs_hkp^
frd_vabsaa.obj			: fex_vabsaa^
frd_mtmsweep.obj		: fex_mtmsweep^
frd_basis.obj	      		: fex_basis^
frd_cmdgain.obj	      		: fex_cmdgain^
frd_dgtl_transient.obj		:
frd_dgtl_transient_compress.obj	:
frd_dgtl_transient_fcn.obj	:
frd_nyquist.obj 		: fex_nyquist^
frd_nyquistl.obj 		: fex_nyquistl^
frd_samprate.obj 		: fex_samprate^
frd_reftemps.obj		: fex_reftemps^
frd_apod.obj 			:
frd_apod_fcn.obj 		:
frd_apodl.obj 			:
frd_apod_fcnl.obj 		:
frd_extract_model.obj		: frdbld.tlb
frd_extract_modell.obj		: frdbld.tlb, -
                                  csdr$library:ctuser.inc
frd_parse_model.obj		: frdbld.tlb, -
                                  csdr$library:futlib.tlb, -
                                  csdr$library:csdrlib.tlb
frd_parse_modell.obj		: frdbld.tlb, -
                                  csdr$library:futlib.tlb, -
                                  csdr$library:csdrlib.tlb
frd_read_model.obj 		: frdbld.tlb, -
                                  csdr$library:futlib.tlb, -
                                  fex_mod^
frd_read_modell.obj 		: frdbld.tlb, -
                                  csdr$library:futlib.tlb
frd_write_model.obj 		: frdbld.tlb, -
                                  csdr$library:futlib.tlb
frd_write_modell.obj 		: frdbld.tlb, -
                                  csdr$library:futlib.tlb
frd_vibcorr.obj 		: fex_vibcorr^
frd_vibcorrl.obj 		: csdr$library:futlib.tlb, -
                                  fex_vibcorrl^
frd_glitch_corr.obj		: fex_gltchcor^
frd_variances.obj		: csdr$library:futlib.tlb, -
                                  fex_var^
frd_parse_variances.obj   	: csdr$library:futlib.tlb
frd_read_variances.obj 		: csdr$library:ctuser.inc, -
                                  csdr$library:futlib.tlb, -
                                  fcf_sky^
frd_galactic_cut.obj 		: csdr$library:futlib.tlb
frd_flv.obj	          	: csdr$library:ctuser.inc, -
                                  csdr$library:futlib.tlb, -
                                  fex_flv^
frd_parse_flv.obj       	: csdr$library:futlib.tlb
frd_read_flv.obj 		: csdr$library:futlib.tlb, -
                                  fil_sky^
frd_gal_cut_flv.obj 		: csdr$library:futlib.tlb
frd_extract_pixel.obj		: csdr$library:ctuser.inc, -
                                  csdr$library:futlib.tlb, -
                                  fic_sky^
frd_extract_pixell.obj		: csdr$library:ctuser.inc, -
                                  csdr$library:futlib.tlb, -
                                  fil_sky^
