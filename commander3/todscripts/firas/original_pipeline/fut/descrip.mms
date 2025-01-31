! Builds the FIRAS subsystem utilities libraries, FUTLIB.OLB and FUTLIB.TLB.
!
! Author: Rob Kummerer
!	  October 14, 1986
!	  STI Incorporated
!
! Changes:
!	R. Kummerer, December 31, 1986. Remove reference to FUT_MSG in text
!		library.
!	R. Kummerer, March 12, 1987. FUT_BOL_TEMP name change to
!		FUT_TEMPERATURE_LIST.
!	R. Kummerer, April 20, 1987. Add FUT_GET_RECNUM.
!	R. Kummerer, May 7, 1987. Add FUT_TIMERANGE.
!	R. Kummerer, August 18, 1987. Add FUT_NYQUIST, FUT_FIND_IFG_CTR.
!	F. Shuman,   November 4, 1987. Add FUT_COMBINE_TEMPS.
!       R. Wilson,   November 20, 1987.  Add FUT_RDISPLAY & FUT_HIST, both
!	       being routines developed/modified by J. Dave' of ARC.
!	R. Kummerer, March 10, 1988. Add FUT_ATTITUDE, FUT_SIM_QUAT.
!	R. Kummerer, March 26, 1988. Add FUT_READ_***.
!       Shirley M. Read, May 27, 1988. Add FUT_OPEN_ATT_ARCV.
!	R. Kummerer, June 10, 1988. Add FDQ related text files.
!	R. Kummerer, June 29, 1988. Add FUT_GET_LUN, FUT_FREE_LUN.
!	R. Kummerer, July 12, 1988. Add FUT_SWG_MAKE_TEMPLATE.
!       Shirley M. Read, July 19, 1988. Add FUT_OPEN_ORBIT_ARCV.
!	R. Kummerer, Aug 25, 1988. Split FUT_DISPLAY, FUT_DISPLAY_MODEL
!		routines.
!	R. Kummerer, Sept 14, 1988. Add FUT_ANG_TO_VEC, FUT_VEC_TO_ANG.
!       Shirley M. Read, October 10, 1988. Add FUT_GOODSCI.
!	R. Kummerer, Oct 25, 1988. Remove text files FUT_FIRAS_DATA_STRUCT
!		and FUT_SCI_REC.
!       Shirley M. Read, November 1, 1988. Add FUT_ERASE.
!	R. Kummerer, November 30, 1988. SPR 2891, 2906, move FDQ_GET_QUALFLAG
!		et al to FUT.
!	R. Kummerer, January 9, 1989. SPR 3104, add FUT_DEFAULT_PEAK.
!	R. Kummerer, January 20, 1989. SPR 3100, add FUT_APOD_RECNUM.
!	R. Kummerer, March 23, 1989. SPR 3114, 3398. Add model display
!		routines.
!	R. Kummerer, April 20, 1989. SPR 3525, add FUT_VECPLT.
!	R. Kummerer, May 17, 1989. SPR 3523, add FUT_QUERY_RSE, FUT_READ_RSE.
!	R. Kummerer, July 11, 1989. SPR 4084, add FUT_Orbital_Period.
!	R. Kummerer, August 18, 1989. SER 3493, add FUT_PLOT and FUT_PLOT_LABEL,
!		the new PLT IFG/spectrum plotter.
!	R. Kummerer, September 15, 1989. SER 4567, remove FUT_HIST reference;
!		made obsolete by GLITCHMAP using FUT_PLOT.
!	R. Kummerer, September 27, 1989. SER 4568, remove FUT_RDISPLAY.
!		reference. Add FUT_SCATTER_PLOT.
!       N. Gonzales, June 6, 1990. SER 6914, add FUT_LONLAT.
!       N. Gonzales, January 9, 1991. SPR 2805, add fut_write_error.
!	S. Alexander, March 28, 1991. Name change to new FUT_TEMP_LIST from old
!		FUT_TEMPERATURE_LIST
!	N. Gonzales, April 10, 1991. Add include file FUT_CALRS.
!       N. Gonzales, April 15, 1991. Add include file FUT_GMT_CONV
!	S. Alexander, April 25, 1991. Add attitude related routines
!		FUT_ATTOFFSET and FUT_VABSAAFLAG.
!	S. Alexander, April 30, 1991. Add procedures for FDQ rewrite.
!       N. Gonzales, September 10, 1991. Add fut_qget_facr, fut_qget_faker,
!                    fut_qget_fakel, fut_qget_gainr, fut_qget_gainl,
!                    fut_qget_gtran, fut_qget_grawt
!       F. Shuman, 1991 September 27.  Add FUT_Combine_HiLo, formerly
!           FIT_Combine_Temps.
!       N. Gonzales, 1991 October 3. Add fut_qget_mtmsweep, fut_qget_vabsaa.
!       N. Gonzales, 1991 December 12. Add fut_qget_mincoadd, ref. SPR 9099.
!       N. Gonzales, 1992 Jan 1992. Add fut_get_rse, SPR 9430.
!       S. Alexander, 1992 May 13. Remove obsolete routines. SPR 9696.
!       S. Alexander, 1992 June 26. Remove more obsolete routines. SPR 9799.
!	S. Alexander, 1992 August 5. Remove routine FUT_NYQUIST. SPR 9860.
!	S. Alexander, 1993 November 29, SER 11417. Add routines
!                     FUT_FCS_AVG, FUT_FCS_MINMAX, FUT_FCS_NON_AVG, and
!                     FUT_FCS_SUM; include file FUT_FCS_INCLUDE.
!	S. Alexander, 1994 May 2, SPR 11741. Remove routines FUT_QGET_GRTSWT 
!                     and FUT_TEMPERATURE_LIST.
!       S. Brodd, HSTX, 4/95, SER 12244. Add routines FUT_APOD_ROTL and
!                     FUT_APOD_RECNUML for pipeline revision.
!       S. Brodd, HSTX, 1/96, SPR 12288. Add routine FUT_SETXAXL.
!
fut : futlib.tlb, futlib.olb
  ! Facility FUT is up to date.

! Define symbols.

text_libs = futlib.tlb/lib + csdr$library:csdrlib.tlb/lib

.for.obj
  $(fort)$(fflags)/EXTEND_SOURCE $(mms$source) + $(text_libs)

fut_object = fut_setxax.obj, -
             fut_setxaxl.obj, -
	     fut_error.obj, -
	     fut_write_error.obj, -
	     fut_gmt_conv.obj, -
	     fut_combine_hilo.obj, -
	     fut_temp_list.obj, -
	     fut_temp_correct.obj, -
	     fut_symcenter.obj, -
	     fut_default_peak.obj, -
	     fut_lonlat.obj, -
	     fut_get_rse.obj, -
	     fut_apod_recnum.obj, -
	     fut_apod_recnuml.obj, -
	     fut_apod_rotl.obj, -
	     fut_orbital_period.obj,-
	     fut_clean_catalog.obj,-
	     fut_plot_title.obj,-
	     fut_scatter_plot.obj,-
	     fut_plot.obj
fut_object_1 = fut_get_plot_parms.obj, -
	     fut_getbb.obj, -
	     fut_get_recnum.obj, -
	     fut_ave_bin_times.obj, -
	     fut_timerange.obj, -
	     fut_find_ifg_ctr.obj, -
	     fut_get_lun.obj, -
	     fut_free_lun.obj, -
	     fut_smooth_data.obj, -
	     fut_goodsci.obj, -
	     fut_erase.obj
fut_object_2 =  fut_open_att_arcv.obj, -
	     fut_open_orbit_arcv.obj, -
	     fut_attitude.obj, -
	     fut_attoffset.obj, -
	     fut_vabsaaflag.obj, -
	     fut_sim_quat.obj, -
	     fut_qget_eng.obj, -
	     fut_qget_hkp.obj, -
	     fut_qget_idx.obj, -
	     fut_qget_sci.obj, -
	     fut_qget_emf.obj, -
	     fut_qget_nos.obj
fut_object_3 = fut_qget_ssc.obj, -
	     fut_qget_anc.obj, -
	     fut_qget_ext.obj, -
	     fut_qget_cth.obj, -
	     fut_qget_facr.obj, -
	     fut_qget_faker.obj, -
	     fut_qget_fakel.obj, -
	     fut_qget_gainr.obj, -
	     fut_qget_gainl.obj, -
	     fut_qget_gtran.obj, -
	     fut_qget_grawt.obj, -
	     fut_qget_mtmsweep.obj, -
	     fut_qget_vabsaa.obj, -
	     fut_qget_mincoadd.obj, -
	     fut_qget_englim.obj, -
	     fut_qget_idxflags.obj, -
	     fut_qget_idxtols.obj, -
	     fut_qget_limflags.obj, -
	     fut_qget_scilim.obj
fut_object_4 = fut_field_attributes.obj, -
	     fut_ttinterp.obj, -
	     fut_format_time.obj, -
	     fut_swg_make_template.obj, -
	     fut_average_angles.obj, -
	     fut_get_qualflags.obj, -
	     fut_checksum.obj, -
	     fut_mjf_change.obj, -
	     fut_sum_flg.obj, -
	     fut_qual_summary.obj, -
	     fut_counts_to_ohms.obj, -
	     fut_xtalk_access.obj, -
	     fut_xtalk_check.obj, -
	     fut_trigger_plot.obj, -
	     fut_write_engplot.obj, -
	     fut_vecplt.obj, -
	     fut_query_rse.obj, -
	     fut_read_rse.obj
fut_object_5 = fut_fcs_avg.obj, -
	     fut_fcs_minmax.obj,-
	     fut_fcs_non_avg.obj,-
	     fut_fcs_sum.obj,-
             fut_fcs_reference.obj,-
	     fut_msg.obj

inc_files = fut_error.txt, -
	    fut_invoc.txt, -
	    fut_params.txt, -
	    fut_text_item.txt, -
	    fut_convbuff.txt, -
	    fut_qualflag_names.txt, -
	    fut_firnames.txt, -
	    fut_grtpars.txt, -
	    fut_lckeypars.txt, -
	    fut_qualflags.txt, -
	    fut_plot_names.txt, -
	    fut_vecplt_labels.txt, -
	    fut_calrs.txt, -
            fut_fcs_include.txt

! Source code dependencies.

fut_setxax.obj		:
fut_setxaxl.obj		: futlib.tlb(fut_params), fex_nyquistl^
fut_erase.obj		:
fut_goodsci.obj         : nfs_sdf^
fut_error.obj		: fut_error.for, futlib.tlb(fut_error)
fut_write_error.obj     : fut_write_error.for, futlib.tlb(fut_error)
fut_gmt_conv.obj        : fut_gmt_conv.for, futlib.tlb(fut_params)
fut_combine_hilo.obj	: fut_combine_hilo.for
fut_temp_list.obj	: fut_temp_list.for, futlib.tlb(fut_error) -
			   futlib.tlb(fut_params), fut_enganlg^
fut_temp_correct.obj	: fut_enganlg^
fut_symcenter.obj	: fut_symcenter.for
fut_default_peak.obj	:
fut_lonlat.obj          :
fut_get_rse.obj         :
fut_apod_recnum.obj	:
fut_apod_recnuml.obj	:
fut_apod_rotl.obj	: futlib.tlb(fut_params)
fut_get_plot_parms.obj		: futlib.tlb(fut_params)
fut_getbb.obj		: fut_getbb.for
fut_get_recnum.obj	: fut_get_recnum.for, futlib.tlb(fut_error)
fut_ave_bin_times.obj	: fut_ave_bin_times.for
fut_timerange.obj	: fut_timerange.for, futlib.tlb(fut_error)
fut_find_ifg_ctr.obj	: fut_find_ifg_ctr.for, futlib.tlb(fut_error)
fut_open_att_arcv.obj   :
fut_open_orbit_arcv.obj :
fut_attitude.obj	: fut_attitude.for, nfs_sdf^, futlib.tlb(fut_params)
fut_attoffset.obj	:
fut_vabsaaflag.obj	: fex_vabsaa^
fut_sim_quat.obj	: fut_sim_quat.for
fut_qget_eng.obj	: fdq_eng^, ct$library:ctuser.inc
fut_qget_hkp.obj	: nfs_hkp^, ct$library:ctuser.inc
fut_qget_sci.obj	: nfs_sdf^, ct$library:ctuser.inc
fut_qget_emf.obj	: nfs_emf^, ct$library:ctuser.inc
fut_qget_idx.obj	: fdq_idx^, ct$library:ctuser.inc
fut_qget_nos.obj	: fnt_noise^, ct$library:ctuser.inc
fut_qget_ssc.obj	: fec_sscal^, ct$library:ctuser.inc
fut_qget_anc.obj	: nfs_anc^, ct$library:ctuser.inc
fut_qget_ext.obj	: fxt_eng_xtrm^, ct$library:ctuser.inc
fut_qget_cth.obj	: fex_cth^, ct$library:ctuser.inc
fut_qget_facr.obj	: fex_av_calrs^, ct$library:ctuser.inc
fut_qget_faker.obj	: fex_fakeit^, ct$library:ctuser.inc
fut_qget_fakel.obj	: fex_fakeit^, ct$library:ctuser.inc
fut_qget_gainr.obj	: fex_gain^, ct$library:ctuser.inc
fut_qget_gainl.obj	: fex_gain^, ct$library:ctuser.inc
fut_qget_gtran.obj	: fex_grttrans^, ct$library:ctuser.inc
fut_qget_grawt.obj	: fex_grtrawwt^, ct$library:ctuser.inc
fut_qget_mtmsweep.obj	: fex_mtmsweep^, ct$library:ctuser.inc
fut_qget_vabsaa.obj	: fex_vabsaa^, ct$library:ctuser.inc
fut_qget_mincoadd.obj	: fex_mincoadd^, ct$library:ctuser.inc
fut_qget_englim.obj	: fex_englim^, ct$library:ctuser.inc
fut_qget_idxflags.obj	: fex_idx_flag^, ct$library:ctuser.inc
fut_qget_idxtols.obj	: fex_idx_tols^, ct$library:ctuser.inc
fut_qget_limflags.obj	: fex_limflags^, ct$library:ctuser.inc
fut_qget_scilim.obj	: fex_scilim^, ct$library:ctuser.inc
fut_field_attributes.obj : fut_field_attributes.for, futlib.tlb(fut_error)
fut_ttinterp.obj	:
fut_format_time.obj	: ct$library:ctuser.inc
fut_get_lun.obj		:
fut_free_lun.obj	:
fut_smooth_data.obj	: futlib.tlb(fut_params)
fut_swg_make_template.obj :
fut_average_angles.obj	: futlib.tlb(fut_params)
fut_get_qualflags.obj	: futlib.tlb(fut_qualflags), -
			  fdq_eng^
fut_counts_to_ohms.obj	:
fut_xtalk_access.obj	:
fut_xtalk_check.obj	:
fut_checksum.obj	:
fut_mjf_change.obj	:
fut_sum_flg.obj		: futlib.tlb(fut_qualflags), -
			  futlib.tlb(fut_params)
fut_qual_summary.obj	: futlib.tlb(fut_qualflags), -
			  futlib.tlb(fut_params)
fut_trigger_plot.obj	: futlib.tlb(fut_qualflags)
fut_write_engplot.obj	: futlib.tlb(fut_plot_names), -
			  futlib.tlb(fut_params)
fut_vecplt.obj          : futlib.tlb(fut_vecplt_labels)
fut_query_rse.obj	:
fut_read_rse.obj	: ct$library:ctuser.inc
fut_orbital_period.obj	: ct$library:ctuser.inc, -
			  futlib.tlb(fut_params), -
			  nor_orb_hdr^
fut_clean_catalog.obj	: futlib.tlb(fut_params)
fut_plot_title.obj	: futlib.tlb(fut_params)
fut_scatter_plot.obj	: futlib.tlb(fut_params)
fut_plot.obj		: futlib.tlb(fut_params)
fut_fcs_avg.obj		: futlib.tlb(fut_params), -
                          futlib.tlb(fut_fcs_include), -
                          fcs_sky^
fut_fcs_minmax.obj 	: fcf_sky^, -
                          fcs_sky^
fut_fcs_non_avg.obj 	: fcf_sky^, -
                          fcs_sky^
fut_fcs_reference.obj	: ccm_cme_catalog_entry^, -
                          fex_gltchcor^
fut_fcs_sum.obj		: futlib.tlb(fut_params), -
                          futlib.tlb(fut_fcs_include), -
                          fcf_sky^, -
                          fcs_sky^

fut_msg.obj : fut_msg.msg
  MESSAGE fut_msg


! Build the FUT facility.

futlib.tlb : $(inc_files)
  LIBRARY/TEXT/CREATE futlib $(inc_files)
  ! FUTLIB.TLB has been built.

futlib.olb : $(fut_object), $(fut_object_1), $(fut_object_2), $(fut_object_3), -
	     $(fut_object_4), $(fut_object_5), futlib.tlb
  LIBRARY/CREATE futlib $(fut_object)
  LIBRARY/REPLACE futlib $(fut_object_1)
  LIBRARY/REPLACE futlib $(fut_object_2)
  LIBRARY/REPLACE futlib $(fut_object_3)
  LIBRARY/REPLACE futlib $(fut_object_4)
  LIBRARY/REPLACE futlib $(fut_object_5)
  ! FUTLIB.OLB has been built.
