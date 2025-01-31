!==============================================================================
!
!   FIRAS Regression Test facility (FRT)
!
!   MMS Revision History:
!
!    Date   	Version   SPR#   Programmer   Comments
!    ----   	-------   ----   ----------   --------
!
!  10/08/91       8.8            A.Raugh      Creation
!  10/28/91       8.8     9207   K.Jensen     Add procedures for FEX_GRTRANS
!                                             and FEX_GRTRAWWT
!  05/02/94      11.9     11741  S.Alexander  Remove procedure for FEX_GRTSWT
!
!==============================================================================

!  > Initialize suffix list and define rule for building IDL text libraries:

.suffixes

.suffixes .tlb .pro

.PRO.TLB
   $(LIBR) $(LIBRFLAGS)/TEXT $(MMS$TARGET) $(MMS$SOURCE)

!  > Create local copy of text library if it does not already exist:

.FIRST
   IF "''F$SEARCH("FRT_IDL.TLB")'" .EQS. "" -
      THEN $(LIBR)/CREATE/TEXT FRT_IDL.TLB


! **************************  FRT Facility Build  *****************************


!  > FRT procedure file list:

FRT_FILES = ave_cal_res.pro,     -
            conv_adt_sec.pro,    -
            conv_gmt_sec.pro,    -
            conv_jdt_jds.pro,    -
            conv_sec_jdt.pro,    -
            ddl_abbrev.pro,      -
            ddl_st.pro,          -
            eng_idx_plot.pro,    -
            eng_plot.pro,        -
            fdq_eng_st.pro,      -
            fdq_etr_st.pro,      -
            fdq_idx_st.pro,      -
            fdq_sdf_st.pro,      -
            fex_av_calrs_st.pro, -
            fex_englim_st.pro,   -
            fex_fakeit_st.pro,   -
            fex_gain_st.pro,     -
	    fex_grtrawwt_st.pro, -
	    fex_grttrans_st.pro, -
            fex_idx_flag_st.pro, -
            fex_idx_tols_st.pro, -
            fex_limflags_st.pro, -
            fex_mtmsweep_st.pro, -
            fex_scilim_st.pro,   -
            fex_vabsaa_st.pro,   -
            fir_fpp_verf.pro,    -
            fir_hkp_ver.pro,     -
            fpp_sdf_st.pro,      -
            get_words.pro,       -
            idx_plot.pro,        -
            nfs_anc_st.pro,      -
            nfs_hkp_st.pro,      -
            nfs_sdf_st.pro,      -
            radplot.pro,         -
            saa.pro,             -
            saa_exclude.pro,     -
            sdf_plot.pro,        -
            vab_north.pro,       -
            vab_n_exclude.pro,   -
            vab_south.pro,       -
            vab_s_exclude.pro

!   > Build the FRT_IDL text library:

FRT : frt_idl.tlb($(FRT_FILES))
   ! ***                                           ***
   ! ***   FRT procedures loaded into FRT_IDL.TLB  ***
   ! ***                                           ***
