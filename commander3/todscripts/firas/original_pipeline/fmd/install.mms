!
!   Purpose: Install facility fmd_megadestriper.
!
!   Author: K. Jensen, HSTX, 5/02/97, SPR 12351
!
!   Modifications: 
!

.suffixes
.suffixes .pro

.default
      copy/log $(mms$source) $(mms$target);0

fmd : csdr$idl:fmd_adt_sec.pro,-
      csdr$idl:fmd_apply.pro,-
      csdr$idl:fmd_avc_hif2.pro,-
      csdr$idl:fmd_avc_hif3.pro,-
      csdr$idl:fmd_avc_hif4.pro,-
      csdr$idl:fmd_avc_high.pro,-
      csdr$idl:fmd_avc_hres.pro,-
      csdr$idl:fmd_avc_lowf.pro,-
      csdr$idl:fmd_avc_st.pro,-
      csdr$idl:fmd_avec_high.pro,-
      csdr$idl:fmd_avec_high_2.pro,-
      csdr$idl:fmd_avec_high_3.pro,-
      csdr$idl:fmd_avec_high_4.pro,-
      csdr$idl:fmd_avec_hres.pro,-
      csdr$idl:fmd_avec_lowf.pro,-
      csdr$idl:fmd_bad2_hi.pro,-
      csdr$idl:fmd_bad2_hr.pro,-
      csdr$idl:fmd_bad2_lhf.pro,-
      csdr$idl:fmd_bad2_lhs.pro,-
      csdr$idl:fmd_bad2_llf.pro,-
      csdr$idl:fmd_bad2_lls.pro,-
      csdr$idl:fmd_bad2_lo.pro,-
      csdr$idl:fmd_bad2_lsf.pro,-
      csdr$idl:fmd_bad2_rhf.pro,-
      csdr$idl:fmd_bad2_rhs.pro,-
      csdr$idl:fmd_bad2_rlf.pro,-
      csdr$idl:fmd_bad2_rls.pro,-
      csdr$idl:fmd_bad2_rsf.pro,-
      csdr$idl:fmd_bad_hi.pro,-
      csdr$idl:fmd_bad_hr.pro,-
      csdr$idl:fmd_bad_lhf.pro,-
      csdr$idl:fmd_bad_lhs.pro,-
      csdr$idl:fmd_bad_llf.pro,-
      csdr$idl:fmd_bad_lls.pro,-
      csdr$idl:fmd_bad_lo.pro,-
      csdr$idl:fmd_bad_lsf.pro,-
      csdr$idl:fmd_bad_rhf.pro,-
      csdr$idl:fmd_bad_rhs.pro,-
      csdr$idl:fmd_bad_rlf.pro,-
      csdr$idl:fmd_bad_rls.pro,-
      csdr$idl:fmd_bad_rsf.pro,-
      csdr$idl:fmd_bandspec.pro,-
      csdr$idl:fmd_cal_chi2.pro,-
      csdr$idl:fmd_chi2.pro,-
      csdr$idl:fmd_chi2_high_2.pro,-
      csdr$idl:fmd_chi2_high_3.pro,-
      csdr$idl:fmd_chi2_high_4.pro,-
      csdr$idl:fmd_chi2_hres.pro,-
      csdr$idl:fmd_chi2_lhf_2.pro,-
      csdr$idl:fmd_chi2_lhf_2x.pro,-
      csdr$idl:fmd_chi2_lhf_3.pro,-
      csdr$idl:fmd_chi2_lhf_3x.pro,-
      csdr$idl:fmd_chi2_lhf_4.pro,-
      csdr$idl:fmd_chi2_lhf_4x.pro,-
      csdr$idl:fmd_chi2_lhs_2.pro,-
      csdr$idl:fmd_chi2_lhs_2x.pro,-
      csdr$idl:fmd_chi2_lhs_3.pro,-
      csdr$idl:fmd_chi2_lhs_3x.pro,-
      csdr$idl:fmd_chi2_lhs_4.pro,-
      csdr$idl:fmd_chi2_lhs_4x.pro,-
      csdr$idl:fmd_chi2_llf.pro,-
      csdr$idl:fmd_chi2_llfx.pro,-
      csdr$idl:fmd_chi2_lls.pro,-
      csdr$idl:fmd_chi2_llsx.pro,-
      csdr$idl:fmd_chi2_lowf.pro,-
      csdr$idl:fmd_chi2_lsf.pro,-
      csdr$idl:fmd_chi2_lsfx.pro,-
      csdr$idl:fmd_chi2_rhf_2.pro,-
      csdr$idl:fmd_chi2_rhf_2x.pro,-
      csdr$idl:fmd_chi2_rhf_3.pro,-
      csdr$idl:fmd_chi2_rhf_3x.pro,-
      csdr$idl:fmd_chi2_rhf_4.pro,-
      csdr$idl:fmd_chi2_rhf_4x.pro,-
      csdr$idl:fmd_chi2_rhs_2.pro,-
      csdr$idl:fmd_chi2_rhs_2x.pro,-
      csdr$idl:fmd_chi2_rhs_3.pro,-
      csdr$idl:fmd_chi2_rhs_3x.pro,-
      csdr$idl:fmd_chi2_rhs_4.pro,-
      csdr$idl:fmd_chi2_rhs_4x.pro,-
      csdr$idl:fmd_chi2_rlf.pro,-
      csdr$idl:fmd_chi2_rlfx.pro,-
      csdr$idl:fmd_chi2_rls.pro,-
      csdr$idl:fmd_chi2_rlsx.pro,-
      csdr$idl:fmd_chi2_rsf.pro,-
      csdr$idl:fmd_chi2_rsfx.pro,-
      csdr$idl:fmd_choleskey.pro,-
      csdr$idl:fmd_concat2_hi.pro,-
      csdr$idl:fmd_concat2_hr.pro,-
      csdr$idl:fmd_concat_high.pro,-
      csdr$idl:fmd_concat_hres.pro,-
      csdr$idl:fmd_concat_lowf.pro,-
      csdr$idl:fmd_covar.pro,-
      csdr$idl:fmd_covar_high.pro,-
      csdr$idl:fmd_covar_high2.pro,-
      csdr$idl:fmd_covar_high3.pro,-
      csdr$idl:fmd_covar_high4.pro,-
      csdr$idl:fmd_covar_hi_23.pro,-
      csdr$idl:fmd_covar_hi_24.pro,-
      csdr$idl:fmd_covar_hi_34.pro,-
      csdr$idl:fmd_covar_hres.pro,-
      csdr$idl:fmd_covar_hres1.pro,-
      csdr$idl:fmd_covar_hres2.pro,-
      csdr$idl:fmd_covar_hres3.pro,-
      csdr$idl:fmd_covar_hr_12.pro,-
      csdr$idl:fmd_covar_hr_13.pro,-
      csdr$idl:fmd_covar_hr_23.pro,-
      csdr$idl:fmd_covar_lowf.pro,-
      csdr$idl:fmd_covar_off.pro,-
      csdr$idl:fmd_cov_hif2.pro,-
      csdr$idl:fmd_cov_hif3.pro,-
      csdr$idl:fmd_cov_hif4.pro,-
      csdr$idl:fmd_cov_high.pro,-
      csdr$idl:fmd_cov_hres.pro,-
      csdr$idl:fmd_cov_lowf.pro,-
      csdr$idl:fmd_cov_st.pro,-
      csdr$idl:fmd_csq_hif2.pro,-
      csdr$idl:fmd_csq_hif3.pro,-
      csdr$idl:fmd_csq_hif4.pro,-
      csdr$idl:fmd_csq_high.pro,-
      csdr$idl:fmd_csq_hres.pro,-
      csdr$idl:fmd_csq_lowf.pro,-
      csdr$idl:fmd_csq_st.pro,-
      csdr$idl:fmd_cvc_hif2.pro,-
      csdr$idl:fmd_cvc_hif3.pro,-
      csdr$idl:fmd_cvc_hif4.pro,-
      csdr$idl:fmd_cvc_high.pro,-
      csdr$idl:fmd_cvc_hres.pro,-
      csdr$idl:fmd_cvc_lowf.pro,-
      csdr$idl:fmd_cvc_st.pro,-
      csdr$idl:fmd_cvector.pro,-
      csdr$idl:fmd_cvectors_hi.pro,-
      csdr$idl:fmd_cvectors_hr.pro,-
      csdr$idl:fmd_cvectors_lo.pro,-
      csdr$idl:fmd_cvec_lhf_2.pro,-
      csdr$idl:fmd_cvec_lhf_3.pro,-
      csdr$idl:fmd_cvec_lhf_4.pro,-
      csdr$idl:fmd_cvec_lhs_2.pro,-
      csdr$idl:fmd_cvec_lhs_3.pro,-
      csdr$idl:fmd_cvec_lhs_4.pro,-
      csdr$idl:fmd_cvec_llf.pro,-
      csdr$idl:fmd_cvec_lls.pro,-
      csdr$idl:fmd_cvec_lsf.pro,-
      csdr$idl:fmd_cvec_rhf_2.pro,-
      csdr$idl:fmd_cvec_rhf_3.pro,-
      csdr$idl:fmd_cvec_rhf_4.pro,-
      csdr$idl:fmd_cvec_rhs_2.pro,-
      csdr$idl:fmd_cvec_rhs_3.pro,-
      csdr$idl:fmd_cvec_rhs_4.pro,-
      csdr$idl:fmd_cvec_rlf.pro,-
      csdr$idl:fmd_cvec_rls.pro,-
      csdr$idl:fmd_cvec_rsf.pro,-
      csdr$idl:fmd_czm_high.pro,-
      csdr$idl:fmd_czm_st.pro,-
      csdr$idl:fmd_deflt_quals.pro,-
      csdr$idl:fmd_destriper.pro,-
      csdr$idl:fmd_dgk_high.pro,-
      csdr$idl:fmd_dgk_hres.pro,-
      csdr$idl:fmd_dgk_lowf.pro,-
      csdr$idl:fmd_dgk_st.pro,-
      csdr$idl:fmd_dirbe_func.pro,-
      csdr$idl:fmd_dirbe_func0.pro,-
      csdr$idl:fmd_dirbe_funch.pro,-
      csdr$idl:fmd_dirbe_funcl.pro,-
      csdr$idl:fmd_dirbe_funcr.pro,-
      csdr$idl:fmd_dirbe_funcx.pro,-
      csdr$idl:fmd_dsp_high_2.pro,-
      csdr$idl:fmd_dsp_high_3.pro,-
      csdr$idl:fmd_dsp_high_4.pro,-
      csdr$idl:fmd_dsp_hres.pro,-
      csdr$idl:fmd_dsp_lhf_2.pro,-
      csdr$idl:fmd_dsp_lhf_2x.pro,-
      csdr$idl:fmd_dsp_lhf_3.pro,-
      csdr$idl:fmd_dsp_lhf_3x.pro,-
      csdr$idl:fmd_dsp_lhf_4.pro,-
      csdr$idl:fmd_dsp_lhf_4x.pro,-
      csdr$idl:fmd_dsp_lhs_2.pro,-
      csdr$idl:fmd_dsp_lhs_2x.pro,-
      csdr$idl:fmd_dsp_lhs_3.pro,-
      csdr$idl:fmd_dsp_lhs_3x.pro,-
      csdr$idl:fmd_dsp_lhs_4.pro,-
      csdr$idl:fmd_dsp_lhs_4x.pro,-
      csdr$idl:fmd_dsp_llf.pro,-
      csdr$idl:fmd_dsp_llfx.pro,-
      csdr$idl:fmd_dsp_lls.pro,-
      csdr$idl:fmd_dsp_llsx.pro,-
      csdr$idl:fmd_dsp_lowf.pro,-
      csdr$idl:fmd_dsp_lsf.pro,-
      csdr$idl:fmd_dsp_lsfx.pro,-
      csdr$idl:fmd_dsp_rhf_2.pro,-
      csdr$idl:fmd_dsp_rhf_2x.pro,-
      csdr$idl:fmd_dsp_rhf_3.pro,-
      csdr$idl:fmd_dsp_rhf_3x.pro,-
      csdr$idl:fmd_dsp_rhf_4.pro,-
      csdr$idl:fmd_dsp_rhf_4x.pro,-
      csdr$idl:fmd_dsp_rhs_2.pro,-
      csdr$idl:fmd_dsp_rhs_2x.pro,-
      csdr$idl:fmd_dsp_rhs_3.pro,-
      csdr$idl:fmd_dsp_rhs_3x.pro,-
      csdr$idl:fmd_dsp_rhs_4.pro,-
      csdr$idl:fmd_dsp_rhs_4x.pro,-
      csdr$idl:fmd_dsp_rlf.pro,-
      csdr$idl:fmd_dsp_rlfx.pro,-
      csdr$idl:fmd_dsp_rls.pro,-
      csdr$idl:fmd_dsp_rlsx.pro,-
      csdr$idl:fmd_dsp_rsf.pro,-
      csdr$idl:fmd_dsp_rsfx.pro,-
      csdr$idl:fmd_dvector.pro,-
      csdr$idl:fmd_dvec_lhf_2.pro,-
      csdr$idl:fmd_dvec_lhf_3.pro,-
      csdr$idl:fmd_dvec_lhf_4.pro,-
      csdr$idl:fmd_dvec_lhs_2.pro,-
      csdr$idl:fmd_dvec_lhs_3.pro,-
      csdr$idl:fmd_dvec_lhs_4.pro,-
      csdr$idl:fmd_dvec_llf.pro,-
      csdr$idl:fmd_dvec_lls.pro,-
      csdr$idl:fmd_dvec_lsf.pro,-
      csdr$idl:fmd_dvec_rhf_2.pro,-
      csdr$idl:fmd_dvec_rhf_3.pro,-
      csdr$idl:fmd_dvec_rhf_4.pro,-
      csdr$idl:fmd_dvec_rhs_2.pro,-
      csdr$idl:fmd_dvec_rhs_3.pro,-
      csdr$idl:fmd_dvec_rhs_4.pro,-
      csdr$idl:fmd_dvec_rlf.pro,-
      csdr$idl:fmd_dvec_rls.pro,-
      csdr$idl:fmd_dvec_rsf.pro,-
      csdr$idl:fmd_ejg_high_2.pro,-
      csdr$idl:fmd_ejg_high_3.pro,-
      csdr$idl:fmd_ejg_high_4.pro,-
      csdr$idl:fmd_ejg_hres.pro,-
      csdr$idl:fmd_ejg_lhf_2.pro,-
      csdr$idl:fmd_ejg_lhf_2x.pro,-
      csdr$idl:fmd_ejg_lhf_3.pro,-
      csdr$idl:fmd_ejg_lhf_3x.pro,-
      csdr$idl:fmd_ejg_lhf_4.pro,-
      csdr$idl:fmd_ejg_lhf_4x.pro,-
      csdr$idl:fmd_ejg_lhs_2.pro,-
      csdr$idl:fmd_ejg_lhs_2x.pro,-
      csdr$idl:fmd_ejg_lhs_3.pro,-
      csdr$idl:fmd_ejg_lhs_3x.pro,-
      csdr$idl:fmd_ejg_lhs_4.pro,-
      csdr$idl:fmd_ejg_lhs_4x.pro,-
      csdr$idl:fmd_ejg_llf.pro,-
      csdr$idl:fmd_ejg_llfx.pro,-
      csdr$idl:fmd_ejg_lls.pro,-
      csdr$idl:fmd_ejg_llsx.pro,-
      csdr$idl:fmd_ejg_lowf.pro,-
      csdr$idl:fmd_ejg_lsf.pro,-
      csdr$idl:fmd_ejg_lsfx.pro,-
      csdr$idl:fmd_ejg_rhf_2.pro,-
      csdr$idl:fmd_ejg_rhf_2x.pro,-
      csdr$idl:fmd_ejg_rhf_3.pro,-
      csdr$idl:fmd_ejg_rhf_3x.pro,-
      csdr$idl:fmd_ejg_rhf_4.pro,-
      csdr$idl:fmd_ejg_rhf_4x.pro,-
      csdr$idl:fmd_ejg_rhs_2.pro,-
      csdr$idl:fmd_ejg_rhs_2x.pro,-
      csdr$idl:fmd_ejg_rhs_3.pro,-
      csdr$idl:fmd_ejg_rhs_3x.pro,-
      csdr$idl:fmd_ejg_rhs_4.pro,-
      csdr$idl:fmd_ejg_rhs_4x.pro,-
      csdr$idl:fmd_ejg_rlf.pro,-
      csdr$idl:fmd_ejg_rlfx.pro,-
      csdr$idl:fmd_ejg_rls.pro,-
      csdr$idl:fmd_ejg_rlsx.pro,-
      csdr$idl:fmd_ejg_rsf.pro,-
      csdr$idl:fmd_ejg_rsfx.pro,-
      csdr$idl:fmd_errmat.pro,-
      csdr$idl:fmd_func.pro,-
      csdr$idl:fmd_high_err.pro,-
      csdr$idl:fmd_high_skymap.pro,-
      csdr$idl:fmd_leg_func.pro,-
      csdr$idl:fmd_match.pro,-
      csdr$idl:fmd_model_wgt.pro,-
      csdr$idl:fmd_ost_high.pro,-
      csdr$idl:fmd_ost_hif2.pro,-
      csdr$idl:fmd_ost_hif3.pro,-
      csdr$idl:fmd_ost_hif4.pro,-
      csdr$idl:fmd_ost_hres.pro,-
      csdr$idl:fmd_ost_lowf.pro,-
      csdr$idl:fmd_ost_st.pro,-
      csdr$idl:fmd_pixel_wgt.pro,-
      csdr$idl:fmd_pix_sum.pro,-
      csdr$idl:fmd_pst_hif2.pro,-
      csdr$idl:fmd_pst_hif3.pro,-
      csdr$idl:fmd_pst_hif4.pro,-
      csdr$idl:fmd_pst_hres.pro,-
      csdr$idl:fmd_pst_lowf.pro,-
      csdr$idl:fmd_pst_st.pro,-
      csdr$idl:fmd_pzm_high.pro,-
      csdr$idl:fmd_pzm_st.pro,-
      csdr$idl:fmd_quals.pro,-
      csdr$idl:fmd_read.pro,-
      csdr$idl:fmd_readv.pro,-
      csdr$idl:fmd_read_data.pro,-
      csdr$idl:fmd_resid.pro,-
      csdr$idl:fmd_resid_high2.pro,-
      csdr$idl:fmd_resid_high3.pro,-
      csdr$idl:fmd_resid_high4.pro,-
      csdr$idl:fmd_resid_hres.pro,-
      csdr$idl:fmd_resid_lowf.pro,-
      csdr$idl:fmd_save_quals.pro,-
      csdr$idl:fmd_sky_hif2.pro,-
      csdr$idl:fmd_sky_hif3.pro,-
      csdr$idl:fmd_sky_hif4.pro,-
      csdr$idl:fmd_sky_high.pro,-
      csdr$idl:fmd_sky_hres.pro,-
      csdr$idl:fmd_sky_lowf.pro,-
      csdr$idl:fmd_sky_st.pro,-
      csdr$idl:fmd_variance.pro,-
      csdr$idl:fmd_variance2.pro,-
      csdr$idl:fmd_wgts_high.pro,-
      csdr$idl:fmd_wgts_hres.pro,-
      csdr$idl:fmd_wgts_lhf.pro,-
      csdr$idl:fmd_wgts_lhfx.pro,-
      csdr$idl:fmd_wgts_lhs.pro,-
      csdr$idl:fmd_wgts_lhsx.pro,-
      csdr$idl:fmd_wgts_llf.pro,-
      csdr$idl:fmd_wgts_llfx.pro,-
      csdr$idl:fmd_wgts_lls.pro,-
      csdr$idl:fmd_wgts_llsx.pro,-
      csdr$idl:fmd_wgts_lowf.pro,-
      csdr$idl:fmd_wgts_lsf.pro,-
      csdr$idl:fmd_wgts_lsfx.pro,-
      csdr$idl:fmd_wgts_rhf.pro,-
      csdr$idl:fmd_wgts_rhfx.pro,-
      csdr$idl:fmd_wgts_rhs.pro,-
      csdr$idl:fmd_wgts_rhsx.pro,-
      csdr$idl:fmd_wgts_rlf.pro,-
      csdr$idl:fmd_wgts_rlfx.pro,-
      csdr$idl:fmd_wgts_rls.pro,-
      csdr$idl:fmd_wgts_rlsx.pro,-
      csdr$idl:fmd_wgts_rsf.pro,-
      csdr$idl:fmd_wgts_rsfx.pro,-
      csdr$idl:fmd_zodi_high.pro,-
      csdr$idl:fmd_zodi_hi_2.pro,-
      csdr$idl:fmd_zodi_hi_3.pro,-
      csdr$idl:fmd_zodi_hi_4.pro,-
      csdr$idl:fmd_zodi_skymap.pro,-
      csdr$idl:fmd_zsub_high_2.pro,-
      csdr$idl:fmd_zsub_high_3.pro,-
      csdr$idl:fmd_zsub_high_4.pro
 ! Facility fmd has been installed.
 
csdr$idl:fmd_adt_sec.pro      : fmd_adt_sec.pro      
csdr$idl:fmd_apply.pro        : fmd_apply.pro        
csdr$idl:fmd_avc_hif2.pro     : fmd_avc_hif2.pro
csdr$idl:fmd_avc_hif3.pro     : fmd_avc_hif3.pro
csdr$idl:fmd_avc_hif4.pro     : fmd_avc_hif4.pro
csdr$idl:fmd_avc_high.pro     : fmd_avc_high.pro
csdr$idl:fmd_avc_hres.pro     : fmd_avc_hres.pro
csdr$idl:fmd_avc_lowf.pro     : fmd_avc_lowf.pro
csdr$idl:fmd_avc_st.pro       : fmd_avc_st.pro
csdr$idl:fmd_avec_high.pro    : fmd_avec_high.pro    
csdr$idl:fmd_avec_high_2.pro  : fmd_avec_high_2.pro  
csdr$idl:fmd_avec_high_3.pro  : fmd_avec_high_3.pro  
csdr$idl:fmd_avec_high_4.pro  : fmd_avec_high_4.pro  
csdr$idl:fmd_avec_hres.pro    : fmd_avec_hres.pro    
csdr$idl:fmd_avec_lowf.pro    : fmd_avec_lowf.pro    
csdr$idl:fmd_bad2_hi.pro      : fmd_bad2_hi.pro      
csdr$idl:fmd_bad2_hr.pro      : fmd_bad2_hr.pro      
csdr$idl:fmd_bad2_lhf.pro     : fmd_bad2_lhf.pro     
csdr$idl:fmd_bad2_lhs.pro     : fmd_bad2_lhs.pro     
csdr$idl:fmd_bad2_llf.pro     : fmd_bad2_llf.pro
csdr$idl:fmd_bad2_lls.pro     : fmd_bad2_lls.pro     
csdr$idl:fmd_bad2_lo.pro      : fmd_bad2_lo.pro      
csdr$idl:fmd_bad2_lsf.pro     : fmd_bad2_lsf.pro     
csdr$idl:fmd_bad2_rhf.pro     : fmd_bad2_rhf.pro     
csdr$idl:fmd_bad2_rhs.pro     : fmd_bad2_rhs.pro     
csdr$idl:fmd_bad2_rlf.pro     : fmd_bad2_rlf.pro     
csdr$idl:fmd_bad2_rls.pro     : fmd_bad2_rls.pro     
csdr$idl:fmd_bad2_rsf.pro     : fmd_bad2_rsf.pro     
csdr$idl:fmd_bad_hi.pro       : fmd_bad_hi.pro       
csdr$idl:fmd_bad_hr.pro       : fmd_bad_hr.pro       
csdr$idl:fmd_bad_lhf.pro      : fmd_bad_lhf.pro      
csdr$idl:fmd_bad_lhs.pro      : fmd_bad_lhs.pro      
csdr$idl:fmd_bad_llf.pro      : fmd_bad_llf.pro      
csdr$idl:fmd_bad_lls.pro      : fmd_bad_lls.pro      
csdr$idl:fmd_bad_lo.pro       : fmd_bad_lo.pro       
csdr$idl:fmd_bad_lsf.pro      : fmd_bad_lsf.pro      
csdr$idl:fmd_bad_rhf.pro      : fmd_bad_rhf.pro      
csdr$idl:fmd_bad_rhs.pro      : fmd_bad_rhs.pro      
csdr$idl:fmd_bad_rlf.pro      : fmd_bad_rlf.pro      
csdr$idl:fmd_bad_rls.pro      : fmd_bad_rls.pro      
csdr$idl:fmd_bad_rsf.pro      : fmd_bad_rsf.pro      
csdr$idl:fmd_bandspec.pro     : fmd_bandspec.pro     
csdr$idl:fmd_cal_chi2.pro     : fmd_cal_chi2.pro     
csdr$idl:fmd_chi2.pro         : fmd_chi2.pro         
csdr$idl:fmd_chi2_high_2.pro  : fmd_chi2_high_2.pro  
csdr$idl:fmd_chi2_high_3.pro  : fmd_chi2_high_3.pro  
csdr$idl:fmd_chi2_high_4.pro  : fmd_chi2_high_4.pro  
csdr$idl:fmd_chi2_hres.pro    : fmd_chi2_hres.pro    
csdr$idl:fmd_chi2_lhf_2.pro   : fmd_chi2_lhf_2.pro   
csdr$idl:fmd_chi2_lhf_2x.pro  : fmd_chi2_lhf_2x.pro  
csdr$idl:fmd_chi2_lhf_3.pro   : fmd_chi2_lhf_3.pro   
csdr$idl:fmd_chi2_lhf_3x.pro  : fmd_chi2_lhf_3x.pro  
csdr$idl:fmd_chi2_lhf_4.pro   : fmd_chi2_lhf_4.pro   
csdr$idl:fmd_chi2_lhf_4x.pro  : fmd_chi2_lhf_4x.pro  
csdr$idl:fmd_chi2_lhs_2.pro   : fmd_chi2_lhs_2.pro   
csdr$idl:fmd_chi2_lhs_2x.pro  : fmd_chi2_lhs_2x.pro  
csdr$idl:fmd_chi2_lhs_3.pro   : fmd_chi2_lhs_3.pro   
csdr$idl:fmd_chi2_lhs_3x.pro  : fmd_chi2_lhs_3x.pro  
csdr$idl:fmd_chi2_lhs_4.pro   : fmd_chi2_lhs_4.pro   
csdr$idl:fmd_chi2_lhs_4x.pro  : fmd_chi2_lhs_4x.pro  
csdr$idl:fmd_chi2_llf.pro     : fmd_chi2_llf.pro     
csdr$idl:fmd_chi2_llfx.pro    : fmd_chi2_llfx.pro    
csdr$idl:fmd_chi2_lls.pro     : fmd_chi2_lls.pro     
csdr$idl:fmd_chi2_llsx.pro    : fmd_chi2_llsx.pro    
csdr$idl:fmd_chi2_lowf.pro    : fmd_chi2_lowf.pro    
csdr$idl:fmd_chi2_lsf.pro     : fmd_chi2_lsf.pro     
csdr$idl:fmd_chi2_lsfx.pro    : fmd_chi2_lsfx.pro    
csdr$idl:fmd_chi2_rhf_2.pro   : fmd_chi2_rhf_2.pro   
csdr$idl:fmd_chi2_rhf_2x.pro  : fmd_chi2_rhf_2x.pro  
csdr$idl:fmd_chi2_rhf_3.pro   : fmd_chi2_rhf_3.pro   
csdr$idl:fmd_chi2_rhf_3x.pro  : fmd_chi2_rhf_3x.pro  
csdr$idl:fmd_chi2_rhf_4.pro   : fmd_chi2_rhf_4.pro   
csdr$idl:fmd_chi2_rhf_4x.pro  : fmd_chi2_rhf_4x.pro  
csdr$idl:fmd_chi2_rhs_2.pro   : fmd_chi2_rhs_2.pro   
csdr$idl:fmd_chi2_rhs_2x.pro  : fmd_chi2_rhs_2x.pro  
csdr$idl:fmd_chi2_rhs_3.pro   : fmd_chi2_rhs_3.pro   
csdr$idl:fmd_chi2_rhs_3x.pro  : fmd_chi2_rhs_3x.pro  
csdr$idl:fmd_chi2_rhs_4.pro   : fmd_chi2_rhs_4.pro   
csdr$idl:fmd_chi2_rhs_4x.pro  : fmd_chi2_rhs_4x.pro  
csdr$idl:fmd_chi2_rlf.pro     : fmd_chi2_rlf.pro     
csdr$idl:fmd_chi2_rlfx.pro    : fmd_chi2_rlfx.pro    
csdr$idl:fmd_chi2_rls.pro     : fmd_chi2_rls.pro     
csdr$idl:fmd_chi2_rlsx.pro    : fmd_chi2_rlsx.pro    
csdr$idl:fmd_chi2_rsf.pro     : fmd_chi2_rsf.pro     
csdr$idl:fmd_chi2_rsfx.pro    : fmd_chi2_rsfx.pro    
csdr$idl:fmd_choleskey.pro    : fmd_choleskey.pro    
csdr$idl:fmd_concat2_hi.pro   : fmd_concat2_hi.pro
csdr$idl:fmd_concat2_hr.pro   : fmd_concat2_hr.pro
csdr$idl:fmd_concat_high.pro  : fmd_concat_high.pro  
csdr$idl:fmd_concat_hres.pro  : fmd_concat_hres.pro  
csdr$idl:fmd_concat_lowf.pro  : fmd_concat_lowf.pro  
csdr$idl:fmd_covar.pro        : fmd_covar.pro        
csdr$idl:fmd_covar_high.pro   : fmd_covar_high.pro   
csdr$idl:fmd_covar_high2.pro  : fmd_covar_high2.pro  
csdr$idl:fmd_covar_high3.pro  : fmd_covar_high3.pro  
csdr$idl:fmd_covar_high4.pro  : fmd_covar_high4.pro  
csdr$idl:fmd_covar_hi_23.pro  : fmd_covar_hi_23.pro  
csdr$idl:fmd_covar_hi_24.pro  : fmd_covar_hi_24.pro  
csdr$idl:fmd_covar_hi_34.pro  : fmd_covar_hi_34.pro  
csdr$idl:fmd_covar_hres.pro   : fmd_covar_hres.pro   
csdr$idl:fmd_covar_hres1.pro  : fmd_covar_hres1.pro
csdr$idl:fmd_covar_hres2.pro  : fmd_covar_hres2.pro
csdr$idl:fmd_covar_hres3.pro  : fmd_covar_hres3.pro
csdr$idl:fmd_covar_hr_12.pro  : fmd_covar_hr_12.pro
csdr$idl:fmd_covar_hr_13.pro  : fmd_covar_hr_13.pro
csdr$idl:fmd_covar_hr_23.pro  : fmd_covar_hr_23.pro
csdr$idl:fmd_covar_lowf.pro   : fmd_covar_lowf.pro   
csdr$idl:fmd_covar_off.pro    : fmd_covar_off.pro    
csdr$idl:fmd_cov_hif2.pro     : fmd_cov_hif2.pro
csdr$idl:fmd_cov_hif3.pro     : fmd_cov_hif3.pro
csdr$idl:fmd_cov_hif4.pro     : fmd_cov_hif4.pro
csdr$idl:fmd_cov_high.pro     : fmd_cov_high.pro
csdr$idl:fmd_cov_hres.pro     : fmd_cov_hres.pro
csdr$idl:fmd_cov_lowf.pro     : fmd_cov_lowf.pro
csdr$idl:fmd_cov_st.pro       : fmd_cov_st.pro
csdr$idl:fmd_csq_hif2.pro     : fmd_csq_hif2.pro
csdr$idl:fmd_csq_hif3.pro     : fmd_csq_hif3.pro
csdr$idl:fmd_csq_hif4.pro     : fmd_csq_hif4.pro
csdr$idl:fmd_csq_high.pro     : fmd_csq_high.pro
csdr$idl:fmd_csq_hres.pro     : fmd_csq_hres.pro
csdr$idl:fmd_csq_lowf.pro     : fmd_csq_lowf.pro
csdr$idl:fmd_csq_st.pro       : fmd_csq_st.pro
csdr$idl:fmd_cvc_hif2.pro     : fmd_cvc_hif2.pro
csdr$idl:fmd_cvc_hif3.pro     : fmd_cvc_hif3.pro
csdr$idl:fmd_cvc_hif4.pro     : fmd_cvc_hif4.pro
csdr$idl:fmd_cvc_high.pro     : fmd_cvc_high.pro
csdr$idl:fmd_cvc_hres.pro     : fmd_cvc_hres.pro
csdr$idl:fmd_cvc_lowf.pro     : fmd_cvc_lowf.pro
csdr$idl:fmd_cvc_st.pro       : fmd_cvc_st.pro
csdr$idl:fmd_cvector.pro      : fmd_cvector.pro      
csdr$idl:fmd_cvectors_hi.pro  : fmd_cvectors_hi.pro  
csdr$idl:fmd_cvectors_hr.pro  : fmd_cvectors_hr.pro  
csdr$idl:fmd_cvectors_lo.pro  : fmd_cvectors_lo.pro  
csdr$idl:fmd_cvec_lhf_2.pro   : fmd_cvec_lhf_2.pro   
csdr$idl:fmd_cvec_lhf_3.pro   : fmd_cvec_lhf_3.pro   
csdr$idl:fmd_cvec_lhf_4.pro   : fmd_cvec_lhf_4.pro   
csdr$idl:fmd_cvec_lhs_2.pro   : fmd_cvec_lhs_2.pro   
csdr$idl:fmd_cvec_lhs_3.pro   : fmd_cvec_lhs_3.pro   
csdr$idl:fmd_cvec_lhs_4.pro   : fmd_cvec_lhs_4.pro   
csdr$idl:fmd_cvec_llf.pro     : fmd_cvec_llf.pro     
csdr$idl:fmd_cvec_lls.pro     : fmd_cvec_lls.pro     
csdr$idl:fmd_cvec_lsf.pro     : fmd_cvec_lsf.pro     
csdr$idl:fmd_cvec_rhf_2.pro   : fmd_cvec_rhf_2.pro   
csdr$idl:fmd_cvec_rhf_3.pro   : fmd_cvec_rhf_3.pro   
csdr$idl:fmd_cvec_rhf_4.pro   : fmd_cvec_rhf_4.pro   
csdr$idl:fmd_cvec_rhs_2.pro   : fmd_cvec_rhs_2.pro   
csdr$idl:fmd_cvec_rhs_3.pro   : fmd_cvec_rhs_3.pro   
csdr$idl:fmd_cvec_rhs_4.pro   : fmd_cvec_rhs_4.pro   
csdr$idl:fmd_cvec_rlf.pro     : fmd_cvec_rlf.pro     
csdr$idl:fmd_cvec_rls.pro     : fmd_cvec_rls.pro     
csdr$idl:fmd_cvec_rsf.pro     : fmd_cvec_rsf.pro     
csdr$idl:fmd_czm_high.pro     : fmd_czm_high.pro
csdr$idl:fmd_czm_st.pro       : fmd_czm_st.pro
csdr$idl:fmd_deflt_quals.pro  : fmd_deflt_quals.pro
csdr$idl:fmd_destriper.pro    : fmd_destriper.pro    
csdr$idl:fmd_dgk_high.pro     : fmd_dgk_high.pro
csdr$idl:fmd_dgk_hres.pro     : fmd_dgk_hres.pro
csdr$idl:fmd_dgk_lowf.pro     : fmd_dgk_lowf.pro
csdr$idl:fmd_dgk_st.pro       : fmd_dgk_st.pro
csdr$idl:fmd_dirbe_func.pro   : fmd_dirbe_func.pro   
csdr$idl:fmd_dirbe_func0.pro  : fmd_dirbe_func0.pro   
csdr$idl:fmd_dirbe_funch.pro  : fmd_dirbe_funch.pro  
csdr$idl:fmd_dirbe_funcl.pro  : fmd_dirbe_funcl.pro  
csdr$idl:fmd_dirbe_funcr.pro  : fmd_dirbe_funcr.pro  
csdr$idl:fmd_dirbe_funcx.pro  : fmd_dirbe_funcx.pro
csdr$idl:fmd_dsp_high_2.pro   : fmd_dsp_high_2.pro   
csdr$idl:fmd_dsp_high_3.pro   : fmd_dsp_high_3.pro   
csdr$idl:fmd_dsp_high_4.pro   : fmd_dsp_high_4.pro   
csdr$idl:fmd_dsp_hres.pro     : fmd_dsp_hres.pro     
csdr$idl:fmd_dsp_lhf_2.pro    : fmd_dsp_lhf_2.pro    
csdr$idl:fmd_dsp_lhf_2x.pro   : fmd_dsp_lhf_2x.pro   
csdr$idl:fmd_dsp_lhf_3.pro    : fmd_dsp_lhf_3.pro    
csdr$idl:fmd_dsp_lhf_3x.pro   : fmd_dsp_lhf_3x.pro   
csdr$idl:fmd_dsp_lhf_4.pro    : fmd_dsp_lhf_4.pro    
csdr$idl:fmd_dsp_lhf_4x.pro   : fmd_dsp_lhf_4x.pro   
csdr$idl:fmd_dsp_lhs_2.pro    : fmd_dsp_lhs_2.pro    
csdr$idl:fmd_dsp_lhs_2x.pro   : fmd_dsp_lhs_2x.pro   
csdr$idl:fmd_dsp_lhs_3.pro    : fmd_dsp_lhs_3.pro    
csdr$idl:fmd_dsp_lhs_3x.pro   : fmd_dsp_lhs_3x.pro   
csdr$idl:fmd_dsp_lhs_4.pro    : fmd_dsp_lhs_4.pro    
csdr$idl:fmd_dsp_lhs_4x.pro   : fmd_dsp_lhs_4x.pro   
csdr$idl:fmd_dsp_llf.pro      : fmd_dsp_llf.pro      
csdr$idl:fmd_dsp_llfx.pro     : fmd_dsp_llfx.pro     
csdr$idl:fmd_dsp_lls.pro      : fmd_dsp_lls.pro      
csdr$idl:fmd_dsp_llsx.pro     : fmd_dsp_llsx.pro     
csdr$idl:fmd_dsp_lowf.pro     : fmd_dsp_lowf.pro     
csdr$idl:fmd_dsp_lsf.pro      : fmd_dsp_lsf.pro      
csdr$idl:fmd_dsp_lsfx.pro     : fmd_dsp_lsfx.pro     
csdr$idl:fmd_dsp_rhf_2.pro    : fmd_dsp_rhf_2.pro    
csdr$idl:fmd_dsp_rhf_2x.pro   : fmd_dsp_rhf_2x.pro   
csdr$idl:fmd_dsp_rhf_3.pro    : fmd_dsp_rhf_3.pro    
csdr$idl:fmd_dsp_rhf_3x.pro   : fmd_dsp_rhf_3x.pro   
csdr$idl:fmd_dsp_rhf_4.pro    : fmd_dsp_rhf_4.pro    
csdr$idl:fmd_dsp_rhf_4x.pro   : fmd_dsp_rhf_4x.pro   
csdr$idl:fmd_dsp_rhs_2.pro    : fmd_dsp_rhs_2.pro    
csdr$idl:fmd_dsp_rhs_2x.pro   : fmd_dsp_rhs_2x.pro   
csdr$idl:fmd_dsp_rhs_3.pro    : fmd_dsp_rhs_3.pro    
csdr$idl:fmd_dsp_rhs_3x.pro   : fmd_dsp_rhs_3x.pro   
csdr$idl:fmd_dsp_rhs_4.pro    : fmd_dsp_rhs_4.pro    
csdr$idl:fmd_dsp_rhs_4x.pro   : fmd_dsp_rhs_4x.pro   
csdr$idl:fmd_dsp_rlf.pro      : fmd_dsp_rlf.pro      
csdr$idl:fmd_dsp_rlfx.pro     : fmd_dsp_rlfx.pro     
csdr$idl:fmd_dsp_rls.pro      : fmd_dsp_rls.pro      
csdr$idl:fmd_dsp_rlsx.pro     : fmd_dsp_rlsx.pro     
csdr$idl:fmd_dsp_rsf.pro      : fmd_dsp_rsf.pro      
csdr$idl:fmd_dsp_rsfx.pro     : fmd_dsp_rsfx.pro     
csdr$idl:fmd_dvector.pro      : fmd_dvector.pro      
csdr$idl:fmd_dvec_lhf_2.pro   : fmd_dvec_lhf_2.pro   
csdr$idl:fmd_dvec_lhf_3.pro   : fmd_dvec_lhf_3.pro   
csdr$idl:fmd_dvec_lhf_4.pro   : fmd_dvec_lhf_4.pro   
csdr$idl:fmd_dvec_lhs_2.pro   : fmd_dvec_lhs_2.pro   
csdr$idl:fmd_dvec_lhs_3.pro   : fmd_dvec_lhs_3.pro   
csdr$idl:fmd_dvec_lhs_4.pro   : fmd_dvec_lhs_4.pro   
csdr$idl:fmd_dvec_llf.pro     : fmd_dvec_llf.pro     
csdr$idl:fmd_dvec_lls.pro     : fmd_dvec_lls.pro     
csdr$idl:fmd_dvec_lsf.pro     : fmd_dvec_lsf.pro     
csdr$idl:fmd_dvec_rhf_2.pro   : fmd_dvec_rhf_2.pro   
csdr$idl:fmd_dvec_rhf_3.pro   : fmd_dvec_rhf_3.pro   
csdr$idl:fmd_dvec_rhf_4.pro   : fmd_dvec_rhf_4.pro   
csdr$idl:fmd_dvec_rhs_2.pro   : fmd_dvec_rhs_2.pro   
csdr$idl:fmd_dvec_rhs_3.pro   : fmd_dvec_rhs_3.pro   
csdr$idl:fmd_dvec_rhs_4.pro   : fmd_dvec_rhs_4.pro   
csdr$idl:fmd_dvec_rlf.pro     : fmd_dvec_rlf.pro     
csdr$idl:fmd_dvec_rls.pro     : fmd_dvec_rls.pro     
csdr$idl:fmd_dvec_rsf.pro     : fmd_dvec_rsf.pro     
csdr$idl:fmd_ejg_high_2.pro   : fmd_ejg_high_2.pro   
csdr$idl:fmd_ejg_high_3.pro   : fmd_ejg_high_3.pro   
csdr$idl:fmd_ejg_high_4.pro   : fmd_ejg_high_4.pro   
csdr$idl:fmd_ejg_hres.pro     : fmd_ejg_hres.pro     
csdr$idl:fmd_ejg_lhf_2.pro    : fmd_ejg_lhf_2.pro    
csdr$idl:fmd_ejg_lhf_2x.pro   : fmd_ejg_lhf_2x.pro   
csdr$idl:fmd_ejg_lhf_3.pro    : fmd_ejg_lhf_3.pro    
csdr$idl:fmd_ejg_lhf_3x.pro   : fmd_ejg_lhf_3x.pro   
csdr$idl:fmd_ejg_lhf_4.pro    : fmd_ejg_lhf_4.pro    
csdr$idl:fmd_ejg_lhf_4x.pro   : fmd_ejg_lhf_4x.pro   
csdr$idl:fmd_ejg_lhs_2.pro    : fmd_ejg_lhs_2.pro    
csdr$idl:fmd_ejg_lhs_2x.pro   : fmd_ejg_lhs_2x.pro   
csdr$idl:fmd_ejg_lhs_3.pro    : fmd_ejg_lhs_3.pro    
csdr$idl:fmd_ejg_lhs_3x.pro   : fmd_ejg_lhs_3x.pro   
csdr$idl:fmd_ejg_lhs_4.pro    : fmd_ejg_lhs_4.pro    
csdr$idl:fmd_ejg_lhs_4x.pro   : fmd_ejg_lhs_4x.pro   
csdr$idl:fmd_ejg_llf.pro      : fmd_ejg_llf.pro      
csdr$idl:fmd_ejg_llfx.pro     : fmd_ejg_llfx.pro     
csdr$idl:fmd_ejg_lls.pro      : fmd_ejg_lls.pro      
csdr$idl:fmd_ejg_llsx.pro     : fmd_ejg_llsx.pro     
csdr$idl:fmd_ejg_lowf.pro     : fmd_ejg_lowf.pro     
csdr$idl:fmd_ejg_lsf.pro      : fmd_ejg_lsf.pro      
csdr$idl:fmd_ejg_lsfx.pro     : fmd_ejg_lsfx.pro     
csdr$idl:fmd_ejg_rhf_2.pro    : fmd_ejg_rhf_2.pro    
csdr$idl:fmd_ejg_rhf_2x.pro   : fmd_ejg_rhf_2x.pro   
csdr$idl:fmd_ejg_rhf_3.pro    : fmd_ejg_rhf_3.pro    
csdr$idl:fmd_ejg_rhf_3x.pro   : fmd_ejg_rhf_3x.pro   
csdr$idl:fmd_ejg_rhf_4.pro    : fmd_ejg_rhf_4.pro    
csdr$idl:fmd_ejg_rhf_4x.pro   : fmd_ejg_rhf_4x.pro   
csdr$idl:fmd_ejg_rhs_2.pro    : fmd_ejg_rhs_2.pro    
csdr$idl:fmd_ejg_rhs_2x.pro   : fmd_ejg_rhs_2x.pro   
csdr$idl:fmd_ejg_rhs_3.pro    : fmd_ejg_rhs_3.pro    
csdr$idl:fmd_ejg_rhs_3x.pro   : fmd_ejg_rhs_3x.pro   
csdr$idl:fmd_ejg_rhs_4.pro    : fmd_ejg_rhs_4.pro    
csdr$idl:fmd_ejg_rhs_4x.pro   : fmd_ejg_rhs_4x.pro   
csdr$idl:fmd_ejg_rlf.pro      : fmd_ejg_rlf.pro      
csdr$idl:fmd_ejg_rlfx.pro     : fmd_ejg_rlfx.pro     
csdr$idl:fmd_ejg_rls.pro      : fmd_ejg_rls.pro      
csdr$idl:fmd_ejg_rlsx.pro     : fmd_ejg_rlsx.pro     
csdr$idl:fmd_ejg_rsf.pro      : fmd_ejg_rsf.pro      
csdr$idl:fmd_ejg_rsfx.pro     : fmd_ejg_rsfx.pro     
csdr$idl:fmd_errmat.pro       : fmd_errmat.pro       
csdr$idl:fmd_func.pro         : fmd_func.pro         
csdr$idl:fmd_high_err.pro     : fmd_high_err.pro
csdr$idl:fmd_high_skymap.pro  : fmd_high_skymap.pro
csdr$idl:fmd_leg_func.pro     : fmd_leg_func.pro     
csdr$idl:fmd_match.pro        : fmd_match.pro
csdr$idl:fmd_model_wgt.pro    : fmd_model_wgt.pro    
csdr$idl:fmd_ost_high.pro     : fmd_ost_high.pro
csdr$idl:fmd_ost_hif2.pro     : fmd_ost_hif2.pro
csdr$idl:fmd_ost_hif3.pro     : fmd_ost_hif3.pro
csdr$idl:fmd_ost_hif4.pro     : fmd_ost_hif4.pro
csdr$idl:fmd_ost_hres.pro     : fmd_ost_hres.pro
csdr$idl:fmd_ost_lowf.pro     : fmd_ost_lowf.pro
csdr$idl:fmd_ost_st.pro       : fmd_ost_st.pro
csdr$idl:fmd_pixel_wgt.pro    : fmd_pixel_wgt.pro    
csdr$idl:fmd_pix_sum.pro      : fmd_pix_sum.pro      
csdr$idl:fmd_pst_hif2.pro     : fmd_pst_hif2.pro
csdr$idl:fmd_pst_hif3.pro     : fmd_pst_hif3.pro
csdr$idl:fmd_pst_hif4.pro     : fmd_pst_hif4.pro
csdr$idl:fmd_pst_hres.pro     : fmd_pst_hres.pro
csdr$idl:fmd_pst_lowf.pro     : fmd_pst_lowf.pro
csdr$idl:fmd_pst_st.pro       : fmd_pst_st.pro
csdr$idl:fmd_pzm_high.pro     : fmd_pzm_high.pro
csdr$idl:fmd_pzm_st.pro       : fmd_pzm_st.pro
csdr$idl:fmd_quals.pro        : fmd_quals.pro        
csdr$idl:fmd_read.pro         : fmd_read.pro
csdr$idl:fmd_readv.pro        : fmd_readv.pro
csdr$idl:fmd_read_data.pro    : fmd_read_data.pro    
csdr$idl:fmd_resid.pro        : fmd_resid.pro        
csdr$idl:fmd_resid_high2.pro  : fmd_resid_high2.pro  
csdr$idl:fmd_resid_high3.pro  : fmd_resid_high3.pro  
csdr$idl:fmd_resid_high4.pro  : fmd_resid_high4.pro  
csdr$idl:fmd_resid_hres.pro   : fmd_resid_hres.pro   
csdr$idl:fmd_resid_lowf.pro   : fmd_resid_lowf.pro   
csdr$idl:fmd_save_quals.pro   : fmd_save_quals.pro   
csdr$idl:fmd_sky_hif2.pro     : fmd_sky_hif2.pro
csdr$idl:fmd_sky_hif3.pro     : fmd_sky_hif3.pro
csdr$idl:fmd_sky_hif4.pro     : fmd_sky_hif4.pro
csdr$idl:fmd_sky_high.pro     : fmd_sky_high.pro
csdr$idl:fmd_sky_hres.pro     : fmd_sky_hres.pro
csdr$idl:fmd_sky_lowf.pro     : fmd_sky_lowf.pro
csdr$idl:fmd_sky_st.pro       : fmd_sky_st.pro
csdr$idl:fmd_variance.pro     : fmd_variance.pro     
csdr$idl:fmd_variance2.pro    : fmd_variance2.pro    
csdr$idl:fmd_wgts_high.pro    : fmd_wgts_high.pro    
csdr$idl:fmd_wgts_hres.pro    : fmd_wgts_hres.pro    
csdr$idl:fmd_wgts_lhf.pro     : fmd_wgts_lhf.pro     
csdr$idl:fmd_wgts_lhfx.pro    : fmd_wgts_lhfx.pro    
csdr$idl:fmd_wgts_lhs.pro     : fmd_wgts_lhs.pro     
csdr$idl:fmd_wgts_lhsx.pro    : fmd_wgts_lhsx.pro    
csdr$idl:fmd_wgts_llf.pro     : fmd_wgts_llf.pro     
csdr$idl:fmd_wgts_llfx.pro    : fmd_wgts_llfx.pro    
csdr$idl:fmd_wgts_lls.pro     : fmd_wgts_lls.pro     
csdr$idl:fmd_wgts_llsx.pro    : fmd_wgts_llsx.pro    
csdr$idl:fmd_wgts_lowf.pro    : fmd_wgts_lowf.pro    
csdr$idl:fmd_wgts_lsf.pro     : fmd_wgts_lsf.pro     
csdr$idl:fmd_wgts_lsfx.pro    : fmd_wgts_lsfx.pro    
csdr$idl:fmd_wgts_rhf.pro     : fmd_wgts_rhf.pro     
csdr$idl:fmd_wgts_rhfx.pro    : fmd_wgts_rhfx.pro    
csdr$idl:fmd_wgts_rhs.pro     : fmd_wgts_rhs.pro     
csdr$idl:fmd_wgts_rhsx.pro    : fmd_wgts_rhsx.pro    
csdr$idl:fmd_wgts_rlf.pro     : fmd_wgts_rlf.pro     
csdr$idl:fmd_wgts_rlfx.pro    : fmd_wgts_rlfx.pro    
csdr$idl:fmd_wgts_rls.pro     : fmd_wgts_rls.pro     
csdr$idl:fmd_wgts_rlsx.pro    : fmd_wgts_rlsx.pro    
csdr$idl:fmd_wgts_rsf.pro     : fmd_wgts_rsf.pro     
csdr$idl:fmd_wgts_rsfx.pro    : fmd_wgts_rsfx.pro    
csdr$idl:fmd_zodi_high.pro    : fmd_zodi_high.pro    
csdr$idl:fmd_zodi_hi_2.pro    : fmd_zodi_hi_2.pro    
csdr$idl:fmd_zodi_hi_3.pro    : fmd_zodi_hi_3.pro    
csdr$idl:fmd_zodi_hi_4.pro    : fmd_zodi_hi_4.pro    
csdr$idl:fmd_zodi_skymap.pro  : fmd_zodi_skymap.pro  
csdr$idl:fmd_zsub_high_2.pro  : fmd_zsub_high_2.pro  
csdr$idl:fmd_zsub_high_3.pro  : fmd_zsub_high_3.pro  
csdr$idl:fmd_zsub_high_4.pro  : fmd_zsub_high_4.pro  
