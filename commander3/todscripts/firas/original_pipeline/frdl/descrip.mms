! Delete MMS' suffix list and define a suffix list with only three entries.
! With this suffix list, MMS will be able to infer that both a DDL file and 
! an RDF file both depend on an RDL file of the same name.

.SUFFIXES
.SUFFIXES .DDL .RDF .RDL

! Define two default action lines so MMS can use it's new suffixes:
! define the RDL-to-RDF action and the RDL-to-DDL action.

.RDL.RDF
 RDC $(MMS$SOURCE)

.RDL.DDL
 RDC/DDL /OUTPUT=$(MMS$TARGET) $(MMS$SOURCE)

! Make sure the Record Definition Compiler's verb is defined.

.FIRST
 SET COMMAND CSDR$CLD:CRD

!===============================================================================
FRDL : DDL_TARGET, RDF_TARGET
 ! $(MMS$TARGET) is up to date

DDL_TARGET : -
      fad_sky.ddl,-
      fcc_cov.ddl,-
      fcf_cal.ddl,-
      fcf_dcl.ddl,-
      fcf_dsk.ddl,-
      fcf_sky.ddl,-
      fcf_vcl.ddl,-
      fcf_vsk.ddl,-
      fcl_cov.ddl,-
      fcs_sky.ddl,-
      fdq_eng.ddl,-
      fdq_etr.ddl,-
      fdq_idx.ddl,-
      fdq_sdf.ddl,-
      fec_sscal.ddl,-
      fef_spc.ddl,-
      fel_err.ddl,-
      fep_stoldb.ddl,-      
      fer_err.ddl,-      
      fex_apod.ddl,-
      fex_apodl.ddl,-
      fex_av_calrs.ddl,-
      fex_basis.ddl,-
      fex_calres.ddl,-
      fex_cmdgain.ddl,-
      fex_cth.ddl,-
      fex_cvs.ddl,-
      fex_dtf.ddl,-
      fex_dtrf.ddl,-
      fex_ejv.ddl,-
      fex_englim.ddl,-
      fex_etf.ddl,-
      fex_etfl.ddl,-
      fex_fakeit.ddl,-
      fex_flv.ddl,-
      fex_gain.ddl,-
      fex_gltchcor.ddl,-
      fex_gltchpro.ddl,-
      fex_grtcoawt.ddl,-
      fex_grtrawwt.ddl,-
      fex_grttrans.ddl,-
      fex_idx_flag.ddl,-
      fex_idx_tols.ddl,-
      fex_limflags.ddl,-
      fex_mcs.ddl,-
      fex_mincoadd.ddl,-
      fex_mod.ddl,-
      fex_mtmsweep.ddl,-
      fex_nyquist.ddl,-
      fex_nyquistl.ddl,-
      fex_reftemps.ddl,-
      fex_samprate.ddl,-
      fex_scilim.ddl,-
      fex_vabsaa.ddl,-
      fex_var.ddl,-
      fex_vibcorr.ddl,-
      fex_vibcorrl.ddl,-
      ffp_ifg.ddl,-
      ffp_spc.ddl,-
      fic_cal.ddl,-
      fic_ccc.ddl,-
      fic_cco.ddl,-
      fic_cdg.ddl,-
      fic_cov.ddl,-
      fic_cpt.ddl,-
      fic_cst.ddl,-
      fic_scc.ddl,-
      fic_sco.ddl,-
      fic_sdg.ddl,-
      fic_sky.ddl,-
      fic_spt.ddl,-
      fic_sst.ddl,-
      fil_cov.ddl,-
      fil_scc.ddl,-
      fil_sky.ddl,-
      fip_ccl.ddl,-
      fip_cov.ddl,-
      fip_csk.ddl,-
      fip_csq.ddl,-
      fip_cvs.ddl,-
      fip_dcl.ddl,-
      fip_dse.ddl,-
      fip_dsk.ddl,-
      fip_dsq.ddl,-
      fip_dst.ddl,-
      fip_ejv.ddl,-
      fip_err.ddl,-
      fip_fef.ddl,-
      fip_hlp.ddl,-
      fip_icl.ddl,-
      fip_isk.ddl,-
      fip_llp.ddl,-
      fip_lmh.ddl,-
      fip_lml.ddl,-
      fip_mod.ddl,-
      fip_sky.ddl,-
      fip_tcb.ddl,-
      fla_dst.ddl,-
      fms_cvs.ddl,-
      fms_sky.ddl,-
      fms_var.ddl,-
      fnt_noise.ddl,-
      fpp_sdf.ddl,-
      fsd_sky.ddl,-
      fsl_sky.ddl,-
      fss_sssky.ddl,-
      fut_attit.ddl,-
      fut_enganlg.ddl,-
      fut_engstat.ddl,-
      fxt_eng_xtrm.ddl
 ! $(MMS$TARGET) is up to date

RDF_TARGET : -
      fad_sky.rdf,-
      fcc_cov.rdf,-
      fcf_cal.rdf,-
      fcf_dcl.rdf,-
      fcf_dsk.rdf,-
      fcf_sky.rdf,-
      fcf_vcl.rdf,-
      fcf_vsk.rdf,-
      fcl_cov.rdf,-
      fcs_sky.rdf,-
      fdq_eng.rdf,-
      fdq_etr.rdf,-
      fdq_idx.rdf,-
      fdq_sdf.rdf,-
      fec_sscal.rdf,-
      fef_spc.rdf,-
      fel_err.rdf,-
      fep_stoldb.rdf,-      
      fer_err.rdf,-
      fex_apod.rdf,-
      fex_apodl.rdf,-
      fex_av_calrs.rdf,-
      fex_basis.rdf,-
      fex_calres.rdf,-
      fex_cmdgain.rdf,-
      fex_cth.rdf,-
      fex_cvs.rdf,-
      fex_dtf.rdf,-
      fex_dtrf.rdf,-
      fex_ejv.rdf,-
      fex_englim.rdf,-
      fex_etf.rdf,-
      fex_etfl.rdf,-
      fex_fakeit.rdf,-
      fex_flv.rdf,-
      fex_gain.rdf,-
      fex_gltchcor.rdf,-
      fex_gltchpro.rdf,-
      fex_grtcoawt.rdf,-
      fex_grtrawwt.rdf,-
      fex_grttrans.rdf,-
      fex_idx_flag.rdf,-
      fex_idx_tols.rdf,-
      fex_limflags.rdf,-
      fex_mcs.rdf,-
      fex_mincoadd.rdf,-
      fex_mod.rdf,-
      fex_mtmsweep.rdf,-
      fex_nyquist.rdf,-
      fex_nyquistl.rdf,-
      fex_reftemps.rdf,-
      fex_samprate.rdf,-
      fex_scilim.rdf,-
      fex_vabsaa.rdf,-
      fex_var.rdf,-
      fex_vibcorr.rdf,-
      fex_vibcorrl.rdf,-
      ffp_ifg.rdf,-
      ffp_spc.rdf,-
      fic_cal.rdf,-
      fic_ccc.rdf,-
      fic_cco.rdf,-
      fic_cdg.rdf,-
      fic_cov.rdf,-
      fic_cpt.rdf,-
      fic_cst.rdf,-
      fic_scc.rdf,-
      fic_sco.rdf,-
      fic_sdg.rdf,-
      fic_sky.rdf,-
      fic_spt.rdf,-
      fic_sst.rdf,-
      fil_cov.rdf,-
      fil_scc.rdf,-
      fil_sky.rdf,-
      fip_ccl.rdf,-
      fip_cov.rdf,-
      fip_csk.rdf,-
      fip_csq.rdf,-
      fip_cvs.rdf,-
      fip_dcl.rdf,-
      fip_dse.rdf,-
      fip_dsk.rdf,-
      fip_dsq.rdf,-
      fip_dst.rdf,-
      fip_ejv.rdf,-
      fip_err.rdf,-
      fip_fef.rdf,-
      fip_hlp.rdf,-
      fip_icl.rdf,-
      fip_isk.rdf,-
      fip_llp.rdf,-
      fip_lmh.rdf,-
      fip_lml.rdf,-
      fip_mod.rdf,-
      fip_sky.rdf,-
      fip_tcb.rdf,-
      fla_dst.rdf,-
      fms_cvs.rdf,-
      fms_sky.rdf,-
      fms_var.rdf,-
      fnt_noise.rdf,-
      fpp_sdf.rdf,-
      fsd_sky.rdf,-
      fsl_sky.rdf,-
      fss_sssky.rdf,-
      fut_attit.rdf,-
      fut_enganlg.rdf,-
      fut_engstat.rdf,-
      fxt_eng_xtrm.rdf
 ! $(MMS$TARGET) is up to date

! There are some include-file dependencies to be considered.

fad_sky.ddl, fad_sky.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_cal.ddl, fcf_cal.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_dcl.ddl, fcf_dcl.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_dsk.ddl, fcf_dsk.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_sky.ddl, fcf_sky.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_vcl.ddl, fcf_vcl.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcf_vsk.ddl, fcf_vsk.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fcs_sky.ddl, fcs_sky.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fdq_eng.ddl, fdq_eng.rdf : ct_head.rdl,-
                           eng_head.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_xcal.rdl,-
                           eng_setup.rdl,-
                           eng_tempdiff.rdl,-
                           eng_tail.rdl
fdq_etr.ddl, fdq_etr.rdf : ct_head.rdl,-
                           eng_head.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_xcal.rdl,-
                           eng_setup.rdl,-
                           eng_tempdiff.rdl,-
                           eng_tail.rdl
fdq_idx.ddl, fdq_idx.rdf : index_head.rdl,-
                           idx_xcal.rdl,-
                           idx_chan_data.rdl,-
                           idx_grt_data.rdl,-
                           temp_ctrl.rdl,-
                           status_monitor.rdl, -
                           idx_sync.rdl,-
                           ipdu_power_relay.rdl,-
                           ipdu_status.rdl,-
                           misc_stat.rdl,-
                           idx_analog.rdl,-
                           idx_tail.rdl
fdq_sdf.ddl, fdq_sdf.rdf : ct_head.rdl,-
                           sci_head.rdl,-
                           ifg_data.rdl,-
                           dq_data.rdl,-
                           collect_time.rdl,-
                           attitude.rdl
fex_av_calrs.ddl, fex_av_calrs.rdf : ct_head.rdl
fex_calres.ddl, fex_calres.rdf : ct_head.rdl
fex_cmdgain.ddl, fex_cmdgain.rdf : ct_head.rdl
fex_cth.ddl, fex_cth.rdf : ct_head.rdl
fex_ejv.ddl, fex_ejv.rdf : ct_head.rdl
fex_englim.ddl, fex_englim.rdf : ct_head.rdl,-
                                 eng_head.rdl,-
                                 eng_status.rdl,-
                                 eng_analog.rdl,-
                                 eng_xcal.rdl,-
                                 eng_setup.rdl,-
                                 eng_tempdiff.rdl,-
                                 eng_tail.rdl
fex_fakeit.ddl, fex_fakeit.rdf : ct_head.rdl
fex_gain.ddl, fex_gain.rdf : ct_head.rdl
fex_grtcoawt.ddl, fex_grtcoawt.rdf : ct_head.rdl
fex_grtrawwt.ddl, fex_grtrawwt.rdf : ct_head.rdl
fex_grttrans.ddl, fex_grttrans.rdf : ct_head.rdl
fex_idx_flag.ddl, fex_idx_flag.rdf : ct_head.rdl,-
                                     idx_flags.rdl
fex_idx_tols.ddl, fex_idx_tols.rdf : ct_head.rdl,-
                                     idx_tols.rdl
fex_limflags.ddl, fex_limflags.rdf : ct_head.rdl,-
                                     lim_flags.rdl
fex_mincoadd.ddl, fex_mincoadd.rdf : ct_head.rdl
fex_mod.ddl, fex_mod.rdf : ct_head.rdl,-
                           mod_head.rdl
fex_mtmsweep.ddl, fex_mtmsweep.rdf : ct_head.rdl
fex_scilim.ddl, fex_scilim.rdf : ct_head.rdl,-
                                 sci_head.rdl,-
                                 dq_data.rdl,-
                                 collect_time.rdl,-
                                 attitude.rdl
fex_vabsaa.ddl, fex_vabsaa.rdf : ct_head.rdl
fex_vibcorr.ddl, fex_vibcorr.rdf : ct_head.rdl
fex_vibcorrl.ddl, fex_vibcorrl.rdf : ct_head.rdl
fic_cal.ddl, fic_cal.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_ccc.ddl, fic_ccc.rdf : con_check.rdl
fic_cco.ddl, fic_cco.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_cdg.ddl, fic_cdg.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_cov.ddl, fic_cov.rdf : ct_head.rdl
fic_cpt.ddl, fic_cpt.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_cst.ddl, fic_cst.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_scc.ddl, fic_scc.rdf : con_check.rdl
fic_sco.ddl, fic_sco.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_sdg.ddl, fic_sdg.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_sky.ddl, fic_sky.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_spt.ddl, fic_spt.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fic_sst.ddl, fic_sst.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           coad_data.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fil_cov.ddl, fil_cov.rdf : ct_head.rdl
fil_scc.ddl, fil_scc.rdf : con_checkl.rdl
fil_sky.ddl, fil_sky.rdf : ct_head.rdl,-
                           coad_spec_headl.rdl,-
                           coad_spec_datal.rdl,-
                           coad_datal.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fms_sky.ddl, fms_sky.rdf : ct_head.rdl,-
                           coad_spec_head.rdl,-
                           coad_spec_data.rdl,-
                           comb_spec_data.rdl,-
                           combinations.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fnt_noise.ddl, fnt_noise.rdf : ct_head.rdl,-
                               spec_head.rdl,-
                               chan_data.rdl,-
                               attitude.rdl
fpp_sdf.ddl, fpp_sdf.rdf : ct_head.rdl,-
                           sci_head.rdl,-
                           ifg_data.rdl,-
                           dq_data.rdl,-
                           collect_time.rdl,-
                           attitude.rdl
fsl_sky.ddl, fsl_sky.rdf : ct_head.rdl,-
                           coad_spec_headl.rdl,-
                           coad_spec_datal.rdl,-
                           spec_datal.rdl,-
                           eng_status.rdl,-
                           eng_analog.rdl,-
                           eng_sigma.rdl,-
                           eng_tempdiff.rdl,-
                           attitude.rdl
fut_attit.ddl, fut_attit.rdf : attitude.rdl
fut_enganlg.ddl, fut_enganlg.rdf : eng_analog.rdl
fut_engstat.ddl, fut_engstat.rdf : eng_status.rdl
fxt_eng_xtrm.ddl, fxt_eng_xtrm.rdf : ct_head.rdl
