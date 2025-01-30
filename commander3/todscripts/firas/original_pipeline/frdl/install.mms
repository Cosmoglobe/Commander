.DEFAULT
 COPY $(MMS$SOURCE) $(MMS$TARGET);0

FRDL : RDF_TARGET, CDD_TARGET
 ! $(MMS$TARGET) is up to date.

RDF_TARGET : -
      csdr$rdf:fad_sky.rdf,-
      csdr$rdf:fcc_cov.rdf,-
      csdr$rdf:fcf_cal.rdf,-
      csdr$rdf:fcf_dcl.rdf,-
      csdr$rdf:fcf_dsk.rdf,-
      csdr$rdf:fcf_sky.rdf,-
      csdr$rdf:fcf_vcl.rdf,-
      csdr$rdf:fcf_vsk.rdf,-
      csdr$rdf:fcl_cov.rdf,-
      csdr$rdf:fcs_sky.rdf,-
      csdr$rdf:fdq_eng.rdf,-
      csdr$rdf:fdq_etr.rdf,-
      csdr$rdf:fdq_idx.rdf,-
      csdr$rdf:fdq_sdf.rdf,-
      csdr$rdf:fec_sscal.rdf,-
      csdr$rdf:fef_spc.rdf,-
      csdr$rdf:fel_err.rdf,-
      csdr$rdf:fep_stoldb.rdf,-
      csdr$rdf:fer_err.rdf,-
      csdr$rdf:fex_apod.rdf,-
      csdr$rdf:fex_apodl.rdf,-
      csdr$rdf:fex_av_calrs.rdf,-
      csdr$rdf:fex_basis.rdf,-
      csdr$rdf:fex_calres.rdf,-
      csdr$rdf:fex_cmdgain.rdf,-
      csdr$rdf:fex_cth.rdf,-
      csdr$rdf:fex_cvs.rdf,-
      csdr$rdf:fex_dtf.rdf,-
      csdr$rdf:fex_dtrf.rdf,-
      csdr$rdf:fex_ejv.rdf,-
      csdr$rdf:fex_englim.rdf,-
      csdr$rdf:fex_etf.rdf,-
      csdr$rdf:fex_etfl.rdf,-
      csdr$rdf:fex_fakeit.rdf,-
      csdr$rdf:fex_flv.rdf,-
      csdr$rdf:fex_gain.rdf,-
      csdr$rdf:fex_gltchcor.rdf,-
      csdr$rdf:fex_gltchpro.rdf,-
      csdr$rdf:fex_grtcoawt.rdf,-
      csdr$rdf:fex_grtrawwt.rdf,-
      csdr$rdf:fex_grttrans.rdf,-
      csdr$rdf:fex_idx_flag.rdf,-
      csdr$rdf:fex_idx_tols.rdf,-
      csdr$rdf:fex_limflags.rdf,-
      csdr$rdf:fex_mcs.rdf,-
      csdr$rdf:fex_mincoadd.rdf,-
      csdr$rdf:fex_mod.rdf,-
      csdr$rdf:fex_mtmsweep.rdf,-
      csdr$rdf:fex_nyquist.rdf,-
      csdr$rdf:fex_nyquistl.rdf,-
      csdr$rdf:fex_reftemps.rdf,-
      csdr$rdf:fex_samprate.rdf,-
      csdr$rdf:fex_scilim.rdf,-
      csdr$rdf:fex_vabsaa.rdf,-
      csdr$rdf:fex_var.rdf,-
      csdr$rdf:fex_vibcorr.rdf,-
      csdr$rdf:fex_vibcorrl.rdf,-
      csdr$rdf:ffp_ifg.rdf,-
      csdr$rdf:ffp_spc.rdf,-
      csdr$rdf:fic_cal.rdf,-
      csdr$rdf:fic_ccc.rdf,-
      csdr$rdf:fic_cco.rdf,-
      csdr$rdf:fic_cdg.rdf,-
      csdr$rdf:fic_cov.rdf,-
      csdr$rdf:fic_cpt.rdf,-
      csdr$rdf:fic_cst.rdf,-
      csdr$rdf:fic_scc.rdf,-
      csdr$rdf:fic_sco.rdf,-
      csdr$rdf:fic_sdg.rdf,-
      csdr$rdf:fic_sky.rdf,-
      csdr$rdf:fic_spt.rdf,-
      csdr$rdf:fic_sst.rdf,-
      csdr$rdf:fil_cov.rdf,-
      csdr$rdf:fil_scc.rdf,-
      csdr$rdf:fil_sky.rdf,-
      csdr$rdf:fip_ccl.rdf,-
      csdr$rdf:fip_cov.rdf,-
      csdr$rdf:fip_csk.rdf,-
      csdr$rdf:fip_csq.rdf,-
      csdr$rdf:fip_cvs.rdf,-
      csdr$rdf:fip_dcl.rdf,-
      csdr$rdf:fip_dse.rdf,-
      csdr$rdf:fip_dsk.rdf,-
      csdr$rdf:fip_dsq.rdf,-
      csdr$rdf:fip_dst.rdf,-
      csdr$rdf:fip_ejv.rdf,-
      csdr$rdf:fip_err.rdf,-
      csdr$rdf:fip_fef.rdf,-
      csdr$rdf:fip_hlp.rdf,-
      csdr$rdf:fip_icl.rdf,-
      csdr$rdf:fip_isk.rdf,-
      csdr$rdf:fip_llp.rdf,-
      csdr$rdf:fip_lmh.rdf,-
      csdr$rdf:fip_lml.rdf,-
      csdr$rdf:fip_mod.rdf,-
      csdr$rdf:fip_sky.rdf,-
      csdr$rdf:fip_tcb.rdf,-
      csdr$rdf:fla_dst.rdf,-
      csdr$rdf:fms_cvs.rdf,-      
      csdr$rdf:fms_sky.rdf,-
      csdr$rdf:fms_var.rdf,-      
      csdr$rdf:fnt_noise.rdf,-
      csdr$rdf:fpp_sdf.rdf,-
      csdr$rdf:fsd_sky.rdf,-
      csdr$rdf:fsl_sky.rdf,-
      csdr$rdf:fss_sssky.rdf,-
      csdr$rdf:fut_attit.rdf,-
      csdr$rdf:fut_enganlg.rdf,-
      csdr$rdf:fut_engstat.rdf,-
      csdr$rdf:fxt_eng_xtrm.rdf
 ! $(MMS$TARGET) is up to date.

CDD_TARGET : -
      fad_sky^,-
      fcc_cov^,-
      fcf_cal^,-
      fcf_dcl^,-
      fcf_dsk^,-
      fcf_sky^,-
      fcf_vcl^,-
      fcf_vsk^,-
      fcl_cov^,-
      fcs_sky^,-
      fdq_eng^,-
      fdq_etr^,-
      fdq_idx^,-
      fdq_sdf^,-
      fec_sscal^,-
      fef_spc^,-
      fel_err^,-
      fep_stoldb^,-
      fer_err^,-
      fex_apod^,-
      fex_apodl^,-
      fex_av_calrs^,-
      fex_basis^,-
      fex_calres^,-
      fex_cmdgain^,-
      fex_cth^,-
      fex_cvs^,-
      fex_dtf^,-
      fex_dtrf^,-
      fex_ejv^,-
      fex_englim^,-
      fex_etf^,-
      fex_etfl^,-
      fex_fakeit^,-
      fex_flv^,-
      fex_gain^,-
      fex_gltchcor^,-
      fex_gltchpro^,-
      fex_grtcoawt^,-
      fex_grtrawwt^,-
      fex_grttrans^,-
      fex_idx_flag^,-
      fex_idx_tols^,-
      fex_limflags^,-
      fex_mcs^,-
      fex_mincoadd^,-
      fex_mod^,-
      fex_mtmsweep^,-
      fex_nyquist^,-
      fex_nyquistl^,-
      fex_reftemps^,-
      fex_samprate^,-
      fex_scilim^,-
      fex_vabsaa^,-
      fex_var^,-
      fex_vibcorr^,-
      fex_vibcorrl^,-
      ffp_ifg^,-
      ffp_spc^,-
      fic_cal^,-
      fic_ccc^,-
      fic_cco^,-
      fic_cdg^,-
      fic_cov^,-
      fic_cpt^,-
      fic_cst^,-
      fic_scc^,-
      fic_sco^,-
      fic_sdg^,-
      fic_sky^,-
      fic_spt^,-
      fic_sst^,-
      fil_cov^,-
      fil_scc^,-
      fil_sky^,-
      fip_ccl^,-
      fip_cov^,-
      fip_csk^,-
      fip_csq^,-
      fip_cvs^,-
      fip_dcl^,-
      fip_dse^,-
      fip_dsk^,-
      fip_dsq^,-
      fip_dst^,-
      fip_ejv^,-
      fip_err^,-
      fip_fef^,-
      fip_hlp^,-
      fip_icl^,-
      fip_isk^,-
      fip_llp^,-
      fip_lmh^,-
      fip_lml^,-
      fip_mod^,-
      fip_sky^,-
      fip_tcb^,-
      fla_dst^,-
      fms_cvs^,-
      fms_sky^,-
      fms_var^,-
      fnt_noise^,-
      fpp_sdf^,-
      fsd_sky^,-
      fsl_sky^,-
      fss_sssky^,-
      fut_attit^,-
      fut_enganlg^,-
      fut_engstat^,-
      fxt_eng_xtrm^
 ! $(MMS$TARGET) is up to date.

csdr$rdf:fad_sky.rdf : fad_sky.rdf
csdr$rdf:fcc_cov.rdf : fcc_cov.rdf
csdr$rdf:fcf_cal.rdf : fcf_cal.rdf
csdr$rdf:fcf_dcl.rdf : fcf_dcl.rdf
csdr$rdf:fcf_dsk.rdf : fcf_dsk.rdf
csdr$rdf:fcf_sky.rdf : fcf_sky.rdf
csdr$rdf:fcf_vcl.rdf : fcf_vcl.rdf
csdr$rdf:fcf_vsk.rdf : fcf_vsk.rdf
csdr$rdf:fcl_cov.rdf : fcl_cov.rdf
csdr$rdf:fcs_sky.rdf : fcs_sky.rdf
csdr$rdf:fdq_eng.rdf : fdq_eng.rdf
csdr$rdf:fdq_etr.rdf : fdq_etr.rdf
csdr$rdf:fdq_idx.rdf : fdq_idx.rdf
csdr$rdf:fdq_sdf.rdf : fdq_sdf.rdf
csdr$rdf:fec_sscal.rdf : fec_sscal.rdf
csdr$rdf:fef_spc.rdf : fef_spc.rdf
csdr$rdf:fel_err.rdf : fel_err.rdf
csdr$rdf:fep_stoldb.rdf : fep_stoldb.rdf
csdr$rdf:fer_err.rdf : fer_err.rdf
csdr$rdf:fex_apod.rdf : fex_apod.rdf
csdr$rdf:fex_apodl.rdf : fex_apodl.rdf
csdr$rdf:fex_av_calrs.rdf : fex_av_calrs.rdf
csdr$rdf:fex_basis.rdf : fex_basis.rdf
csdr$rdf:fex_calres.rdf : fex_calres.rdf
csdr$rdf:fex_cmdgain.rdf : fex_cmdgain.rdf
csdr$rdf:fex_cth.rdf : fex_cth.rdf
csdr$rdf:fex_cvs.rdf : fex_cvs.rdf
csdr$rdf:fex_dtf.rdf : fex_dtf.rdf
csdr$rdf:fex_dtrf.rdf : fex_dtrf.rdf
csdr$rdf:fex_ejv.rdf : fex_ejv.rdf
csdr$rdf:fex_englim.rdf : fex_englim.rdf
csdr$rdf:fex_etf.rdf : fex_etf.rdf
csdr$rdf:fex_etfl.rdf : fex_etfl.rdf
csdr$rdf:fex_fakeit.rdf : fex_fakeit.rdf
csdr$rdf:fex_flv.rdf : fex_flv.rdf
csdr$rdf:fex_gain.rdf : fex_gain.rdf
csdr$rdf:fex_gltchcor.rdf : fex_gltchcor.rdf
csdr$rdf:fex_gltchpro.rdf : fex_gltchpro.rdf
csdr$rdf:fex_grtcoawt.rdf : fex_grtcoawt.rdf
csdr$rdf:fex_grtrawwt.rdf : fex_grtrawwt.rdf
csdr$rdf:fex_grttrans.rdf : fex_grttrans.rdf
csdr$rdf:fex_idx_flag.rdf : fex_idx_flag.rdf
csdr$rdf:fex_idx_tols.rdf : fex_idx_tols.rdf
csdr$rdf:fex_limflags.rdf : fex_limflags.rdf
csdr$rdf:fex_mcs.rdf : fex_mcs.rdf
csdr$rdf:fex_mincoadd.rdf : fex_mincoadd.rdf
csdr$rdf:fex_mod.rdf : fex_mod.rdf
csdr$rdf:fex_mtmsweep.rdf : fex_mtmsweep.rdf
csdr$rdf:fex_nyquist.rdf : fex_nyquist.rdf
csdr$rdf:fex_nyquistl.rdf : fex_nyquistl.rdf
csdr$rdf:fex_reftemps.rdf : fex_reftemps.rdf
csdr$rdf:fex_samprate.rdf : fex_samprate.rdf
csdr$rdf:fex_scilim.rdf : fex_scilim.rdf
csdr$rdf:fex_vabsaa.rdf : fex_vabsaa.rdf
csdr$rdf:fex_var.rdf : fex_var.rdf
csdr$rdf:fex_vibcorr.rdf : fex_vibcorr.rdf
csdr$rdf:fex_vibcorrl.rdf : fex_vibcorrl.rdf
csdr$rdf:ffp_ifg.rdf : ffp_ifg.rdf
csdr$rdf:ffp_spc.rdf : ffp_spc.rdf
csdr$rdf:fic_cal.rdf : fic_cal.rdf
csdr$rdf:fic_ccc.rdf : fic_ccc.rdf
csdr$rdf:fic_cco.rdf : fic_cco.rdf
csdr$rdf:fic_cdg.rdf : fic_cdg.rdf
csdr$rdf:fic_cov.rdf : fic_cov.rdf
csdr$rdf:fic_cpt.rdf : fic_cpt.rdf
csdr$rdf:fic_cst.rdf : fic_cst.rdf
csdr$rdf:fic_scc.rdf : fic_scc.rdf
csdr$rdf:fic_sco.rdf : fic_sco.rdf
csdr$rdf:fic_sdg.rdf : fic_sdg.rdf
csdr$rdf:fic_sky.rdf : fic_sky.rdf
csdr$rdf:fic_spt.rdf : fic_spt.rdf
csdr$rdf:fic_sst.rdf : fic_sst.rdf
csdr$rdf:fil_cov.rdf : fil_cov.rdf
csdr$rdf:fil_scc.rdf : fil_scc.rdf
csdr$rdf:fil_sky.rdf : fil_sky.rdf
csdr$rdf:fip_ccl.rdf : fip_ccl.rdf
csdr$rdf:fip_cov.rdf : fip_cov.rdf
csdr$rdf:fip_csk.rdf : fip_csk.rdf
csdr$rdf:fip_csq.rdf : fip_csq.rdf
csdr$rdf:fip_cvs.rdf : fip_cvs.rdf
csdr$rdf:fip_dcl.rdf : fip_dcl.rdf
csdr$rdf:fip_dse.rdf : fip_dse.rdf
csdr$rdf:fip_dsk.rdf : fip_dsk.rdf
csdr$rdf:fip_dsq.rdf : fip_dsq.rdf
csdr$rdf:fip_dst.rdf : fip_dst.rdf
csdr$rdf:fip_ejv.rdf : fip_ejv.rdf
csdr$rdf:fip_err.rdf : fip_err.rdf
csdr$rdf:fip_fef.rdf : fip_fef.rdf
csdr$rdf:fip_hlp.rdf : fip_hlp.rdf
csdr$rdf:fip_icl.rdf : fip_icl.rdf
csdr$rdf:fip_isk.rdf : fip_isk.rdf
csdr$rdf:fip_llp.rdf : fip_llp.rdf
csdr$rdf:fip_lmh.rdf : fip_lmh.rdf
csdr$rdf:fip_lml.rdf : fip_lml.rdf
csdr$rdf:fip_mod.rdf : fip_mod.rdf
csdr$rdf:fip_sky.rdf : fip_sky.rdf
csdr$rdf:fip_tcb.rdf : fip_tcb.rdf
csdr$rdf:fla_dst.rdf : fla_dst.rdf
csdr$rdf:fms_cvs.rdf : fms_cvs.rdf
csdr$rdf:fms_sky.rdf : fms_sky.rdf
csdr$rdf:fms_var.rdf : fms_var.rdf
csdr$rdf:fnt_noise.rdf : fnt_noise.rdf
csdr$rdf:fpp_sdf.rdf : fpp_sdf.rdf
csdr$rdf:fsd_sky.rdf : fsd_sky.rdf
csdr$rdf:fsl_sky.rdf : fsl_sky.rdf
csdr$rdf:fss_sssky.rdf : fss_sssky.rdf
csdr$rdf:fut_attit.rdf : fut_attit.rdf
csdr$rdf:fut_enganlg.rdf : fut_enganlg.rdf
csdr$rdf:fut_engstat.rdf : fut_engstat.rdf
csdr$rdf:fxt_eng_xtrm.rdf : fxt_eng_xtrm.rdf

fad_sky^ : fad_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcc_cov^ : fcc_cov.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_cal^ : fcf_cal.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_dcl^ : fcf_dcl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_dsk^ : fcf_dsk.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_sky^ : fcf_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_vcl^ : fcf_vcl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcf_vsk^ : fcf_vsk.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcl_cov^ : fcl_cov.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fcs_sky^ : fcs_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fdq_eng^ : fdq_eng.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fdq_etr^ : fdq_etr.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fdq_idx^ : fdq_idx.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fdq_sdf^ : fdq_sdf.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fec_sscal^ : fec_sscal.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fef_spc^ : fef_spc.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fel_err^ : fel_err.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fep_stoldb^ : fep_stoldb.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fer_err^ : fer_err.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_apod^ : fex_apod.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_apodl^ : fex_apodl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_av_calrs^ : fex_av_calrs.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_basis^ : fex_basis.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_calres^ : fex_calres.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_cmdgain^ : fex_cmdgain.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_cth^ : fex_cth.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_cvs^ : fex_cvs.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_dtf^ : fex_dtf.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_dtrf^ : fex_dtrf.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_ejv^ : fex_ejv.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_englim^ : fex_englim.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_etf^ : fex_etf.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_etfl^ : fex_etfl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_fakeit^ : fex_fakeit.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_flv^ : fex_flv.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_gain^ : fex_gain.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_gltchcor^ : fex_gltchcor.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_gltchpro^ : fex_gltchpro.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_grtcoawt^ : fex_grtcoawt.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_grtrawwt^ : fex_grtrawwt.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_grttrans^ : fex_grttrans.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_idx_flag^ : fex_idx_flag.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_idx_tols^ : fex_idx_tols.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_limflags^ : fex_limflags.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_mcs^ : fex_mcs.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_mincoadd^ : fex_mincoadd.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_mod^ : fex_mod.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_mtmsweep^ : fex_mtmsweep.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_nyquist^ : fex_nyquist.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_nyquistl^ : fex_nyquistl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_reftemps^ : fex_reftemps.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_samprate^ : fex_samprate.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_scilim^ : fex_scilim.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_vabsaa^ : fex_vabsaa.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_var^ : fex_var.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_vibcorr^ : fex_vibcorr.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fex_vibcorrl^ : fex_vibcorrl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
ffp_ifg^ : ffp_ifg.ddl
 CDDL/REPLACE $(MMS$SOURCE)
ffp_spc^ : ffp_spc.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cal^ : fic_cal.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_ccc^ : fic_ccc.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cco^ : fic_cco.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cdg^ : fic_cdg.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cov^ : fic_cov.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cpt^ : fic_cpt.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_cst^ : fic_cst.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_scc^ : fic_scc.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_sco^ : fic_sco.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_sdg^ : fic_sdg.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_sky^ : fic_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_spt^ : fic_spt.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fic_sst^ : fic_sst.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fil_cov^ : fil_cov.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fil_scc^ : fil_scc.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fil_sky^ : fil_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_ccl^ : fip_ccl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_cov^ : fip_cov.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_csk^ : fip_csk.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_csq^ : fip_csq.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_cvs^ : fip_cvs.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_dcl^ : fip_dcl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_dse^ : fip_dse.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_dsk^ : fip_dsk.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_dsq^ : fip_dsq.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_dst^ : fip_dst.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_ejv^ : fip_ejv.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_err^ : fip_err.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_fef^ : fip_fef.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_hlp^ : fip_hlp.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_icl^ : fip_icl.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_isk^ : fip_isk.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_llp^ : fip_llp.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_lmh^ : fip_lmh.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_lml^ : fip_lml.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_mod^ : fip_mod.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_sky^ : fip_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fip_tcb^ : fip_tcb.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fla_dst^ : fla_dst.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fms_cvs^ : fms_cvs.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fms_sky^ : fms_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fms_var^ : fms_var.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fnt_noise^ : fnt_noise.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fpp_sdf^ : fpp_sdf.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fsd_sky^ : fsd_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fsl_sky^ : fsl_sky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fss_sssky^ : fss_sssky.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fut_attit^ : fut_attit.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fut_enganlg^ : fut_enganlg.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fut_engstat^ : fut_engstat.ddl
 CDDL/REPLACE $(MMS$SOURCE)
fxt_eng_xtrm^ : fxt_eng_xtrm.ddl
 CDDL/REPLACE $(MMS$SOURCE)
