!
!   Purpose: Build facility FSL_Spectra_Long
!
!   Author: Shirley M. Read, Hughes STX Corporation, 8/15/95, SPR 12288.
!
!   Modifications:
!

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fslbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fslbld.tlb")'" .eqs. "" then $(libr)/create/text fslbld
 if "''f$search("fslbld.olb")'" .eqs. "" then $(libr)/create fslbld

fslbld_txt = fsl_config,-
             fsl_display,-
             fsl_invoc,-
             fsl_model

fslbld_obj = fsl,-
             fsl_apply_model,-
             fsl_autophase_correct,-
             fsl_calc_responsivity,-
             fsl_calibrate_spectra,-
             fsl_close_config,-
             fsl_close_spectra,-
             fsl_compute_constants,-
             fsl_display_model,-
             fsl_display_spectra,-
             fsl_doppler_shift,-
             fsl_initialize_report,-
             fsl_open_cal_coadd,-
             fsl_open_coadd,-
             fsl_open_config,-
             fsl_open_sky_coadd,-
             fsl_open_spectra,-
             fsl_parse,-
             fsl_pcalib_variances,-
             fsl_produce_spectra,-
             fsl_read_cal_coadd,-
             fsl_read_coadd,-
             fsl_read_dvector,-
	     fsl_read_flv,-
             fsl_read_model,-
             fsl_read_reference,-
             fsl_read_sky_coadd,-
             fsl_read_sky_coadd_list,-
             fsl_temporal_drift,-
             fsl_update_report,-
             fsl_write_cal_spec,-
             fsl_write_sky_spec,-
             fsl_msg

fsl : fsl.exe
 ! Facility fsl has been built.

fsl.exe : fslbld.tlb($(fslbld_txt)),-
          fslbld.olb($(fslbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          xanadu:[lib]vilib.olb,-
          xanadu:[lib]xlib.olb,-
          graphics:grpshr.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fslbld/lib/inc=fsl,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      xanadu:[lib]vilib/lib, -
                      xanadu:[lib]xlib/lib, -
                      graphics:grpshr/lib, -
                      csdr$library:v5.opt/option

fslbld.tlb(fsl_config) : fsl_config.txt
 $(libr) $(librflags)/text fslbld fsl_config

fslbld.tlb(fsl_display) : fsl_display.txt
 $(libr) $(librflags)/text fslbld fsl_display

fslbld.tlb(fsl_invoc) : fsl_invoc.txt
 $(libr) $(librflags)/text fslbld fsl_invoc

fslbld.tlb(fsl_model) : fsl_model.txt
 $(libr) $(librflags)/text fslbld fsl_model

fsl.obj : fslbld.tlb(fsl_config),-
          fslbld.tlb(fsl_invoc),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctuser.inc,-
          fex_grtcoawt^,-
          fex_grttrans^,-
          fex_nyquistl^,-
          fex_vibcorrl^,-
          fil_sky^,-
          fsl_sky^

fsl_apply_model.obj : fslbld.tlb(fsl_config),-
                      fslbld.tlb(fsl_display),-
                      fslbld.tlb(fsl_invoc),-
                      fslbld.tlb(fsl_model),-
                      csdr$library:futlib.tlb(fut_params),-
                      fex_grtcoawt^,-
                      fex_grttrans^,-
                      fex_nyquistl^,-
                      fex_vibcorrl^

fsl_autophase_correct.obj : fslbld.tlb(fsl_model)

fsl_calc_responsivity.obj : fslbld.tlb(fsl_display),-
                            fslbld.tlb(fsl_model)

fsl_calibrate_spectra.obj : fslbld.tlb(fsl_config),-
                            fslbld.tlb(fsl_display),-
                            fslbld.tlb(fsl_invoc),-
                            fslbld.tlb(fsl_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquistl^,-
                            fex_vibcorrl^,-
                            fsl_sky^

fsl_close_config.obj : fslbld.tlb(fsl_config),-
                       fex_grtcoawt^,-
                       fex_grttrans^,-
                       fex_nyquistl^,-
                       fex_vibcorrl^

fsl_close_spectra.obj : fslbld.tlb(fsl_invoc),-
                        csdr$library:futlib.tlb(fut_params),-
                        csdr$library:ctuser.inc


fsl_compute_constants.obj : fslbld.tlb(fsl_config),-
                            fslbld.tlb(fsl_invoc),-
                            fslbld.tlb(fsl_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquistl^,-
                            fex_vibcorrl^

fsl_display_model.obj : fslbld.tlb(fsl_config),-
                        fslbld.tlb(fsl_invoc),-
                        fslbld.tlb(fsl_model),-
                        csdr$library:futlib.tlb(fut_params),-
                        fex_grtcoawt^,-
                        fex_grttrans^,-
                        fex_nyquistl^,-
                        fex_vibcorrl^

fsl_display_spectra.obj : fslbld.tlb(fsl_config),-
                          fslbld.tlb(fsl_display),-
                          fslbld.tlb(fsl_invoc),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fex_nyquistl^,-
                          fex_vibcorrl^

fsl_doppler_shift.obj : fslbld.tlb(fsl_invoc),-
                        fslbld.tlb(fsl_model),-
                        csdr$library:futlib.tlb(fut_params)

fsl_initialize_report.obj : fslbld.tlb(fsl_invoc),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fsl_open_cal_coadd.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                         csdr$library:csdrlib.tlb(cct_query_tod_catalog_record),-
                         ccm_cme_catalog_entry^

fsl_open_coadd.obj : fslbld.tlb(fsl_config),-
                     fslbld.tlb(fsl_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_grtcoawt^,-
                     fex_grttrans^,-
                     fex_nyquistl^,-
                     fex_vibcorrl^

fsl_open_config.obj : fslbld.tlb(fsl_config),-
                       fex_grtcoawt^,-
                       fex_grttrans^,-
                       fex_nyquistl^,-
                       fex_vibcorrl^

fsl_open_sky_coadd.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                         csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                         ccm_cme_catalog_entry^

fsl_open_spectra.obj : fslbld.tlb(fsl_config),-
                       fslbld.tlb(fsl_invoc),-
                       csdr$library:futlib.tlb(fut_params),-
                       fex_grtcoawt^,-
                       fex_grttrans^,-
                       fex_nyquistl^,-
                       fex_vibcorrl^


fsl_parse.obj : fslbld.tlb(fsl_invoc),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg)

fsl_pcalib_variances.obj : fslbld.tlb(fsl_model)

fsl_produce_spectra.obj : fslbld.tlb(fsl_config),-
                          fslbld.tlb(fsl_display),-
                          fslbld.tlb(fsl_invoc),-
                          fslbld.tlb(fsl_model),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fex_nyquistl^,-
                          fex_vibcorrl^,-
                          fil_sky^,-
                          fsl_sky^

fsl_read_cal_coadd.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fil_sky^

fsl_read_coadd.obj : fslbld.tlb(fsl_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:ctuser.inc,-
                     fil_sky^

fsl_read_dvector.obj : fslbld.tlb(fsl_invoc),-
                       fslbld.tlb(fsl_model),-
                       csdr$library:futlib.tlb(fut_params),-
                       fex_var^

fsl_read_flv.obj : fslbld.tlb(fsl_invoc),-
                   fslbld.tlb(fsl_model),-
                   csdr$library:futlib.tlb(fut_params),-
                   fex_flv^

fsl_read_model.obj : fslbld.tlb(fsl_invoc),-
                     fslbld.tlb(fsl_model),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_mod^

fsl_read_reference.obj : fslbld.tlb(fsl_config),-
                         fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         fex_grtcoawt^,-
                         fex_grttrans^,-
                         fex_nyquistl^,-
                         fex_vibcorrl^

fsl_read_sky_coadd.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                         csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                         fil_sky^

fsl_read_sky_coadd_list.obj : fslbld.tlb(fsl_invoc),-
                              csdr$library:futlib.tlb(fut_params),-
                              csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                              csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                              fil_sky^

fsl_temporal_drift.obj : fslbld.tlb(fsl_model)

fsl_update_report.obj : fslbld.tlb(fsl_invoc),-
                        fslbld.tlb(fsl_model),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params)

fsl_write_cal_spec.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fsl_sky^

fsl_write_sky_spec.obj : fslbld.tlb(fsl_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         fsl_sky^

