!
!   Purpose: Build facility fcf_calibrate_firas.
!
!   Author: S. Alexander, STX, 3/19/92, SER 8292
!
!   Modifications:
!
!   Added fcf_read_dvector to allow reading of the FEX_VAR
!       referene dataset.
!   Gene Eplee, GSC, 25 October 1993
!   SER 11397
!

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fcfbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fcfbld.tlb")'" .eqs. "" then $(libr)/create/text fcfbld
 if "''f$search("fcfbld.olb")'" .eqs. "" then $(libr)/create fcfbld

fcfbld_txt = fcf_config,-
             fcf_display,-
             fcf_invoc,-
             fcf_model

fcfbld_obj = fcf,-
             fcf_apodize_and_rotate,-
             fcf_apply_model,-
             fcf_autophase_correct,-
             fcf_calc_responsivity,-
             fcf_calibrate_spectra,-
             fcf_compute_constants,-
             fcf_display_model,-
             fcf_display_spectra,-
             fcf_doppler_shift,-
             fcf_initialize_report,-
             fcf_open_cal_coadd,-
             fcf_open_coadd,-
             fcf_open_sky_coadd,-
             fcf_open_temp_spec,-
             fcf_parse,-
             fcf_pcalib_variances,-
             fcf_produce_spectra,-
             fcf_read_cal_coadd,-
             fcf_read_coadd,-
             fcf_read_dvector,-
             fcf_read_model,-
             fcf_read_reference,-
             fcf_read_sky_coadd,-
             fcf_read_sky_coadd_list,-
             fcf_temporal_drift,-
             fcf_update_report,-
             fcf_write_cal_spec,-
             fcf_write_sky_spec,-
             fcf_write_temp_spec,-
             fcf_msg

fcf : fcf.exe
 ! Facility fcf has been built.

fcf.exe : fcfbld.tlb($(fcfbld_txt)),-
          fcfbld.olb($(fcfbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          xanadu:[lib]vilib.olb,-
          xanadu:[lib]xlib.olb,-
          graphics:grpshr.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fcfbld/lib/inc=fcf,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      xanadu:[lib]vilib/lib, -
                      xanadu:[lib]xlib/lib, -
                      graphics:grpshr/lib, -
                      csdr$library:v5.opt/option

fcfbld.tlb(fcf_config) : fcf_config.txt
 $(libr) $(librflags)/text fcfbld fcf_config

fcfbld.tlb(fcf_display) : fcf_display.txt
 $(libr) $(librflags)/text fcfbld fcf_display

fcfbld.tlb(fcf_invoc) : fcf_invoc.txt
 $(libr) $(librflags)/text fcfbld fcf_invoc

fcfbld.tlb(fcf_model) : fcf_model.txt
 $(libr) $(librflags)/text fcfbld fcf_model

fcf.obj : fcfbld.tlb(fcf_config),-
          fcfbld.tlb(fcf_invoc),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctuser.inc,-
          fcf_sky^,-
          fex_grtcoawt^,-
          fex_grttrans^,-
          fex_nyquist^,-
          fex_vibcorr^,-
          fic_sky^

fcf_apodize_and_rotate.obj : fcfbld.tlb(fcf_config),-
                             fcfbld.tlb(fcf_invoc),-
                             csdr$library:futlib.tlb(fut_params),-
                             fex_grtcoawt^,-
                             fex_grttrans^,-
                             fex_nyquist^,-
                             fex_vibcorr^

fcf_apply_model.obj : fcfbld.tlb(fcf_config),-
                      fcfbld.tlb(fcf_display),-
                      fcfbld.tlb(fcf_invoc),-
                      fcfbld.tlb(fcf_model),-
                      csdr$library:futlib.tlb(fut_params),-
                      fex_grtcoawt^,-
                      fex_grttrans^,-
                      fex_nyquist^,-
                      fex_vibcorr^

fcf_autophase_correct.obj : fcfbld.tlb(fcf_model)

fcf_calc_responsivity.obj : fcfbld.tlb(fcf_display),-
                            fcfbld.tlb(fcf_model)

fcf_calibrate_spectra.obj : fcfbld.tlb(fcf_config),-
                            fcfbld.tlb(fcf_display),-
                            fcfbld.tlb(fcf_invoc),-
                            fcfbld.tlb(fcf_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquist^,-
                            fex_vibcorr^,-
                            fcf_sky^

fcf_compute_constants.obj : fcfbld.tlb(fcf_config),-
                            fcfbld.tlb(fcf_invoc),-
                            fcfbld.tlb(fcf_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquist^,-
                            fex_vibcorr^

fcf_display_model.obj : fcfbld.tlb(fcf_config),-
                        fcfbld.tlb(fcf_invoc),-
                        fcfbld.tlb(fcf_model),-
                        csdr$library:futlib.tlb(fut_params),-
                        fex_grtcoawt^,-
                        fex_grttrans^,-
                        fex_nyquist^,-
                        fex_vibcorr^

fcf_display_spectra.obj : fcfbld.tlb(fcf_config),-
                          fcfbld.tlb(fcf_display),-
                          fcfbld.tlb(fcf_invoc),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fex_nyquist^,-
                          fex_vibcorr^

fcf_doppler_shift.obj : fcfbld.tlb(fcf_model)

fcf_initialize_report.obj : fcfbld.tlb(fcf_invoc),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fcf_open_cal_coadd.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                         csdr$library:csdrlib.tlb(cct_query_tod_catalog_record),-
                         ccm_cme_catalog_entry^

fcf_open_coadd.obj : fcfbld.tlb(fcf_config),-
                     fcfbld.tlb(fcf_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_grtcoawt^,-
                     fex_grttrans^,-
                     fex_nyquist^,-
                     fex_vibcorr^

fcf_open_sky_coadd.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                         csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                         ccm_cme_catalog_entry^

fcf_open_temp_spec.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params)

fcf_parse.obj : fcfbld.tlb(fcf_invoc),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg)

fcf_pcalib_variances.obj : fcfbld.tlb(fcf_model)

fcf_produce_spectra.obj : fcfbld.tlb(fcf_config),-
                          fcfbld.tlb(fcf_display),-
                          fcfbld.tlb(fcf_invoc),-
                          fcfbld.tlb(fcf_model),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fex_nyquist^,-
                          fex_vibcorr^,-
                          fcf_sky^,-
                          fic_sky^

fcf_read_cal_coadd.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fic_sky^

fcf_read_coadd.obj : fcfbld.tlb(fcf_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:ctuser.inc,-
                     fic_sky^

fcf_read_dvector.obj : fcfbld.tlb(fcf_invoc),-
                       fcfbld.tlb(fcf_model),-
                       csdr$library:futlib.tlb(fut_params),-
                       fex_var^

fcf_read_model.obj : fcfbld.tlb(fcf_invoc),-
                     fcfbld.tlb(fcf_model),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_mod^

fcf_read_reference.obj : fcfbld.tlb(fcf_config),-
                         fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         fex_grtcoawt^,-
                         fex_grttrans^,-
                         fex_nyquist^,-
                         fex_vibcorr^

fcf_read_sky_coadd.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                         csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                         fic_sky^

fcf_read_sky_coadd_list.obj : fcfbld.tlb(fcf_invoc),-
                              csdr$library:futlib.tlb(fut_params),-
                              csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                              csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                              fic_sky^

fcf_temporal_drift.obj : fcfbld.tlb(fcf_model)

fcf_update_report.obj : fcfbld.tlb(fcf_invoc),-
                        fcfbld.tlb(fcf_model),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params)

fcf_write_cal_spec.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fcf_sky^

fcf_write_sky_spec.obj : fcfbld.tlb(fcf_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         fcf_sky^

fcf_write_temp_spec.obj : fcfbld.tlb(fcf_invoc),-
                          csdr$library:futlib.tlb(fut_params),-
                          fcf_sky^
