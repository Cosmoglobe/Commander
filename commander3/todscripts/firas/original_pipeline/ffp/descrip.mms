!
!   Purpose: Build facility ffp_final_product
!
!   Author: S. Brodd, HSTX, 3/21/96, SPR 12316
!
!   Modifications:
!

.suffixes
.suffixes .exe .olb .obj .for .msg

ffp : ffp_spectra_coadd.exe
 ! Facility ffp is up to date.

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         ffpbld.tlb/lib + -
                                         csdr$library:futlib.tlb/lib + -
                                         csdr$library:csdrlib.tlb/lib
.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("ffpbld.tlb")'" .eqs. "" then $(libr)/create/text ffpbld
 if "''f$search("ffpbld.olb")'" .eqs. "" then $(libr)/create ffpbld

ffp_msg.obj : ffp_msg.msg
 message ffp_msg

ffpbld.tlb(ffp_invoc_sky) : ffp_invoc_sky.txt
 $(libr) $(librflags)/text ffpbld ffp_invoc_sky

ffpbld.tlb(ffp_frequency) : ffp_frequency.txt
 $(libr) $(librflags)/text ffpbld ffp_frequency
 ! Text library ffpbld.tlb has been built.

ffp_spectra_coadd_txt = ffp_invoc_sky,-
                        ffp_frequency

ffp_spectra_coadd_obj = ffp_spectra_coadd,-
                        ffp_frequency_cut,-
                        ffp_galactic_cut,-
                        ffp_reformat_ifg, -
                        ffp_reformat_spectrum,-
                        ffp_sc_close_cal, -
                        ffp_sc_close_sky, -
                        ffp_sc_init_report, -
                        ffp_sc_open_sky, -
                        ffp_sc_open_cal, -
                        ffp_sc_parse, -
                        ffp_sc_process_cal, -
                        ffp_sc_process_sky, -
                        ffp_sc_update_report, -
                        ffp_msg

ffp_spectra_coadd.exe : ffpbld.tlb($(ffp_spectra_coadd_txt)),-
                        ffpbld.olb($(ffp_spectra_coadd_obj)),-
                        csdr$library:futlib.olb,-
                        csdr$library:csdrlib.olb,-
                        csdr$library:v5.opt
 $(link) $(linkflags) ffpbld/lib/inc=ffp_spectra_coadd,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility ffp_spectra_coadd has been built.

ffp_frequency_cut.obj : ffpbld.tlb(ffp_frequency),-
                        csdr$library:futlib.tlb(fut_params)

ffp_galactic_cut.obj : ffpbld.tlb(ffp_invoc_sky),-
                       csdr$library:futlib.tlb(fut_params)

ffp_reformat_ifg.obj : ffpbld.tlb(ffp_invoc_sky),-
                       ffpbld.tlb(ffp_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       csdr$library:csdrlib.tlb(uoe_constants),-
                       fil_sky^,-
                       ffp_ifg^

ffp_reformat_spectrum.obj : ffpbld.tlb(ffp_invoc_sky),-
                            ffpbld.tlb(ffp_frequency),-
                            csdr$library:futlib.tlb(fut_params),-
                            csdr$library:csdrlib.tlb(uoe_constants),-
                            fsl_sky^,-
                            ffp_spc^

ffp_sc_close_cal.obj : ffpbld.tlb(ffp_invoc_sky),-
                       csdr$library:futlib.tlb(fut_params),-
                       csdr$library:ctuser.inc

ffp_sc_close_sky.obj : ffpbld.tlb(ffp_invoc_sky),-
                       csdr$library:futlib.tlb(fut_params)

ffp_sc_init_report.obj : ffpbld.tlb(ffp_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params)

ffp_sc_open_sky.obj : ffpbld.tlb(ffp_invoc_sky),-
                      csdr$library:futlib.tlb(fut_params)

ffp_sc_open_cal.obj : ffpbld.tlb(ffp_invoc_sky),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                      ccm_cme_catalog_entry^

ffp_sc_parse.obj : ffpbld.tlb(ffp_invoc_sky),-
                   ffpbld.tlb(ffp_frequency),-
                   csdr$library:futlib.tlb(fut_params),-
                   csdr$library:csdrlib.tlb(upm_stat_msg)

ffp_sc_process_cal.obj : ffpbld.tlb(ffp_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fsl_sky^,-
                         fil_sky^,-
                         ffp_spc^,-
                         ffp_ifg^

ffp_sc_process_sky.obj : ffpbld.tlb(ffp_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                         csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                         csdr$library:ctuser.inc,-
                         fsl_sky^,-
                         fil_sky^,-
                         ffp_spc^,-
                         ffp_ifg^

ffp_sc_update_report.obj : ffpbld.tlb(ffp_invoc_sky),-
                           ffpbld.tlb(ffp_frequency),-
                           csdr$library:futlib.tlb(fut_error),-
                           csdr$library:futlib.tlb(fut_params)

ffp_spectra_coadd.obj : ffpbld.tlb(ffp_invoc_sky),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params),-
                        csdr$library:ctuser.inc
