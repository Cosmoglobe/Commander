!
!   Purpose: Build facility fip_initial_product
!
!   Author: Gene Eplee, GSC, 7/2/93, SER 11058
!
!   Modifications: Gene Eplee, GSC, 11/1/93, SER 11414.
!                  Added new executable FIP_LINES.
!                  Larry Rosen, HSTX, 12/6/93, SER 11259.
!                  Added new executable FIP_COVAR.
!		Gene Eplee, GSC, 03/11/94
!		Converted FIP_CONFIG_LINES.TXT to FIP_CONIFG_FREQ.TXT
!		Added new executable FIP_EJGN.
!		   Gene Eplee, GSC, 05/09/94
!			Removed FIP_REFORMAT_GN.FOR.
!			Added new executable FIP_ERR.  SER 11704.
!              Replace FIP_EJGN executable with FIP_EJV executable.
!                  Gene Eplee, GSC, 07/08/94  SER 11618
!		Added new executable FIP_FEF.
!		   Gene Eplee, GSC, 09/15/94
!		Added new executable FIP_DUST.
!		   Gene Eplee, GSC, 10/05/94, SER 11936
!		Added new executable FIP_CVS.
!		   Gene Eplee, GSC, 10/05/94, SER 11936
!               Added new executable FIP_SPECTRA_COADD (FIPA verb)
!               that processes FCF_SKY, FCF_DSK, FCF_CAL, FCF_DCL, FIC_SKY,
!               and FIC_CAL.
!                  Larry P. Rosen, HSTX, 12/15,16/94, SER 11936

.suffixes
.suffixes .exe .olb .obj .for .msg

fip : fip_sky.exe, fip_model.exe, fip_lines.exe, fip_covar.exe, fip_ejv.exe, fip_err.exe, fip_fef.exe, fip_dust.exe, fip_cvs.exe, fip_spectra_coadd.exe
    !Facility fip is up to date.

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fipbld.tlb/lib + -
                                         csdr$library:futlib.tlb/lib + -
                                         csdr$library:csdrlib.tlb/lib
.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fipbld.tlb")'" .eqs. "" then $(libr)/create/text fipbld
 if "''f$search("fipbld.olb")'" .eqs. "" then $(libr)/create fipbld


! message file

fip_msg.obj : fip_msg.msg
  message fip_msg

!  text library

fipbld.tlb(fip_invoc_sky) : fip_invoc_sky.txt
 $(libr) $(librflags)/text fipbld fip_invoc_sky

fipbld.tlb(fip_frequency) : fip_frequency.txt
 $(libr) $(librflags)/text fipbld fip_frequency

fipbld.tlb(fip_config_model) : fip_config_model.txt
 $(libr) $(librflags)/text fipbld fip_config_model

fipbld.tlb(fip_invoc_model) : fip_invoc_model.txt
 $(libr) $(librflags)/text fipbld fip_invoc_model

fipbld.tlb(fip_model) : fip_model.txt
 $(libr) $(librflags)/text fipbld fip_model

fipbld.tlb(fip_config_freq) : fip_config_freq.txt
 $(libr) $(librflags)/text fipbld fip_config_freq

fipbld.tlb(fip_invoc_lines) : fip_invoc_lines.txt
 $(libr) $(librflags)/text fipbld fip_invoc_lines

fipbld.tlb(fip_invoc_ejv) : fip_invoc_ejv.txt
 $(libr) $(librflags)/text fipbld fip_invoc_ejv

fipbld.tlb(fip_invoc_err) : fip_invoc_err.txt
 $(libr) $(librflags)/text fipbld fip_invoc_err

fipbld.tlb(fip_invoc_fef) : fip_invoc_fef.txt
 $(libr) $(librflags)/text fipbld fip_invoc_fef

fipbld.tlb(fip_invoc_dust) : fip_invoc_dust.txt
 $(libr) $(librflags)/text fipbld fip_invoc_dust

fipbld.tlb(fip_invoc_cvs) : fip_invoc_cvs.txt
 $(libr) $(librflags)/text fipbld fip_invoc_cvs
  ! fipbld.tlb has been built


! Build fip_sky

fip_sky_txt = fip_invoc_sky,-
              fip_frequency

fip_sky_obj = fip_sky,-
              fip_close_skymaps,-
              fip_frequency_cut,-
              fip_galactic_cut,-
              fip_init_sky_report,-
              fip_open_skymaps,-
              fip_parse_sky,-
              fip_reformat_spectrum,-
              fip_update_sky_report,-
              fip_msg

fip_sky.exe : fipbld.tlb($(fip_sky_txt)),-
              fipbld.olb($(fip_sky_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_sky,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_sky has been built.


! Build fip_model

fip_model_txt = fip_config_model,-
                fip_frequency,-
                fip_invoc_model,-
                fip_model

fip_model_obj = fip_model,-
                fip_calc_responsivity,-
                fip_frequency_cut,-
                fip_init_model_report,-
                fip_parse_model,-
                fip_read_model,-
                fip_read_reference,-
                fip_update_model_report,-
                fip_write_model,-
                fip_msg

fip_model.exe : fipbld.tlb($(fip_model_txt)),-
                fipbld.olb($(fip_model_obj)),-
                csdr$library:futlib.olb,-
                csdr$library:csdrlib.olb,-
                csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_model,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_model has been built.


! Build fip_lines

fip_lines_txt = fip_config_freq,-
                fip_frequency,-
                fip_invoc_lines

fip_lines_obj = fip_lines,-
                fip_frequency_cut,-
                fip_init_lines_report,-
                fip_parse_lines,-
                fip_read_nyquist,-
                fip_reformat_hilines,-
                fip_reformat_lolines,-
                fip_update_lines_report,-
                fip_msg

fip_lines.exe : fipbld.tlb($(fip_lines_txt)),-
                fipbld.olb($(fip_lines_obj)),-
                csdr$library:futlib.olb,-
                csdr$library:csdrlib.olb,-
                csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_lines,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_lines has been built.


! Build fip_covar

fip_covar_txt = fip_config_freq,-
                fip_frequency

fip_covar_obj = fip_covar,-
                fip_parse_cov,-
                fip_init_report_cov,-
                fip_read_cov,-
                fip_convert_cov,-
                fip_read_nyquist,-
                fip_frequency_cut,-
                fip_write_cov,-
                fip_transcribe_cov,-
                fip_pack_cov,-
                fip_msg

fip_covar.exe : fipbld.tlb($(fip_covar_txt)),-
                fipbld.olb($(fip_covar_obj)),-
                csdr$library:futlib.olb,-
                csdr$library:csdrlib.olb,-
                csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_covar,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_covar has been built.


! Build fip_ejv

fip_ejv_txt = fip_config_freq,-
              fip_frequency,-
              fip_invoc_ejv

fip_ejv_obj = fip_ejv,-
              fip_frequency_cut,-
              fip_init_ejv_report,-
              fip_parse_ejv,-
              fip_read_nyquist,-
              fip_reformat_ejv,-
              fip_update_ejv_report,-
              fip_msg

fip_ejv.exe : fipbld.tlb($(fip_ejv_txt)),-
              fipbld.olb($(fip_ejv_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_ejv,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_ejv has been built.


! Build fip_err

fip_err_txt = fip_config_freq,-
              fip_frequency,-
              fip_invoc_err

fip_err_obj = fip_err,-
              fip_frequency_cut,-
              fip_init_err_report,-
              fip_parse_err,-
              fip_read_nyquist,-
              fip_reformat_err,-
              fip_update_err_report,-
              fip_msg

fip_err.exe : fipbld.tlb($(fip_err_txt)),-
              fipbld.olb($(fip_err_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_err,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_err has been built.


! Build fip_fef

fip_fef_txt = fip_config_freq,-
              fip_frequency,-
              fip_invoc_fef

fip_fef_obj = fip_fef,-
              fip_frequency_cut,-
              fip_init_fef_report,-
              fip_parse_fef,-
              fip_read_fef,-
              fip_read_nyquist,-
              fip_reformat_fef,-
              fip_update_fef_report,-
              fip_write_fef,-
              fip_msg

fip_fef.exe : fipbld.tlb($(fip_fef_txt)),-
              fipbld.olb($(fip_fef_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_fef,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_fef has been built.


! Build fip_dust

fip_dust_txt = fip_config_freq,-
               fip_frequency,-
               fip_invoc_dust

fip_dust_obj = fip_dust,-
               fip_frequency_cut,-
               fip_init_dust_report,-
               fip_parse_dust,-
               fip_read_nyquist,-
               fip_reformat_dust,-
               fip_update_dust_report,-
               fip_msg

fip_dust.exe : fipbld.tlb($(fip_dust_txt)),-
               fipbld.olb($(fip_dust_obj)),-
               csdr$library:futlib.olb,-
               csdr$library:csdrlib.olb,-
               csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_dust,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_dust has been built.


! Build fip_cvs

fip_cvs_txt = fip_config_freq,-
              fip_frequency,-
              fip_invoc_cvs

fip_cvs_obj = fip_cvs,-
              fip_frequency_cut,-
              fip_init_cvs_report,-
              fip_parse_cvs,-
              fip_read_nyquist,-
              fip_reformat_cvs,-
              fip_update_cvs_report,-
              fip_msg

fip_cvs.exe : fipbld.tlb($(fip_cvs_txt)),-
              fipbld.olb($(fip_cvs_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_cvs,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_cvs has been built.


fip_spectra_coadd_txt = fip_invoc_sky,-
                        fip_frequency

fip_spectra_coadd_obj = fip_spectra_coadd,-
                        fip_frequency_cut,-
                        fip_galactic_cut,-
                        fip_reformat_ifg, -
                        fip_reformat_spectrum,-
                        fip_sc_close_cal, -
                        fip_sc_close_sky, -
                        fip_sc_init_report, -
                        fip_sc_open_sky, -
                        fip_sc_open_cal, -
                        fip_sc_parse, -
                        fip_sc_process_cal, -
                        fip_sc_process_sky, -
                        fip_sc_update_report, -
                        fip_msg

fip_spectra_coadd.exe : fipbld.tlb($(fip_spectra_coadd_txt)),-
              fipbld.olb($(fip_spectra_coadd_obj)),-
              csdr$library:futlib.olb,-
              csdr$library:csdrlib.olb,-
              csdr$library:v5.opt
 $(link) $(linkflags) fipbld/lib/inc=fip_spectra_coadd,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fip_spectra_coadd has been built.

! Dependencies

fip_calc_responsivity.obj : fipbld.tlb(fip_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_mod^,-
                            fip_mod^

fip_close_skymaps.obj : fipbld.tlb(fip_invoc_sky),-
                        csdr$library:futlib.tlb(fut_params)

fip_convert_cov.obj : fipbld.tlb(fip_config_freq),-
                      fipbld.tlb(fip_frequency),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                      ccm_cme_catalog_entry^,-
                      fcc_cov^,-
                      fex_nyquist^

fip_covar.obj : csdr$library:futlib.tlb(fut_error),-
                csdr$library:ctparams.inc,-
                fcc_cov^

fip_cvs.obj : fipbld.tlb(fip_invoc_cvs),-
              fipbld.tlb(fip_config_freq),-
              csdr$library:futlib.tlb(fut_error),-
              csdr$library:futlib.tlb(fut_params),-
              fex_nyquist^

fip_dust.obj : fipbld.tlb(fip_invoc_dust),-
               fipbld.tlb(fip_config_freq),-
               csdr$library:futlib.tlb(fut_error),-
               csdr$library:futlib.tlb(fut_params),-
               fla_dst^,-
               fip_dst^,-
               fex_nyquist^

fip_ejv.obj : fipbld.tlb(fip_invoc_ejv),-
              fipbld.tlb(fip_config_freq),-
              csdr$library:futlib.tlb(fut_error),-
              csdr$library:futlib.tlb(fut_params),-
              fex_nyquist^

fip_err.obj : fipbld.tlb(fip_invoc_err),-
              fipbld.tlb(fip_config_freq),-
              csdr$library:futlib.tlb(fut_error),-
              csdr$library:futlib.tlb(fut_params),-
              fex_nyquist^

fip_fef.obj : fipbld.tlb(fip_invoc_fef),-
              fipbld.tlb(fip_config_freq),-
              csdr$library:futlib.tlb(fut_error),-
              csdr$library:futlib.tlb(fut_params),-
              fef_spc^,-
              fex_nyquist^,-
              fip_fef^

fip_frequency_cut.obj : fipbld.tlb(fip_frequency),-
                        csdr$library:futlib.tlb(fut_params)

fip_galactic_cut.obj : fipbld.tlb(fip_invoc_sky),-
                       csdr$library:futlib.tlb(fut_params)

fip_init_cvs_report.obj : fipbld.tlb(fip_invoc_cvs),-
                          csdr$library:futlib.tlb(fut_error),-
                          csdr$library:futlib.tlb(fut_params)

fip_init_dust_report.obj : fipbld.tlb(fip_invoc_dust),-
                           csdr$library:futlib.tlb(fut_error),-
                           csdr$library:futlib.tlb(fut_params)

fip_init_ejv_report.obj : fipbld.tlb(fip_invoc_ejv),-
                          csdr$library:futlib.tlb(fut_error),-
                          csdr$library:futlib.tlb(fut_params)

fip_init_err_report.obj : fipbld.tlb(fip_invoc_err),-
                          csdr$library:futlib.tlb(fut_error),-
                          csdr$library:futlib.tlb(fut_params)

fip_init_fef_report.obj : fipbld.tlb(fip_invoc_fef),-
                          csdr$library:futlib.tlb(fut_error),-
                          csdr$library:futlib.tlb(fut_params)

fip_init_lines_report.obj : fipbld.tlb(fip_invoc_lines),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_init_model_report.obj : fipbld.tlb(fip_invoc_model),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_init_report_cov.obj : csdr$library:futlib.tlb(fut_error)

fip_init_sky_report.obj : fipbld.tlb(fip_invoc_sky),-
                          csdr$library:futlib.tlb(fut_error),-
                          csdr$library:futlib.tlb(fut_params)

fip_lines.obj : fipbld.tlb(fip_invoc_lines),-
                fipbld.tlb(fip_config_freq),-
                csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params),-
                fex_nyquist^

fip_model.obj : fipbld.tlb(fip_invoc_model),-
                fipbld.tlb(fip_config_model),-
                csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params),-
                fex_nyquist^

fip_open_skymaps.obj : fipbld.tlb(fip_invoc_sky),-
                       csdr$library:futlib.tlb(fut_params)

fip_open_cal.obj : fipbld.tlb(fip_invoc_sky),-
                   csdr$library:futlib.tlb(fut_params),-
                   csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                   ccm_cme_catalog_entry^

fip_parse_cov.obj : csdr$library:csdrlib.tlb(upm_stat_msg)

fip_parse_cvs.obj : fipbld.tlb(fip_invoc_cvs),-
                    fipbld.tlb(fip_config_freq),-
                    fipbld.tlb(fip_frequency),-
                    csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(upm_stat_msg),-
                    fex_nyquist^

fip_parse_dust.obj : fipbld.tlb(fip_invoc_dust),-
                     fipbld.tlb(fip_config_freq),-
                     fipbld.tlb(fip_frequency),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:csdrlib.tlb(upm_stat_msg),-
                     fex_nyquist^

fip_parse_ejv.obj : fipbld.tlb(fip_invoc_ejv),-
                    fipbld.tlb(fip_config_freq),-
                    fipbld.tlb(fip_frequency),-
                    csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(upm_stat_msg),-
                    fex_nyquist^

fip_parse_err.obj : fipbld.tlb(fip_invoc_err),-
                    fipbld.tlb(fip_config_freq),-
                    fipbld.tlb(fip_frequency),-
                    csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(upm_stat_msg),-
                    fex_nyquist^

fip_parse_fef.obj : fipbld.tlb(fip_invoc_fef),-
                    fipbld.tlb(fip_config_freq),-
                    fipbld.tlb(fip_frequency),-
                    csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(upm_stat_msg),-
                    fex_nyquist^

fip_parse_lines.obj : fipbld.tlb(fip_invoc_lines),-
                      fipbld.tlb(fip_config_freq),-
                      fipbld.tlb(fip_frequency),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(upm_stat_msg),-
                      fex_nyquist^

fip_parse_model.obj : fipbld.tlb(fip_invoc_model),-
                      fipbld.tlb(fip_config_model),-
                      fipbld.tlb(fip_frequency),-
                      fipbld.tlb(fip_model),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(upm_stat_msg),-
                      fex_mod^,-
                      fip_mod^,-
                      fex_nyquist^

fip_parse_sky.obj : fipbld.tlb(fip_invoc_sky),-
                    fipbld.tlb(fip_frequency),-
                    csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(upm_stat_msg)

fip_read_cov.obj : csdr$library:ctparams.inc,-
                   fcc_cov^

fip_read_fef.obj : fipbld.tlb(fip_invoc_fef),-
                   csdr$library:futlib.tlb(fut_params),-
                   fef_spc^

fip_read_model.obj : fipbld.tlb(fip_invoc_model),-
                     fipbld.tlb(fip_model),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_mod^,-
                     fip_mod^

fip_read_nyquist.obj : fipbld.tlb(fip_config_freq),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:ctuser.inc,-
                       fex_nyquist^

fip_read_reference.obj : fipbld.tlb(fip_invoc_model),-
                         fipbld.tlb(fip_config_model),-
                         csdr$library:ctuser.inc,-
                         csdr$library:futlib.tlb(fut_params),-
                         fex_nyquist^

fip_reformat_cvs.obj : fipbld.tlb(fip_invoc_cvs),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       fms_cvs^,-
                       fip_cvs^

fip_reformat_dust.obj : fipbld.tlb(fip_invoc_dust),-
                        fipbld.tlb(fip_frequency),-
                        csdr$library:futlib.tlb(fut_params),-
                        fla_dst^,-
                        fip_dst^

fip_reformat_ejv.obj : fipbld.tlb(fip_invoc_ejv),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       csdr$library:csdrlib.tlb(uoe_constants),-
                       fex_ejv^,-
                       fip_ejv^

fip_reformat_err.obj : fipbld.tlb(fip_invoc_err),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       fer_err^,-
                       fip_err^

fip_reformat_err.obj : fipbld.tlb(fip_invoc_err),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       fer_err^,-
                       fip_err^

fip_reformat_fef.obj : fipbld.tlb(fip_invoc_fef),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       csdr$library:csdrlib.tlb(uoe_constants),-
                       fef_spc^,-
                       fip_fef^

fip_reformat_hilines.obj : fipbld.tlb(fip_invoc_lines),-
                           fipbld.tlb(fip_frequency),-
                           csdr$library:futlib.tlb(fut_params),-
                           fip_hlp^

fip_reformat_ifg.obj : fipbld.tlb(fip_invoc_sky),-
                       fipbld.tlb(fip_frequency),-
                       csdr$library:futlib.tlb(fut_params),-
                       csdr$library:csdrlib.tlb(uoe_constants),-
                       fic_sky^,-
                       fip_isk^

fip_reformat_lolines.obj : fipbld.tlb(fip_invoc_lines),-
                           fipbld.tlb(fip_frequency),-
                           csdr$library:futlib.tlb(fut_params),-
                           fip_llp^

fip_reformat_spectrum.obj : fipbld.tlb(fip_invoc_sky),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_params),-
                            csdr$library:csdrlib.tlb(uoe_constants),-
                            fcf_sky^,-
                            fip_sky^

fip_sc_close_cal.obj : fipbld.tlb(fip_invoc_sky),-
                        csdr$library:futlib.tlb(fut_params),-
                        csdr$library:ctuser.inc

fip_sc_close_sky.obj : fipbld.tlb(fip_invoc_sky),-
                        csdr$library:futlib.tlb(fut_params)

fip_sc_init_report.obj : fipbld.tlb(fip_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params)

fip_sc_open_sky.obj : fipbld.tlb(fip_invoc_sky),-
                      csdr$library:futlib.tlb(fut_params)

fip_sc_open_cal.obj : fipbld.tlb(fip_invoc_sky),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                      ccm_cme_catalog_entry^

fip_sc_parse.obj : fipbld.tlb(fip_invoc_sky),-
                   fipbld.tlb(fip_frequency),-
                   csdr$library:futlib.tlb(fut_params),-
                   csdr$library:csdrlib.tlb(upm_stat_msg)

fip_sc_process_cal.obj : fipbld.tlb(fip_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:ctuser.inc,-
                         fcf_cal^,-
                         fic_cal^,-
                         fip_ccl^,-
                         fip_icl^

fip_sc_process_sky.obj : fipbld.tlb(fip_invoc_sky),-
                         csdr$library:futlib.tlb(fut_error),-
                         csdr$library:futlib.tlb(fut_params),-
                         csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                         csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                         csdr$library:ctuser.inc,-
                         fcf_sky^,-
                         fic_sky^,-
                         fip_csk^,-
                         fip_isk^

fip_sc_update_report.obj : fipbld.tlb(fip_invoc_cvs),-
                           fipbld.tlb(fip_frequency),-
                           csdr$library:futlib.tlb(fut_error),-
                           csdr$library:futlib.tlb(fut_params)

fip_sky.obj : fipbld.tlb(fip_invoc_sky),-
              csdr$library:futlib.tlb(fut_error),-
              csdr$library:futlib.tlb(fut_params),-
              csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
              csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
              csdr$library:ctuser.inc,-
              fcf_sky^,-
              fip_sky^

fip_spectra_coadd.obj : fipbld.tlb(fip_invoc_sky),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params),-
                        csdr$library:ctuser.inc

fip_update_cvs_report.obj : fipbld.tlb(fip_invoc_cvs),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_update_dust_report.obj : fipbld.tlb(fip_invoc_dust),-
                             fipbld.tlb(fip_frequency),-
                             csdr$library:futlib.tlb(fut_error),-
                             csdr$library:futlib.tlb(fut_params)

fip_update_ejv_report.obj : fipbld.tlb(fip_invoc_ejv),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_update_err_report.obj : fipbld.tlb(fip_invoc_err),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_update_fef_report.obj : fipbld.tlb(fip_invoc_fef),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_update_lines_report.obj : fipbld.tlb(fip_invoc_lines),-
                              fipbld.tlb(fip_frequency),-
                              csdr$library:futlib.tlb(fut_error),-
                              csdr$library:futlib.tlb(fut_params)

fip_update_model_report.obj : fipbld.tlb(fip_invoc_model),-
                              fipbld.tlb(fip_frequency),-
                              csdr$library:futlib.tlb(fut_error),-
                              csdr$library:futlib.tlb(fut_params)

fip_update_sky_report.obj : fipbld.tlb(fip_invoc_sky),-
                            fipbld.tlb(fip_frequency),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

fip_write_cov.obj : csdr$library:ctparams.inc,-
                    fcc_cov^,-
                    fip_cov^

fip_write_fef.obj : fipbld.tlb(fip_invoc_fef),-
                    csdr$library:futlib.tlb(fut_params),-
                    fip_fef^

fip_write_model.obj : fipbld.tlb(fip_invoc_model),-
                      fipbld.tlb(fip_config_model),-
                      fipbld.tlb(fip_frequency),-
                      fipbld.tlb(fip_model),-
                      csdr$library:futlib.tlb(fut_params),-
                      fex_nyquist^,-
                      fex_mod^,-
                      fip_mod^
