!
!   Purpose: Build facility fms_merge_skymaps.
!
!   Author: S. Alexander, HSTX, 4/20/94, SER 11418
!
!   Modifications: L. Rosen, HSTX, 12/20/94, SER 11977, New merge capabilities.

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fmsbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fmsbld.tlb")'" .eqs. "" then $(libr)/create/text fmsbld
 if "''f$search("fmsbld.olb")'" .eqs. "" then $(libr)/create fmsbld

fmsbld_obj = fms,-
             fms_close,-
             fms_combine_c_d_vectors,-
             fms_combine_other,-
             fms_compute_weights,-
             fms_correction_spectrum,-
             fms_init_report,-
             fms_open_skymaps,-
             fms_parse,-
             fms_read_c_d_vectors,-
             fms_read_skymaps,-
             fms_spectra_variances,-
             fms_tracking_copy,-
             fms_tracking_info,-
             fms_verify_merge,-
             fms_write_skymap,-
             fms_msg

fms : fms.exe
 ! Facility fms has been built.

fms.exe : fmsbld.tlb(fms_msg),-
          fmsbld.olb($(fmsbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fmsbld/lib/inc=fms,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option

fmsbld.tlb(fms_msg) : fms_msg.txt
 $(libr) $(librflags)/text fmsbld fms_msg

fms.obj : fmsbld.tlb(fms_msg),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctparams.inc,-
          fex_cvs^,-
          fex_var^,-
          fms_sky^

fms_close.obj : fmsbld.tlb(fms_msg),-
                csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params)

fms_combine_c_d_vectors : fmsbld.tlb(fms_msg),-
                          csdr$library:futlib.tlb(fut_params),-
                          csdr$library:ctuser.inc,-
                          fex_cvs^,-
                          fex_var^,-
                          fms_sky^

fms_combine_other : csdr$library:futlib.tlb(fut_fcs_include),-
                    csdr$library:futlib.tlb(fut_params),-
                    fms_sky^

fms_compute_weights : fmsbld.tlb(fms_msg),-
                      csdr$library:futlib.tlb(fut_params),-
                      csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                      csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                      fex_cvs^,-
                      fex_var^,-
                      fms_sky^

fms_correction_spectrum : fmsbld.tlb(fms_msg),-
                          csdr$library:futlib.tlb(fut_params),-
                          csdr$library:csdrlib.tlb(cct_get_config),-
                          csdr$library:ctparams.inc,-
                          fex_mcs^

fms_init_report : fmsbld.tlb(fms_msg),-
                  csdr$library:futlib.tlb(fut_error),-
                  csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                  ccm_cme_catalog_entry^

fms_open_skymaps : fmsbld.tlb(fms_msg),-
                   csdr$library:futlib.tlb(fut_params)

fms_parse : fmsbld.tlb(fms_msg),-
            csdr$library:csdrlib.tlb(upm_stat_msg)

fms_read_c_d_vectors : fmsbld.tlb(fms_msg),-
                       csdr$library:csdrlib.tlb(cct_get_config),-
                       csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                       csdr$library:ctparams.inc,-
                       fex_cvs^,-
                       fex_var^,-
                       ccm_cme_catalog_entry^

fms_read_skymaps : fmsbld.tlb(fms_msg),-
                   csdr$library:futlib.tlb(fut_params),-
                   csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                   csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                   fms_sky^

fms_spectra_variances : fms_sky^

fms_tracking_copy : fms_sky^

fms_tracking_info : fms_sky^

fms_verify_merge : fmsbld.tlb(fms_msg)

fms_write_skymap : fmsbld.tlb(fms_msg),-
                   fms_sky^
