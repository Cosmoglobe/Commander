!
!   Purpose: Build facility ffi_fish_input
!
!   Author: Gene Eplee, GSC, 3/9/93, SER 10763
!
!   Modifications:
!
!   Changes to recover low frequency short fast data.
!   Gene Eplee, GSC, 25 October 1993
!   SER 11690
!

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         ffibld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("ffibld.tlb")'" .eqs. "" then $(libr)/create/text ffibld
 if "''f$search("ffibld.olb")'" .eqs. "" then $(libr)/create ffibld

ffibld_txt = ffi_config,-
             ffi_invoc,-
             ffi_spec

ffibld_obj = ffi,-
             ffi_apodize_and_rotate,-
             ffi_bracket_coadd,-
             ffi_compute_constants,-
             ffi_get_reference,-
             ffi_initialize_report,-
             ffi_open_coadd,-
             ffi_parse,-
             ffi_produce_spectra,-
             ffi_read_coadd,-
             ffi_read_hybrid_coadd,-
             ffi_read_reference,-
             ffi_update_report,-
             ffi_write_spec,-
             ffi_msg

ffi : ffi.exe
 ! Facility ffi has been built.

ffi.exe : ffibld.tlb($(ffibld_txt)),-
          ffibld.olb($(ffibld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) ffibld/lib/inc=ffi,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      csdr$library:v5.opt/option

ffibld.tlb(ffi_config) : ffi_config.txt
 $(libr) $(librflags)/text ffibld ffi_config

ffibld.tlb(ffi_invoc) : ffi_invoc.txt
 $(libr) $(librflags)/text ffibld ffi_invoc

ffibld.tlb(ffi_spec) : ffi_spec.txt
 $(libr) $(librflags)/text ffibld ffi_spec

ffi.obj : ffibld.tlb(ffi_config),-
          ffibld.tlb(ffi_invoc),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctuser.inc,-
          fex_grtcoawt^,-
          fex_grttrans^,-
          fex_nyquist^,-
          fic_sky^

ffi_apodize_and_rotate.obj : ffibld.tlb(ffi_config),-
                             ffibld.tlb(ffi_invoc),-
                             ffibld.tlb(ffi_spec),-
                             csdr$library:futlib.tlb(fut_params),-
                             fex_grtcoawt^,-
                             fex_grttrans^,-
                             fex_nyquist^

ffi_bracket_coadd.obj : ffibld.tlb(ffi_invoc),-
                        csdr$library:ctuser.inc,-
                        fic_sky^

ffi_compute_constants.obj : ffibld.tlb(ffi_config),-
                            ffibld.tlb(ffi_invoc),-
                            ffibld.tlb(ffi_spec),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquist^

ffi_get_reference.obj : ffibld.tlb(ffi_config),-
                        ffibld.tlb(ffi_invoc),-
                        csdr$library:futlib.tlb(fut_params),-
                        fex_grtcoawt^,-
                        fex_grttrans^,-
                        fex_nyquist^

ffi_initialize_report.obj : ffibld.tlb(ffi_invoc),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

ffi_open_coadd.obj : ffibld.tlb(ffi_invoc),-
                     ffibld.tlb(ffi_spec),-
                     csdr$library:futlib.tlb(fut_params)

ffi_parse.obj : ffibld.tlb(ffi_invoc),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg)

ffi_produce_spectra.obj : ffibld.tlb(ffi_config),-
                          ffibld.tlb(ffi_invoc),-
                          ffibld.tlb(ffi_spec),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fex_nyquist^,-
                          fic_sky^

ffi_read_coadd.obj : ffibld.tlb(ffi_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:ctuser.inc,-
                     fic_sky^

ffi_read_hybrid_coadd.obj : ffibld.tlb(ffi_config),-
                            ffibld.tlb(ffi_invoc),-
                            csdr$library:futlib.tlb(fut_params),-
                            csdr$library:ctuser.inc,-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fex_nyquist^,-
                            fic_sky^

ffi_read_reference.obj : ffibld.tlb(ffi_config),-
                         ffibld.tlb(ffi_invoc),-
                         csdr$library:futlib.tlb(fut_params),-
                         fex_grtcoawt^,-
                         fex_grttrans^,-
                         fex_nyquist^

ffi_update_report.obj : ffibld.tlb(ffi_invoc),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params)

ffi_write_spec.obj : ffibld.tlb(ffi_invoc),-
                     ffibld.tlb(ffi_spec),-
                     csdr$library:futlib.tlb(fut_params)
