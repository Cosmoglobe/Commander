! DESCRIP.MMS
!
!   Purpose: Build facility FFL_Fishinput_Long
!
!   Author: Fred Shuman, HSTX, 1995 Jun 30, SPR 12272
!-------------------------------------------------------------------------------
!   Modifications:
!	This file is a descendent of DESCRIP.MMS for facility FFI_Fish_Input,
!	Gene Eplee, GSC, 1993 Mar 09, SER 10763
!
!	Subsumed the FFL_SPEC.TXT include file into FFL_CONFIG.TXT .
!	Fred Shuman, HSTX, 1995 July 21.
!-------------------------------------------------------------------------------

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fflbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fflbld.tlb")'" .eqs. "" then $(libr)/create/text fflbld
 if "''f$search("fflbld.olb")'" .eqs. "" then $(libr)/create fflbld

fflbld_txt = ffl_config,-
             ffl_invoc

fflbld_obj = ffl,-
             ffl_parse,-
             ffl_initialize_report,-
             ffl_get_reference,-
             ffl_open_coadd,-
             ffl_read_hybrid_coadd,-
             ffl_bracket_coadd,-
             ffl_read_coadd,-
             ffl_produce_spectra,-
             ffl_update_report,-
             ffl_msg

ffl : ffl.exe
 ! Facility FFL has been built.

ffl.exe : fflbld.tlb($(fflbld_txt)),-
          fflbld.olb($(fflbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fflbld/lib/inc=ffl,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      csdr$library:v5.opt/option

fflbld.tlb(ffl_config) : ffl_config.txt
 $(libr) $(librflags)/text fflbld ffl_config

fflbld.tlb(ffl_invoc) : ffl_invoc.txt
 $(libr) $(librflags)/text fflbld ffl_invoc

ffl.obj : fflbld.tlb(ffl_config),-
          fflbld.tlb(ffl_invoc),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctuser.inc,-
          fex_grtcoawt^,-
          fex_grttrans^,-
          fil_sky^

ffl_parse.obj : fflbld.tlb(ffl_invoc),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg)

ffl_initialize_report.obj : fflbld.tlb(ffl_invoc),-
                            csdr$library:futlib.tlb(fut_error),-
                            csdr$library:futlib.tlb(fut_params)

ffl_get_reference.obj : fflbld.tlb(ffl_config),-
                        fflbld.tlb(ffl_invoc),-
                        csdr$library:futlib.tlb(fut_params),-
                        fex_grtcoawt^,-
                        fex_grttrans^

ffl_open_coadd.obj : fflbld.tlb(ffl_invoc),-
                     fflbld.tlb(ffl_config),-
                     csdr$library:futlib.tlb(fut_params)

ffl_read_hybrid_coadd.obj : fflbld.tlb(ffl_invoc),-
                            csdr$library:futlib.tlb(fut_params),-
                            csdr$library:ctuser.inc,-
                            fex_grtcoawt^,-
                            fex_grttrans^,-
                            fil_sky^

ffl_bracket_coadd.obj : fflbld.tlb(ffl_invoc),-
                        csdr$library:ctuser.inc,-
                        fil_sky^

ffl_read_coadd.obj : fflbld.tlb(ffl_invoc),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:ctuser.inc,-
                     fil_sky^

ffl_produce_spectra.obj : fflbld.tlb(ffl_config),-
                          fflbld.tlb(ffl_invoc),-
                          csdr$library:futlib.tlb(fut_params),-
                          fex_grtcoawt^,-
                          fex_grttrans^,-
                          fil_sky^

ffl_update_report.obj : fflbld.tlb(ffl_invoc),-
                        csdr$library:futlib.tlb(fut_error),-
                        csdr$library:futlib.tlb(fut_params)
