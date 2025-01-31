!
!   Purpose: Build facility fad_apply_destriper
!
!   Author: S. Alexander, STX, 5/5/93, SER 10798
!
!   Modifications: L. Rosen, HSTX, 4/13/94, SER 11688, New fex_gn records.
!                  L. Rosen, HSTX, 6/20/94, SER 11796, Change fex_ej to fex_ejv.

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fadbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fadbld.tlb")'" .eqs. "" then $(libr)/create/text fadbld
 if "''f$search("fadbld.olb")'" .eqs. "" then $(libr)/create fadbld

fadbld_txt = fad_msg

fadbld_obj = fad,-
             fad_apply_correction,-
             fad_close_archives,-
             fad_close_report,-
             fad_init_report,-
             fad_open_archives,-
             fad_parse_command,-
             fad_process_spectra,-
             fad_read_reference_data,-
             fad_msg

fad : fad.exe
 ! Facility fad has been built.

fad.exe : fadbld.tlb($(fadbld_txt)),-
          fadbld.olb($(fadbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fadbld/lib/inc=fad,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option

fadbld.tlb(fad_msg) : fad_msg.txt
 $(libr) $(librflags)/text fadbld fad_msg

fad.obj : fadbld.tlb(fad_msg),-
          csdr$library:futlib.tlb(fut_error),-
          fex_ejv^

fad_apply_correction.obj : fadbld.tlb(fad_msg),-
                           csdr$library:futlib.tlb(fut_params),-
                           fex_ejv^

fad_close_archives.obj : fadbld.tlb(fad_msg),-
                         csdr$library:futlib.tlb(fut_params)

fad_close_report.obj : fadbld.tlb(fad_msg)

fad_init_report.obj : fadbld.tlb(fad_msg),-
                      csdr$library:futlib.tlb(fut_error)

fad_open_archives.obj : fadbld.tlb(fad_msg),-
                        csdr$library:futlib.tlb(fut_params),-
                        csdr$library:ctparams.inc

fad_parse_command.obj : fadbld.tlb(fad_msg),-
                        csdr$library:csdrlib.tlb(upm_stat_msg)

fad_process_spectra.obj : fadbld.tlb(fad_msg),-
                          csdr$library:futlib.tlb(fut_params),-
                          csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
                          csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
                          fex_ejv^,-
                          fad_sky^,-
                          fcf_sky^

fad_read_reference_data.obj : fadbld.tlb(fad_msg),-
                              fex_ejv^
