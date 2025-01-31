!
!   Purpose: Build facility fla_lines_analysis
!
!   Author: Gene Eplee, GSC, 7/1/93, SER 11413
!
!   Modifications:

.suffixes
.suffixes .exe .olb .obj .for .msg

fla : fla_gain.exe
    !Facility fla is up to date.

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         flabld.tlb/lib + -
                                         csdr$library:futlib.tlb/lib + -
                                         csdr$library:csdrlib.tlb/lib
.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("flabld.tlb")'" .eqs. "" then $(libr)/create/text flabld
 if "''f$search("flabld.olb")'" .eqs. "" then $(libr)/create flabld


! message file

fla_msg.obj : fla_msg.msg
  message fla_msg

!  text library

flabld.tlb(fla_invoc_gain) : fla_invoc_gain.txt
 $(libr) $(librflags)/text flabld fla_invoc_gain

flabld.tlb(fla_config_gain) : fla_config_gain.txt
 $(libr) $(librflags)/text flabld fla_config_gain

flabld.tlb(fla_model) : fla_model.txt
 $(libr) $(librflags)/text flabld fla_model
  ! flabld.tlb has been built


! Build fla_gain

fla_gain_txt = fla_config_gain,-
               fla_invoc_gain,-
               fla_model

fla_gain_obj = fla_gain,-
               fla_calc_responsivity,-
               fla_init_gain_report,-
               fla_parse_gain,-
               fla_read_model,-
               fla_read_reference,-
               fla_update_gain_report,-
               fla_write_gain,-
               fla_msg

fla_gain.exe : flabld.tlb($(fla_gain_txt)),-
               flabld.olb($(fla_gain_obj)),-
               csdr$library:futlib.olb,-
               csdr$library:csdrlib.olb,-
               csdr$library:v5.opt
 $(link) $(linkflags) flabld/lib/inc=fla_gain,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option
 ! Facility fla_gain has been built.


! Dependencies

fla_calc_responsivity.obj : flabld.tlb(fla_model),-
                            csdr$library:futlib.tlb(fut_params),-
                            fex_mod^

fla_init_gain_report.obj : flabld.tlb(fla_invoc_gain),-
                           csdr$library:futlib.tlb(fut_error),-
                           csdr$library:futlib.tlb(fut_params)

fla_gain.obj : flabld.tlb(fla_invoc_gain),-
               csdr$library:futlib.tlb(fut_error),-
               csdr$library:futlib.tlb(fut_params)

fla_parse_gain.obj : flabld.tlb(fla_invoc_gain),-
                     flabld.tlb(fla_config_gain),-
                     flabld.tlb(fla_model),-
                     csdr$library:futlib.tlb(fut_params),-
                     csdr$library:csdrlib.tlb(upm_stat_msg),-
                     fex_mod^,-
                     fex_nyquist^

fla_read_model.obj : flabld.tlb(fla_invoc_gain),-
                     flabld.tlb(fla_model),-
                     csdr$library:futlib.tlb(fut_params),-
                     fex_mod^

fla_read_reference.obj : flabld.tlb(fla_invoc_gain),-
                         flabld.tlb(fla_config_gain),-
                         csdr$library:ctuser.inc,-
                         csdr$library:futlib.tlb(fut_params),-
                         fex_nyquist^

fla_update_gain_report.obj : flabld.tlb(fla_invoc_gain),-
                             csdr$library:futlib.tlb(fut_error),-
                             csdr$library:futlib.tlb(fut_params)

fla_write_gain.obj : flabld.tlb(fla_invoc_gain),-
                     flabld.tlb(fla_config_gain),-
                     flabld.tlb(fla_model),-
                     fex_nyquist^,-
                     fex_mod^
