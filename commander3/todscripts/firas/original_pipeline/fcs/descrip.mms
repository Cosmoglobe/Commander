!
!   Purpose: Build facility fcs_combine_spectra.
!
!   Author: S. Alexander, STX, 5/21/92, SER 10960
!
!   Modifications: S. Alexander, HSTX, 11/29/93, SER 11417. Remove multiple
!                  scan modes from fcs.
!   Modifications: L. Rosen, HSTX, 8/16/94, fcs_include should be
!                  fut_fcs_include
!
.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         fcsbld/lib + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fcsbld.olb")'" .eqs. "" then $(libr)/create fcsbld
 if "''f$search("fcsbld.tlb")'" .eqs. "" then $(libr)/create/text fcsbld

fcsbld_obj = fcs,-
             fcs_avg,-
             fcs_close,-
             fcs_open,-
             fcs_parse,-
             fcs_report_check,-
             fcs_spec,-
             fcs_msg

fcsbld_txt = fcs_msg

fcs : fcs.exe
 ! Facility fcs has been built.

fcs.exe : fcsbld.olb($(fcsbld_obj)),-
          fcsbld.tlb($(fcsbld_txt)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fcsbld/lib/inc=fcs,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option

fcsbld.tlb(fcs_msg) : fcs_msg.txt
 $(libr) $(librflags)/text fcsbld fcs_msg

fcs.obj : fcsbld.tlb(fcs_msg),-
          csdr$library:futlib.tlb(fut_fcs_include),-
          csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:csdrlib.tlb(upm_stat_msg),-
          csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
          csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
          csdr$library:ctparams.inc,-
          fcf_sky^,-
          fcs_sky^

fcs_avg.obj : csdr$library:futlib.tlb(fut_fcs_include),-
              fcsbld.tlb(fcs_msg),-
              csdr$library:futlib.tlb(fut_params),-
              fcs_sky^

fcs_close.obj : fcsbld.tlb(fcs_msg),-
                csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params)
 
fcs_open.obj : fcsbld.tlb(fcs_msg),-
               csdr$library:futlib.tlb(fut_error),-
               csdr$library:futlib.tlb(fut_params)

fcs_parse.obj : fcsbld.tlb(fcs_msg),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg),-
                csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                ccm_cme_catalog_entry^

fcs_report_check.obj : fcsbld.tlb(fcs_msg),-
                       csdr$library:futlib.tlb(fut_error),-
                       csdr$library:futlib.tlb(fut_params)

fcs_spec.obj : fcsbld.tlb(fcs_msg),-
               csdr$library:futlib.tlb(fut_fcs_include),-
               fcf_sky^,-
               fcs_sky^
