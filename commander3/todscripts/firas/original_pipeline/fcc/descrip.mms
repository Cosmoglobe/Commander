!   FIRAS_Calibrate_Covariances (FCC) Descrip MMS File
!
!   Purpose: Build facility fcc_calibrate_covariances.
!
!   Author: S. Alexander, HSTX, 7/93, SER 11189
!
!   Modifications:
!

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("fccbld.olb")'" .eqs. "" then $(libr)/create fccbld

fccbld_obj = fcc,-
             fcc_calibrate,-
             fcc_extract,-
             fcc_gain,-
             fcc_insert,-
             fcc_parse,-
             fcc_read,-
             fcc_reference,-
             fcc_report,-
             fcc_write,-
             fcc_msg

fcc : fcc.exe
 ! Facility fcc has been built.

fcc.exe : fccbld.olb($(fccbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fccbld/lib/inc=fcc,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      csdr$library:v5.opt/option

fcc.obj : csdr$library:futlib.tlb(fut_error),-
          csdr$library:ctparams.inc,-
          fcc_cov^,-
          fex_mod^,-
          fic_cov^

fcc_extract : csdr$library:futlib.tlb(fut_error),-
              fcc_cov^,-
              fex_mod^,-
              fic_cov^

fcc_gain : csdr$library:futlib.tlb(fut_params),-
           fex_mod^

fcc_insert : fcc_cov^

fcc_parse.obj : csdr$library:csdrlib.tlb(upm_stat_msg)

fcc_read : csdr$library:csdrlib.tlb(cct_query_catalog_record),-
           csdr$library:ctparams.inc,-
           ccm_cme_catalog_entry^,-
           fic_cov^

fcc_reference : csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(cct_get_config),-
                fex_apod^,-
                fex_etf^,-
                fex_mod^,-
                fex_nyquist^

fcc_report : csdr$library:futlib.tlb(fut_error),-
             csdr$library:futlib.tlb(fut_params)

fcc_write : csdr$library:ctparams.inc,-
            fcc_cov^
