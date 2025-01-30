!   FIRAS_CalibrateCovariances_Long (FCL) Descrip MMS File
!
!   Purpose: Build facility fcl_calibratecovariances_long.
!
!   Author: S. Brodd, HSTX, 12/95, SPR 12291
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
 if "''f$search("fclbld.olb")'" .eqs. "" then $(libr)/create fclbld

fclbld_obj = fcl,-
             fcl_calibrate,-
             fcl_extract,-
             fcl_gain,-
             fcl_insert,-
             fcl_parse,-
             fcl_read,-
             fcl_reference,-
             fcl_report,-
             fcl_write,-
             fcl_msg

fcl : fcl.exe
 ! Facility fcl has been built.

fcl.exe : fclbld.olb($(fclbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) fclbld/lib/inc=fcl,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      csdr$library:v5.opt/option

fcl.obj : csdr$library:futlib.tlb(fut_error),-
          csdr$library:ctparams.inc,-
          fcl_cov^,-
          fex_mod^,-
          fil_cov^

fcl_extract.obj : csdr$library:futlib.tlb(fut_error),-
                  fcl_cov^,-
                  fex_mod^,-
                  fil_cov^

fcl_gain.obj : csdr$library:futlib.tlb(fut_params),-
               fex_mod^

fcl_insert.obj : fcl_cov^

fcl_parse.obj : csdr$library:csdrlib.tlb(upm_stat_msg)

fcl_read.obj : csdr$library:csdrlib.tlb(cct_query_catalog_record),-
               csdr$library:ctparams.inc,-
               ccm_cme_catalog_entry^,-
               fil_cov^

fcl_reference.obj : csdr$library:futlib.tlb(fut_params),-
                    csdr$library:csdrlib.tlb(cct_get_config),-
                    fex_mod^,-
                    fex_nyquistl^

fcl_report.obj : csdr$library:futlib.tlb(fut_error),-
                 csdr$library:futlib.tlb(fut_params)

fcl_write.obj : csdr$library:ctparams.inc,-
                fcl_cov^
