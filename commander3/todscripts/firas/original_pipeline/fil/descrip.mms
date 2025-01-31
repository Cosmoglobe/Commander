!   FIRAS_Interferogram_Long (FIL) Descrip MMS File
!
!   Purpose: Build facility fil_interferogram_long.
!
!   Author: S. Brodd, HSTX, 4/95, SER 12244
!
!   Modifications:

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("filbld.olb")'" .eqs. "" then $(libr)/create filbld

filbld_obj = fil,-
             fil_baseline_sub,-
             fil_close,-
             fil_close_summary,-
             fil_coadd,-
             fil_convert,-
             fil_covar,-
             fil_deglitch,-
             fil_find_glitch,-
             fil_inst_state,-
             fil_noise,-
             fil_open_input,-
             fil_open_output,-
             fil_open_report,-
             fil_open_summary,-
             fil_parse,-
             fil_qual_check,-
             fil_read,-
             fil_read_covar,-
             fil_reference,-
             fil_sec_template,-
             fil_shape_check,-
             fil_short,-
             fil_sort,-
             fil_template,-
             fil_transient,-
             fil_write,-
             fil_write_covar,-
             fil_write_report,-
             fil_msg

fil : fil.exe
 ! Facility fil has been built.

fil.exe : filbld.olb($(filbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) filbld/lib/inc=fil,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      csdr$library:v5.opt/option

fil.obj : csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctparams.inc,-
          fdq_eng^,-
          fdq_sdf^,-
          fex_basis^,-
          fex_cmdgain^,-
          fex_cth^,-
          fex_gltchcor^,-
          fex_grtcoawt^,-
          fex_grtrawwt^,-
          fex_grttrans^,-
          fex_mincoadd^,-
          fex_nyquistl^,-
          fex_samprate^,-
          fil_cov^,-
          fil_scc^,-
          fil_sky^

fil_baseline_sub : fex_basis^,-
                   fil_sky^

fil_close : csdr$library:futlib.tlb(fut_params),-
            csdr$library:ctparams.inc

fil_close_summary : csdr$library:futlib.tlb(fut_error),-
                    csdr$library:futlib.tlb(fut_params),-
                    fex_cth^,-
                    fex_mincoadd^

fil_coadd : csdr$library:futlib.tlb(fut_params),-
            fdq_eng^,-
            fdq_sdf^,-
            fex_cth^,-
            fex_gltchcor^,-
            fex_etfl^,-
            fex_grtcoawt^,-
            fex_grttrans^,-
            fex_nyquistl^,-
            fil_cov^,-
            fil_scc^,-
            fil_sky^

fil_convert : fdq_sdf^,-
              fil_scc^

fil_covar : fil_cov^

fil_deglitch : csdr$library:futlib.tlb(fut_params),-
               fex_gltchpro^,-
               fex_mincoadd^,-
               fil_scc^,-
               fil_sky^

fil_inst_state : csdr$library:futlib.tlb(fut_params),-
                 fdq_eng^,-
                 fdq_sdf^,-
                 fex_cth^,-
                 fex_grtrawwt^,-
                 fex_grttrans^,-
                 fex_mincoadd^,-
                 fil_scc^,-
                 fut_enganlg^

fil_noise : csdr$library:futlib.tlb(fut_params),-
            fil_scc^,-
            fil_sky^

fil_open_input : csdr$library:futlib.tlb(fut_params),-
                 csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                 csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                 csdr$library:ctparams.inc,-
                 ccm_cme_catalog_entry^

fil_open_output : csdr$library:futlib.tlb(fut_params)

fil_open_report : csdr$library:futlib.tlb(fut_params),-
                  fex_cth^,-
                  fex_mincoadd^

fil_open_summary : csdr$library:futlib.tlb(fut_error),-
                   csdr$library:futlib.tlb(fut_params)

fil_parse.obj : csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg)

fil_qual_check : csdr$library:futlib.tlb(fut_params),-
                 fdq_eng^,-
                 fdq_sdf^,-
                 fex_cmdgain^,-
                 fex_grtcoawt^,-
                 fex_grttrans^,-
                 fex_mincoadd^,-
                 fex_nyquistl^,-
                 fex_samprate^,-
                 fil_scc^,-
                 fil_sky^

fil_read : csdr$library:futlib.tlb(fut_params),-
           csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
           csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
           csdr$library:ctparams.inc,-
           fdq_eng^,-
           fdq_sdf^,-
           fss_sssky^

fil_read_covar : csdr$library:futlib.tlb(fut_params),-
                 csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                 csdr$library:ctparams.inc,-
                 ccm_cme_catalog_entry^,-
                 fil_cov^

fil_reference : csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(cct_get_config),-
                fex_basis^,-
                fex_cmdgain^,-
                fex_cth^,-
                fex_gltchcor^,-
                fex_grtcoawt^,-
                fex_grtrawwt^,-
                fex_grttrans^,-
                fex_mincoadd^,-
                fex_nyquistl^,-
                fex_samprate^

fil_sec_template : csdr$library:futlib.tlb(fut_params),-
                   fex_cth^,-
                   fex_reftemps^,-
                   fil_scc^,-
                   fil_sky^

fil_shape_check : fex_cth^,-
                  fil_scc^

fil_short : csdr$library:futlib.tlb(fut_params),-
            fex_nyquistl^,-
            fil_sky^

fil_sort : csdr$library:futlib.tlb(fut_params),-
           fdq_sdf^,-
           fex_cth^,-
           fil_scc^

fil_template : csdr$library:futlib.tlb(fut_params),-
               fil_scc^,-
               fil_sky^

fil_transient : fex_dtrf^,-
                fil_scc^,-
                fil_sky^

fil_write : csdr$library:futlib.tlb(fut_params),-
            csdr$library:ctparams.inc,-
            fil_cov^,-
            fil_scc^,-
            fil_sky^

fil_write_covar : csdr$library:futlib.tlb(fut_params),-
                  fil_cov^,-
                  fil_scc^,-
                  fil_sky^

fil_write_report : csdr$library:futlib.tlb(fut_params),-
                   fex_cth^,-
                   fex_mincoadd^,-
                   fil_scc^
