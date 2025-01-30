!   FIRAS_Interferogram_Coaddition (FIC) Descrip MMS File
!
!   Purpose: Build facility fic_interferogram_coadd.
!
!   Author: J. Durachta, STI, 10/17/86
!
!   Modifications:
!
!   SPR 6004, R. Kummerer: VMS 5.2 options file changes.
!   SER 7985, S. Alexander, 9/1/91: Implement new fic requirements.

.suffixes
.suffixes .exe .olb .obj .for .msg

.for.obj
 $(fort) $(fflags)/extend_source/cont=99 $(mms$source) + -
                                         csdr$library:futlib/lib + -
                                         csdr$library:csdrlib/lib

.obj.olb
 $(libr) $(librflags) $(mms$target) $(mms$source)

.first
 if "''f$search("ficbld.olb")'" .eqs. "" then $(libr)/create ficbld

ficbld_obj = fic,-
             fic_baseline_sub,-
             fic_close,-
             fic_close_summary,-
             fic_coadd,-
             fic_convert,-
             fic_covar,-
             fic_deglitch,-
             fic_find_glitch,-
             fic_inst_state,-
             fic_noise,-
             fic_open_input,-
             fic_open_output,-
             fic_open_report,-
             fic_open_summary,-
             fic_parse,-
             fic_qual_check,-
             fic_read,-
             fic_read_covar,-
             fic_reference,-
             fic_sec_template,-
             fic_shape_check,-
             fic_template,-
             fic_transient,-
             fic_write,-
             fic_write_covar,-
             fic_write_report,-
             fic_msg

fic : fic.exe
 ! Facility fic has been built.

fic.exe : ficbld.olb($(ficbld_obj)),-
          csdr$library:futlib.olb,-
          csdr$library:csdrlib.olb,-
          imsl$dir:imsl.olb,-
          csdr$library:v5.opt
 $(link) $(linkflags) ficbld/lib/inc=fic,-
                      csdr$library:futlib/lib,-
                      csdr$library:csdrlib/lib,-
                      imsl$dir:imsl/lib,-
                      csdr$library:v5.opt/option

fic.obj : csdr$library:futlib.tlb(fut_error),-
          csdr$library:futlib.tlb(fut_params),-
          csdr$library:ctparams.inc,-
          fdq_eng^,-
          fdq_sdf^,-
          fex_basis^,-
          fex_cmdgain^,-
          fex_cth^,-
          fex_grtcoawt^,-
          fex_grtrawwt^,-
          fex_grttrans^,-
          fex_mincoadd^,-
          fex_nyquist^,-
          fex_samprate^,-
          fic_cov^,-
          fic_scc^,-
          fic_sky^

fic_baseline_sub : fex_basis^,-
                   fic_sky^

fic_close : csdr$library:futlib.tlb(fut_params),-
            csdr$library:ctparams.inc

fic_close_summary : csdr$library:futlib.tlb(fut_error),-
                    csdr$library:futlib.tlb(fut_params),-
                    fex_cth^,-
                    fex_mincoadd^

fic_coadd : csdr$library:futlib.tlb(fut_params),-
            fdq_eng^,-
            fdq_sdf^,-
            fex_cth^,-
            fex_etf^,-
            fex_grtcoawt^,-
            fex_grttrans^,-
            fic_cov^,-
            fic_scc^,-
            fic_sky^

fic_convert : fdq_sdf^,-
              fic_scc^

fic_covar : fic_cov^

fic_deglitch : csdr$library:futlib.tlb(fut_params),-
               fex_gltchpro^,-
               fex_mincoadd^,-
               fic_scc^,-
               fic_sky^

fic_inst_state : csdr$library:futlib.tlb(fut_params),-
                 fdq_eng^,-
                 fdq_sdf^,-
                 fex_cth^,-
                 fex_grtrawwt^,-
                 fex_grttrans^,-
                 fex_mincoadd^,-
                 fic_scc^,-
                 fut_enganlg^

fic_noise : fic_scc^,-
            fic_sky^

fic_open_input : csdr$library:futlib.tlb(fut_params),-
                 csdr$library:csdrlib.tlb(cct_query_catalog_record),-
                 csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                 csdr$library:ctparams.inc,-
                 ccm_cme_catalog_entry^

fic_open_output : csdr$library:futlib.tlb(fut_params)

fic_open_report : csdr$library:futlib.tlb(fut_params),-
                  fex_cth^,-
                  fex_mincoadd^

fic_open_summary : csdr$library:futlib.tlb(fut_error),-
                   csdr$library:futlib.tlb(fut_params)

fic_parse.obj : csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(upm_stat_msg),-
                csdr$library:ctparams.inc

fic_qual_check : csdr$library:futlib.tlb(fut_params),-
                 fdq_eng^,-
                 fdq_sdf^,-
                 fex_cmdgain^,-
                 fex_grtcoawt^,-
                 fex_grttrans^,-
                 fex_mincoadd^,-
                 fex_nyquist^,-
                 fex_samprate^,-
                 fic_scc^,-
                 fic_sky^

fic_read : csdr$library:futlib.tlb(fut_params),-
           csdr$library:csdrlib.tlb(csa_pixel_input_rec),-
           csdr$library:csdrlib.tlb(csa_pixel_output_rec),-
           csdr$library:ctparams.inc,-
           fdq_eng^,-
           fdq_sdf^,-
           fss_sssky^

fic_read_covar : csdr$library:futlib.tlb(fut_params),-
                 csdr$library:csdrlib.tlb(cct_query_ttg_catalog_record),-
                 csdr$library:ctparams.inc,-
                 ccm_cme_catalog_entry^,-
                 fic_cov^

fic_reference : csdr$library:futlib.tlb(fut_error),-
                csdr$library:futlib.tlb(fut_params),-
                csdr$library:csdrlib.tlb(cct_get_config),-
                fex_basis^,-
                fex_cmdgain^,-
                fex_cth^,-
                fex_grtcoawt^,-
                fex_grtrawwt^,-
                fex_grttrans^,-
                fex_mincoadd^,-
                fex_nyquist^,-
                fex_samprate^

fic_sec_template : csdr$library:futlib.tlb(fut_params),-
                   fex_cth^,-
                   fex_reftemps^,-
                   fic_scc^,-
                   fic_sky^

fic_shape_check : fex_cth^,-
                  fic_scc^

fic_template : csdr$library:futlib.tlb(fut_params),-
               fic_scc^,-
               fic_sky^

fic_transient : fex_dtrf^,-
                fic_scc^,-
                fic_sky^

fic_write : csdr$library:futlib.tlb(fut_params),-
            csdr$library:ctparams.inc,-
            fic_cov^,-
            fic_scc^,-
            fic_sky^

fic_write_covar : csdr$library:futlib.tlb(fut_params),-
                  fic_cov^,-
                  fic_scc^,-
                  fic_sky^

fic_write_report : csdr$library:futlib.tlb(fut_params),-
                   fex_cth^,-
                   fex_mincoadd^,-
                   fic_scc^
