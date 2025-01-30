function fla_fms_sky_st, Nstruct
 
  common fla_fms_sky_st, defined
 
  if N_elements( defined ) NE 1 then begin
 
       ct_head = { CT_HEAD                     , $
                   GMT: bytarr( 14 )           , $
                   TIME: Lonarr( 2 )           , $
                   SPACE_TIME: bytarr( 6 )     , $
                   MJR_FRM_NO: 0L              , $
                   ORBIT: 0L                   , $
                   HSKP1_TLM_FMT: 0B           , $
                   HSKP2_TLM_FMT: 0B           , $
                   INGEST_SPARES: bytarr( 18 ) , $
                   DATASET_ID: 0               , $
                   INSTR_SPARES: bytarr( 6 )       }
 
       coad_spec_head = { COAD_SPEC_HEAD                , $
                          FIRST_GMT: bytarr( 14 )       , $
                          FIRST_TIME: Lonarr( 2 )       , $
                          FIRST_SPACE_TIME: bytarr( 6 ) , $
                          FIRST_MJR_FRM_NO: 0L          , $
                          LAST_GMT: bytarr( 14 )        , $
                          LAST_TIME: Lonarr( 2 )        , $
                          LAST_SPACE_TIME: bytarr( 6 )  , $
                          LAST_MJR_FRM_NO: 0L           , $
                          COADD_NO: 0L                  , $
                          NUM_IFGS: 0                   , $
                          LABEL: bytarr( 60 )           , $
                          COMB_NUM_IFGS: 0.0            , $
                          COMB_CHI_SQUARE: 0.0          , $
                          IMAG_COMB_CSQ: 0.0            , $
                          SPARES: bytarr( 4 )              }
 
       prim_template = { PRIM_TEMPLATE  , $
                         SUBTRACTED: 0B , $
                         AMPLITUDE: 0.0 , $
                         SNR: 0.0           }
 
       sec_template = { SEC_TEMPLATE    , $
                        SUBTRACTED: 0B  , $
                        AMPLITUDE: 0.0  , $
                        SNR: 0.0        , $
                        VARIANCE: 0.0   , $
                        B_AVERAGE: 0.0  , $
                        B_VARIANCE: 0.0     }
 
       transient = { TRANSIENT           , $
                     C_AVERAGE: 0.0      , $
                     C_VARIANCE: 0.0     , $
                     BL_TRANS_COEFF: 0.0     }
 
       deglitch = { DEGLITCH           , $
                    GLITCH_ITER: 0     , $
                    GLITCH_SIGNAL: 0.0     }
 
       coad_spec_data = { COAD_SPEC_DATA              , $
                          CHAN_ID: 0B                 , $
                          MTM_SPEED: 0B               , $
                          MTM_LENGTH: 0B              , $
                          FAKEIT: 0B                  , $
                          SCI_MODE: 0B                , $
                          ADDS_PER_GROUP: 0B          , $
                          XCAL_POS: 0B                , $
                          SWEEPS: 0B                  , $
                          GAIN_SUM: 0L                , $
                          GLITCH_RATE: 0.0            , $
                          DQ_SUMMARY_FLAG: 0B         , $
                          ATT_SUMMARY_FLAG: 0B        , $
                          PEAK_POS: 0                 , $
                          NYQUIST_HERTZ: 0.0          , $
                          NYQUIST_ICM: 0.0            , $
                          PRIM_TEMPLATE:prim_template , $
                          SEC_TEMPLATE:sec_template   , $
                          TRANSIENT:transient         , $
                          DEGLITCH:deglitch           , $
                          NOISE: 0.0                  , $
                          BIN_INFO: intarr( 256 )     , $
                          BL_COEFFS: fltarr( 5 )      , $
                          BOL_CMD_BIAS: 0             , $
                          BOL_VOLT: 0.0               , $
                          BOL_VOLT_MIN: fltarr( 4 )   , $
                          BOL_VOLT_MAX: fltarr( 4 )   , $
                          TEMP: fltarr(10)            , $
                          TEMP_SIGMA: fltarr(10)      , $
                          TEMP_MIN: fltarr(10)        , $
                          TEMP_MAX: fltarr(10)        , $
                          SPARES: bytarr( 14 )            }

       spec_data = { SPEC_DATA                        , $
                     MODEL_TTAG: bytarr( 14 )         , $
                     MODEL_LABEL: bytarr( 40 )        , $
                     COMBINE_SCRIPT: bytarr( 80 )     , $
                     RESPONSIVITY: 0.0                , $
                     RESP_SIGMA: 0.0                  , $
                     TIME_CONSTANT: 0.0               , $
                     TC_SIGMA: 0.0                    , $
                     SPEC: complexarr(257)            , $
                     REAL_VAR: fltarr(257)            , $
                     IMAG_VAR: fltarr(257)            , $
                     CVS_VAR: fltarr(257)             , $
                     PHASE_CORR: 0.0                  , $
                     PC_SIGMA: 0.0                    , $
                     QRAD: 0.0                        , $
                     QRAD_SIGMA: 0.0                  , $
                     IR_POWER: 0.0                    , $
                     IR_POWER_SIGMA: 0.0              , $
                     DESTRIPED: 0B                    , $
                     SPARES: bytarr( 135 )                }

       en_stat = { EN_STAT                            , $
                          STAT_WORD_1:0               , $
                          INT_REF_TEMP_A:0            , $
                          REF_HRN_TEMP_A:0            , $
                          STAT_WORD_4:0               , $
                          STAT_WORD_5:0               , $
                          SKY_HRN_TEMP_A:0            , $
                          EXT_CAL_TEMP_A:0            , $
                          STAT_WORD_8:0               , $
                          STAT_WORD_9:0               , $
                          INT_REF_TEMP_B:0            , $
                          REF_HRN_TEMP_B:0            , $
                          STAT_WORD_12:0              , $
                          STAT_WORD_13:0              , $
                          SKY_HRN_TEMP_B:0            , $
                          EXT_CAL_TEMP_B:0            , $
                          STAT_WORD_16:0              , $
                          GRT_ADDR:bytarr(2)          , $
                          MICRO_STAT_BUS:bytarr(4)    , $
                          BOL_CMD_BIAS:bytarr(4)      , $
                          DWELL_STAT:bytarr(2)        , $
                          LVDT_STAT:bytarr(2)         , $
                          HOT_SPOT_CMD:bytarr(2)      , $
                          ENGSTAT_SPARES1:bytarr(10)  , $
                          POWER_A_STATUS:bytarr(2)    , $
                          POWER_B_STATUS:bytarr(2)    , $
                          ENGSTAT_SPARES2:bytarr(6)     }

       en_analog = { EN_ANALOG                  , $
                     A_LO_GRT: fltarr( 16 )     , $
                     A_HI_GRT: fltarr( 16 )     , $
                     B_LO_GRT: fltarr( 16 )     , $
                     B_HI_GRT: fltarr( 16 )     , $
                     TEMP_CTRL: fltarr( 8 )     , $
                     IPDU_TEMP: fltarr( 2 )     , $
                     CNA_TEMP: fltarr( 4 )      , $
                     DBX_TEMP: fltarr( 2 )      , $
                     STAT_MON_TEMP: fltarr( 2 ) , $
                     PAMP_CHAN: 0.0             , $
                     PAMP_OP: 0.0               , $
                     HOT_SPOT: fltarr( 2 )      , $
                     MTM_CAL_MTR: fltarr( 2 )   , $
                     MTM_POS: fltarr( 2 )       , $
                     BMTR_VOLT: fltarr( 4 )     , $
                     IPDU_VOLT: fltarr( 20 )    , $
                     IPDU_AMP: fltarr( 12 )         }

       en_sigma = { EN_SIGMA                        , $
                     SIG_A_LO_GRT: fltarr(16)       , $
                     SIG_A_HI_GRT: fltarr(16)       , $
                     SIG_B_LO_GRT: fltarr(16)       , $
                     SIG_B_HI_GRT: fltarr(16)       , $
                     SIG_TEMP_CTRL: fltarr(8)       , $
                     SIG_IPDU_TEMP: fltarr(2)       , $
                     SIG_CNA_TEMP: fltarr(4)        , $
                     SIG_DBX_TEMP: fltarr(2)        , $
                     SIG_STAT_MON_TEMP: fltarr(2)   , $
                     SIG_PAMP_CHAN:0.0              , $
                     SIG_PAMP_OP:0.0                , $
                     SIG_HOT_SPOT: fltarr(2)        , $
                     SIG_MTM_CAL_MTR: fltarr(2)     , $
                     SIG_MTM_POS: fltarr(2)         , $
                     SIG_BOL_VOLT: fltarr(4)        , $
                     SIG_IPDU_VOLT: fltarr(20)      , $
                     SIG_IPDU_AMP: fltarr(12)         }

       en_tempdiff = { EN_TEMPDIFF              , $
                       XCAL:0                   , $
                       ICAL:0                   , $
                       SKYHORN:0                , $
                       REFHORN:0                , $
                       DIHEDRAL:0               , $
                       COLLIMATOR_MIRROR:0      , $
                       BOL_ASSEM: intarr(4)       }

       attitude = { ATTITUDE                , $
                    PIXEL_NO: 0L            , $
                    CELESTIAL: fltarr( 3 )  , $
                    RA: 0                   , $
                    DEC: 0                  , $
                    TERR_PIXEL_NO: 0L       , $
                    TERR_LAT: 0             , $
                    TERR_LONG: 0            , $
                    EARTH_LIMB: 0           , $
                    EARTH_LIMB_AZIMUTH: 0   , $
                    SUN_ANGLE: 0            , $
                    MOON_ANGLE: 0           , $
                    MOON_AZ_ANGLE: 0        , $
                    MOON_PHASE: 0           , $
                    SUN_MOON_DIST: 0.0      , $
                    COBE_MOON_DIST: 0.0     , $
                    ALTITUDE: 0             , $
                    PROJ_HELIO_VELOC: 0     , $
                    MCILWAIN_L_PARAM: 0     , $
                    GALACTIC_LONG: 0        , $
                    GALACTIC_LAT: 0         , $
                    ECLIPTIC_LONG: 0        , $
                    ECLIPTIC_LAT: 0         , $
                    ORBITAL_PHASE: 0        , $
                    PROJ_GEO_VELOC: 0       , $
                    SCAN_ANGLE: 0           , $
                    SC_ROTATION_ANGLE: 0    , $
                    SOLUTION: 0B            , $
                    PIXEL_DEFINITION: 0B    , $
                    SKYMAP_INDEX: 0         , $
                    EXC_GALACTIC_LAT: 0     , $
                    TERR_RAD_BYTE: 0B       , $
                    OUTSIDE_GALAXY: 0B      , $
                    ATT_SPARES: bytarr( 2 )     }
 
       fms_sky = { FMS_SKY                       , $
                   CT_HEAD:ct_head               , $
                   COAD_SPEC_HEAD:coad_spec_head , $
                   COAD_SPEC_DATA:coad_spec_data , $
                   SPEC_DATA:spec_data           , $
                   EN_STAT:en_stat               , $
                   EN_ANALOG:en_analog           , $
                   EN_SIGMA:en_sigma             , $
                   EN_TEMPDIFF1:en_tempdiff      , $
                   EN_TEMPDIFF2:en_tempdiff      , $
                   ATTITUDE:attitude                 }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMS_SKY}, Nstruct )
end
