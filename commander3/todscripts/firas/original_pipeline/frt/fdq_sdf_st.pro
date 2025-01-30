function FDQ_SDF_st, Nstruct
 
  common FDQ_SDF_st, defined
 
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
 
       sci_head = { SCI_HEAD                , $
                    CHAN_ID: 0              , $
                    GAIN: 0                 , $
                    MTM_SPEED: 0            , $
                    MTM_LENGTH: 0           , $
                    DATA_QUAL: bytarr( 60 ) , $
                    DATA_READY: intarr( 8 ) , $
                    SC_HEAD0: 0             , $
                    SC_HEAD1A: 0B           , $
                    SC_HEAD1B: 0B           , $
                    SC_HEAD2: 0             , $
                    SC_HEAD3: 0             , $
                    SC_HEAD4: 0             , $
                    SC_HEAD5: 0             , $
                    SC_HEAD6: 0             , $
                    SC_HEAD7: 0             , $
                    SC_HEAD8: 0             , $
                    SC_HEAD9: 0             , $
                    SC_HEAD10: 0            , $
                    SC_HEAD11: 0            , $
                    SC_HEAD12: 0            , $
                    SC_HEAD13: 0            , $
                    SC_HEAD14: 0            , $
                    SC_HEAD15: 0            , $
                    SC_HEAD16: 0            , $
                    SC_HEAD17: 0            , $
                    SC_HEAD18: 0            , $
                    SC_HEAD19: 0            , $
                    SC_HEAD20: 0            , $
                    SC_HEAD21: 0            , $
                    SC_HEAD22: 0            , $
                    SC_HEAD23: 0            , $
                    SC_HEAD24: 0            , $
                    SC_HEAD25: 0                }
 
       ifg_data = { IFG_DATA            , $
                    IFG: intarr( 512 )  , $
                    GLTCH: intarr( 32 )     }
 
       dq_data = { DQ_DATA                     , $
                   FAKE: 0                     , $
                   XCAL_POS: 0                 , $
                   IREF_TEMP: 0.0              , $
                   ENG_TIME: Lonarr( 2 )       , $
                   ENG_REC: 0L                 , $
                   DATA_QUALITY: bytarr( 110 ) , $
                   IFG_NO: 0L                  , $
                   DQ_SPARES: bytarr( 24 )         }
 
       collect_time = { COLLECT_TIME               , $
                        MIDPOINT_TIME: Lonarr( 2 ) , $
                        BADTIME_FLAG: 0B           , $
                        FPP_SPARE: 0B                  }
 
       attitude = { ATTITUDE                          , $
                    PIXEL_NO: 0L                      , $
                    EQUATORIAL: fltarr( 3 )           , $
                    RA: 0                             , $
                    DEC: 0                            , $
                    TERR_PIXEL_NO: 0L                 , $
                    TERR_LAT: 0                       , $
                    TERR_LONG: 0                      , $
                    EARTH_LIMB: 0                     , $
                    EARTH_LIMB_AZIM: 0                , $
                    SUN_ANGLE: 0                      , $
                    MOON_ANGLE: 0                     , $
                    MOON_AZ_ANGLE: 0                  , $
                    MOON_PHASE: 0                     , $
                    SUN_MOON_DIST: 0.0                , $
                    COBE_MOON_DIST: 0.0               , $
                    ALTITUDE: 0                       , $
                    PROJ_BARY_VELOC: 0                , $
                    MCILWAIN_L_PARM: 0                , $
                    GALACTIC_LONG: 0                  , $
                    GALACTIC_LAT: 0                   , $
                    ECLIPTIC_LONG: 0                  , $
                    ECLIPTIC_LAT: 0                   , $
                    ORBITAL_PHASE: 0                  , $
                    PROJ_GEO_VELOC: 0                 , $
                    SCAN_ANGLE: 0                     , $
                    SC_ROTATION_ANG: 0                , $
                    SOLUTION: 0B                      , $
                    PIX_DEFINITION: 0B                , $
                    SKYMAP_INDEX: 0                   , $
                    EXC_GALACT_LAT: 0                 , $
                    TERR_RAD_BYTE: 0B                 , $
                    ATT_SPARES: bytarr( 3 )               }
 
       fdq_sdf = { FDQ_SDF                   , $
                   CT_HEAD:ct_head           , $
                   SCI_HEAD:sci_head         , $
                   IFG_DATA:ifg_data         , $
                   DQ_DATA:dq_data           , $
                   COLLECT_TIME:collect_time , $
                   ATTITUDE:attitude             }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FDQ_SDF}, Nstruct )
end
