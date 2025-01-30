function FEX_VABSAA_st, Nstruct
 
  common FEX_VABSAA_st, defined
 
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
 
       vab = { VAB                , $
               LONSTEP: 0.0       , $
               LATMIN: 0.0        , $
               LATMAX: 0.0        , $
               LATN: fltarr( 41 ) , $
               LATS: fltarr( 41 )     }
 
       saa = { SAA                , $
               LONSTEP: 0.0       , $
               LONMIN: 0.0        , $
               LONMAX: 0.0        , $
               LATMIN: 0.0        , $
               LATMAX: 0.0        , $
               LATN: fltarr( 31 ) , $
               LATS: fltarr( 31 )     }
 
       fex_vabsaa = { FEX_VABSAA               , $
                      CT_HEAD:ct_head          , $
                      VAB: replicate( vab, 2 ) , $
                      SAA:saa                  , $
                      SPARES: bytarr( 12 )         }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_VABSAA}, Nstruct )
end
