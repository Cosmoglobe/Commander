function FEX_GRTTRANS_st, Nstruct
 
;+
; Author: K.Jensen (STX) - Oct. 28, 1991 - SPR 9207.
;-

  common FEX_GRTTRANS_st, defined
 
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
 
       fex_grttrans = { FEX_GRTTRANS                   , $
                        CT_HEAD:ct_head                , $
                        NUM_GRT: 0L                    , $
                        GRT_A_TRAN_TEMP: fltarr( 16 ) , $
                        GRT_B_TRAN_TEMP: fltarr( 16 ) , $
                        TEMP_SPARES: 0.0               , $
                        GRT_A_TRAN_HWID: fltarr( 16 ) , $
                        GRT_B_TRAN_HWID: fltarr( 16 ) , $
                        HWID_SPARES: 0.0               , $
                        GRT_SPARES: bytarr( 52 )           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_GRTTRANS}, Nstruct )
end
