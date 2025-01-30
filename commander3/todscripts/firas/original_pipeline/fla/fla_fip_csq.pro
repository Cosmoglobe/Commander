Function FLA_FIP_CSQ, Nstruct
 
  common FLA_FIP_CSQ, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fip_csq = { FIP_CSQ                                 , $
                   OUTMAP: string(' ',format='(A4)')       , $
                   INMAP_1: string(' ',format='(A4)')      , $
                   INMAP_2: string(' ',format='(A4)')      , $
                   PIXEL: 0L                               , $
                   WEIGHT_1: 0.0                           , $
                   WEIGHT_2: 0.0                           , $
                   NUM_FREQ: 0L                            , $
                   COMB_CHI_SQUARE: 0.0                         } 
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FIP_CSQ}, Nstruct )
end
