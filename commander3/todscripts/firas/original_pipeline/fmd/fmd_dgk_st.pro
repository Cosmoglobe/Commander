function FMD_DGK_ST, Nstruct
 
  common FMD_DGK_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_dgk = { FMD_DGK                              , $
                   PIXEL: 0L                            , $
                   CHANSCAN: string(' ',format='(A4)')  , $
                   REC_NUM: 0L                          , $
                   WEIGHT: 0.0                          , $
                   BAND_8: 0.0                          , $
                   BAND_9: 0.0                          , $
                   BAND_10: 0.0                             }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_DGK}, Nstruct )
end
