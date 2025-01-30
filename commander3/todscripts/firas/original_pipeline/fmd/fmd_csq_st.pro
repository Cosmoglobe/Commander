function FMD_CSQ_ST, Nstruct
 
  common FMD_CSQ_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_csq = { FMD_CSQ             , $
                   PIXEL: 0L           , $
                   GAL_LON: 0.0        , $
                   GAL_LAT: 0.0        , $
                   WEIGHT: 0.0         , $
                   DOF: 0L             , $
                   COMB_CSQ: 0.0           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_CSQ}, Nstruct )
end
