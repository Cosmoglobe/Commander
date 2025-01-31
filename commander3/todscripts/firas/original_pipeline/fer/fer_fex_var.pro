function FER_FEX_VAR, Nstruct
 
  common FER_FEX_VAR, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fex_var = { FEX_VAR                   , $
                   GMT: bytarr( 14 )         , $
                   TIME: Lonarr( 2 )         , $
                   CHANNEL: 0B               , $
                   SCAN_LENGTH: 0B           , $
                   SCAN_SPEED: 0B            , $
                   MODEL_LABEL: bytarr( 40 ) , $
                   GALAT_EXC: 0.0            , $
                   NSKY_IFGS: 0.0            , $
                   VARIANCES: dblarr( 257 )  , $
                   SPARES: bytarr( 47 )          }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_VAR}, Nstruct )
end
