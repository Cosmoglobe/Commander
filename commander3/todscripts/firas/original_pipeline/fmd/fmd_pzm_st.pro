function FMD_PZM_ST, Nstruct
 
  common FMD_PZM_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_pzm = { FMD_PZM                            , $
                   PIXEL: 0L                          , $
                   GAL_LON: 0.0                       , $
                   GAL_LAT: 0.0                       , $
                   ZODI_MOD: fltarr( 170 )                }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_PZM}, Nstruct )
end
