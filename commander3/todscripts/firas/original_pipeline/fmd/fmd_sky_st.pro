function FMD_SKY_ST, Nstruct
 
  common FMD_SKY_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_sky = { FMD_SKY                           , $
                   PIXEL: 0L                         , $
                   GAL_LON: 0.0                      , $
                   GAL_LAT: 0.0                      , $
                   WEIGHT: 0.0                       , $
                   STR_USED: 0B                      , $
                   SPECTRUM: fltarr( 182 )               }

       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_SKY}, Nstruct )
end
