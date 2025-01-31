function FMD_OST_ST, Nstruct
 
  common FMD_OST_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_ost = { FMD_OST                      , $
                   STR_GAMA: fltarr(182)        , $
                   STR_BETA: fltarr(6144)           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_OST}, Nstruct )
end
