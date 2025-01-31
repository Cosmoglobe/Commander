function FMD_CVC_ST, Nstruct
 
  common FMD_CVC_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_cvc = { FMD_CVC                                 , $
                   C_VECTOR: fltarr( 182 )                     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_CVC}, Nstruct )
end
