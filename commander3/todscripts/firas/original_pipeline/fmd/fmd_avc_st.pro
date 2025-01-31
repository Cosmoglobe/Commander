function FMD_AVC_ST, Nstruct
 
  common FMD_AVC_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_avc = { FMD_AVC                                 , $
                   A_VECTOR: fltarr( 32 )                      }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_AVC}, Nstruct )
end
