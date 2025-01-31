function FRD_FEX_MCS, Nstruct
 
  common FRD_FEX_MCS, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fex_mcs = { FEX_MCS                             , $
                   CHANSCAN: string(' ',format='(A4)') , $
                   OFFSET_SPEC: complexarr( 257 )      , $
                   GAIN_SPEC: fltarr( 257 )                }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_MCS}, Nstruct )
end
