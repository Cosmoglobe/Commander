function FDS_FIP_DSQ, Nstruct
 
  common FDS_FIP_DSQ, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fip_dsq = { FIP_DSQ                                 , $
                   MODEL_LABEL: string(' ',format='(A40)') , $
                   CHANSCAN: string(' ',format='(A4)')     , $
                   Q_INVERSE: dblarr( 5,5 )                , $
                   STRIPE_CONVERT: dblarr( 5,5 )               }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FIP_DSQ}, Nstruct )
end
