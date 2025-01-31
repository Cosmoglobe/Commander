function FDS_FIP_DSE, Nstruct
 
  common FDS_FIP_DSE, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fip_dse = { FIP_DSE                                 , $
                   PIXEL: 0L                               , $
                   MODEL_LABEL: string(' ',format='(A40)') , $
                   CHANSCAN: string(' ',format='(A4)')     , $
                   DIAG_ELEMENT: 0.D0                      , $
                   RECT_ELEMENT: dblarr( 5 )               , $
                   BETA_ELEMENT: dblarr( 5 )               , $
                   STRIPE_CONTRIB: 0L                          }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FIP_DSE}, Nstruct )
end
