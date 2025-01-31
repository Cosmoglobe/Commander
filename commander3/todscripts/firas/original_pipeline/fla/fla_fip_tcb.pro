function FLA_FIP_TCB, Nstruct
 
  common FLA_FIP_TCB, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fip_tcb = { FIP_TCB                             , $
                   PIXEL: 0L                           , $
                   ECLON: 0.0                          , $
                   ECLAT: 0.0                          , $
                   TEMP: 0.0                           , $
                   TEMP_SIG: 0.0                       , $
                   RESID_TEMP: 0.0                     , $
                   CHANSCAN: string(' ',format='(A4)') , $
                   NUM_IFGS: 0.0                       , $
                   GALON: 0.0                          , $
                   GALAT: 0.0                          , $
                   RA: 0.0                             , $
                   DEC: 0.0                            , $
                   GALATEXC: 0.0                           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FIP_TCB}, Nstruct )
end
