function FDS_FEX_CVS, Nstruct
 
  common FDS_FEX_CVS, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fex_cvs = { FEX_CVS                                 , $
                   GMT: string(' ',format='(A14)')         , $
                   TIME: Lonarr( 2 )                       , $
                   CHANNEL: 0B                             , $
                   SCAN_LENGTH: 0B                         , $
                   SCAN_SPEED: 0B                          , $
                   MODEL_LABEL: string(' ',format='(A40)') , $
                   GALAT_EXC: 0.0                          , $
                   NCAL_IFGS: 0.0                          , $
                   CVECTOR: dblarr( 257 )                  , $
                   SPARES: bytarr( 47 )                        }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEX_CVS}, Nstruct )
end
