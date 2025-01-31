function FEL_ERR_str, Nstruct
 
  common FEL_ERR_str, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fel_err = { FEL_ERR                                 , $
                   PEP_GAIN: fltarr( 210 )                 , $
                   PEP_OFFSET: fltarr( 210 )               , $
                   JCJ_GAIN: fltarr( 13,210 )              , $
                   JCJ_OFFSET: fltarr( 13,210 )            , $
                   PUP_TEMP: 0.0                           , $
                   PUP_SPEC: fltarr( 210 )                 , $
                   PTP_TEMP: 0.0                           , $
                   PTP_SPEC: fltarr( 210 )                     }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FEL_ERR}, Nstruct )
end
