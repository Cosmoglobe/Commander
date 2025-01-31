function FER_ERR_str, Nstruct
 
  common FER_ERR_str, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fer_err = { FER_ERR                                 , $
                   GMT: string(' ',format='(A14)')         , $
                   TIME: Lonarr( 2 )                       , $
                   MODEL_LABEL: string(' ',format='(A40)') , $
                   CHANNEL: 0L                             , $
                   SCAN_LENGTH: 0L                         , $
                   SCAN_SPEED: 0L                          , $
                   GALATEXC: 0.0                           , $
                   D_VECTOR: dblarr( 200 )                 , $
                   C_VECTOR: dblarr( 200 )                 , $
                   PEP_GAIN: fltarr( 200 )                 , $
                   PEP_OFFSET: fltarr( 200 )               , $
                   JCJ_GAIN: fltarr( 14,200 )              , $
                   JCJ_OFFSET: fltarr( 14,200 )            , $
                   PUP_TEMP: 0.0                           , $
                   PUP_SPEC: fltarr( 200 )                 , $
                   PTP_TEMP: 0.0                           , $
                   PTP_SPEC: fltarr( 200 )                 , $
                   SPARES: bytarr( 298 )                       }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FER_ERR}, Nstruct )
end
