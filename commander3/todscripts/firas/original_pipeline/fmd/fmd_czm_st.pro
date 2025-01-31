function FMD_CZM_ST, Nstruct
 
  common FMD_CZM_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_czm = { FMD_CZM                              , $
                   PIXEL: 0L                            , $
                   CHANSCAN: string(' ',format='(A4)')  , $
                   REC_NUM: 0L                          , $
                   WEIGHT: 0.0                          , $
                   ZODI_MOD: fltarr( 170 )                  }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_CZM}, Nstruct )
end
