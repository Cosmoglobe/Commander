function FMD_PST_ST, Nstruct
 
  common FMD_PST_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_pst = { FMD_PST                            , $
                   STR_ID: string(' ',format='(A40)') , $
                   STR_SPEC: fltarr( 182 )            , $
                   STR_RECT: fltarr( 6144 )           , $
                   STR_COVR: fltarr( 29 )                 }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_PST}, Nstruct )
end
