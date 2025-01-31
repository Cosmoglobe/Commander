function FMD_COV_ST, Nstruct
 
  common FMD_COV_ST, defined
 
  if N_elements( defined ) NE 1 then begin
 
       fmd_cov = { FMD_COV                     , $
                   ROW_NUM : 0                 , $
                   ROW_COVR: fltarr(182)           }
 
       defined = 1 
    endif
 
  if N_elements( Nstruct ) NE 1 then Nstruct = 1
 
return, replicate( {FMD_COV}, Nstruct )
end
