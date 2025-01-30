pro fel_pep,chanscan
;
;  Read the data from the ASCII file.
;
readcol,'csdr$firas_in:emis.'+chanscan,fr,a1,a2,a3,a4,a5,a6,pep,gain,rchi,ichi
;
;  Define the frequency array and extract the PEP term.
;
if (strupcase(strmid(chanscan,1,1)) eq 'H') then begin
   f = (5+findgen(210))*144.981*1.00159/320.0
   gain = gain(0:209)
   pep = pep(0:209)
endif else if (strupcase(strmid(chanscan,1,2)) eq 'LL') then begin
   f = (7+findgen(182))*36.245*1.00159/320.0
   gain = gain(0:181)
   pep = pep(0:181)
endif else begin
   f = (5+findgen(43))*144.981*1.00159/320.0
   gain = gain(0:42)
   pep = pep(0:42)
endelse
;
;  Save the PEP and gain error terms.
;
pep = pep/2.99792458e-7
save,file='csdr$firas_out:pep_errors.'+chanscan,f,gain,pep

return
end
