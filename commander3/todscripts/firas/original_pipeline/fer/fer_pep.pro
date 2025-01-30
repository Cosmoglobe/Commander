pro fer_pep,smode

;
;  Read the data from the ASCII file.
;
readcol,'csdr$firas_in:emis.'+smode,fr,a1,a2,a3,a4,a5,a6,pep,gain,rchi,ichi

;
;  Define the frequency array and extract the PEP term.
;
if (strupcase(strmid(smode,1,1)) eq 'H') then begin
   f = (4+findgen(167))*144.981/256.0
   gain = gain(0:166)
   pep = pep(0:166)
endif else if (strupcase(strmid(smode,1,2)) eq 'LS') then begin
   f = (4+findgen(34))*144.981/256.0
   gain = gain(0:33)
   pep = pep(0:33)
endif else if (strupcase(strmid(smode,2,2)) eq 'FL') then begin
   f = (4+findgen(34))*144.981/256.0
   gain = gain(0:33)
   pep = pep(0:33)
endif else if (strupcase(strmid(smode,1,2)) eq 'LL') then begin
   f = (8+findgen(148))*144.981/1024.0
   gain = gain(4:151)
   pep = pep(4:151)
endif

;
;  Save the PEP and gain error terms.
;
save,file='csdr$firas_out:pep_errors.'+smode,f,gain,pep

return
end
