pro fel_pup,chanscan
;
;  Define the frequency array.
;
chanscan = strupcase(chanscan)
if (strmid(chanscan,1,1) eq 'H') then begin
   f = (5+dindgen(210))*144.981*1.00159/320.0
endif else if (strmid(chanscan,1,2) eq 'LL') then begin
   f = (7+dindgen(182))*36.245*1.00159/320.0
endif else begin
   f = (5+dindgen(43))*144.981*1.00159/320.0
endelse
;
;  Restore the Ical temperature saveset.
;
restore,'csdr$firas_ref:fex_icaltemp.dat'
;
;  Compute the PUP error.
;
b = planck(tical,f,dbdt,units='icm',/mjy)
dbdt = dbdt * scale
pup = dbdt * dt
pup = pup/scale
save,file='csdr$firas_out:pup_errors.'+chanscan,dbdt,dt,f,pup

return
end
