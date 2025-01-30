pro fer_pup,smode

;
;  Define the frequency array.
;
smode = strupcase(smode)
if (strmid(smode,1,1) eq 'H') then begin
   f = (4+dindgen(167))*144.981/256.0
endif else if (strmid(smode,1,2) eq 'LS') then begin
   f = (4+dindgen(34))*144.981/256.0
endif else if (strmid(smode,2,2) eq 'FL') then begin
   f = (4+dindgen(34))*144.981/256.0
endif else if (strmid(smode,1,2) eq 'LL') then begin
   f = (8+dindgen(148))*144.981/1024.0
endif

;
;  Restore the Ical temperature saveset.
;
restore,'csdr$firas_ref:fex_icaltemp.dat'

;
;  Compute the PUP error.
;
b = planck(tical,f,dbdt,units='icm',/mjy) * scale
dbdt = dbdt * scale
pup = dbdt * dt
save,file='csdr$firas_out:pup_errors.'+smode,dbdt,dt,f,pup

return
end
