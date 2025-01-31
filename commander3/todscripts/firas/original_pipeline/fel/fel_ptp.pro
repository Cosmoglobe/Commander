pro fel_ptp,chanscan
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
;  Restore the CMBR temperature long version saveset.
;
restore,'csdr$firas_ref:fex_cmbrtempl.dat'
;
;  Compute the PTP error.
;
b = planck(tcmbr,f,dbdt,units='icm',/mjy)
dbdt = dbdt * scale
ptp = dbdt * dt
ptp = ptp/scale
save,file='csdr$firas_out:ptp_errors.'+chanscan,dbdt,dt,f,ptp

return
end
