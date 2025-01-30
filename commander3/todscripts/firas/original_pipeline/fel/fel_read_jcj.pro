pro fel_read_jcj,chanscan,pixel,f,jcj
;
; Read the JCJ file list
;
file = strarr(20)
n_files = 0
;
chanscan = strupcase(chanscan)
if (strmid(chanscan,1,1) eq 'H') then begin
   max_files = 14
endif else begin
   max_files = 12
endelse
while (n_files lt max_files) do begin
   n_str = strcompress(string(n_files),/remove_all)
   file(n_files) = 'csdr$firas_in:fsl_dsk_' + chanscan + '.ep_pass4' + n_str
   n_files = n_files + 1
endwhile
;
file = file(0:n_files-1)
file = strupcase(file)
;
; Build the frequency vector
;
if (strmid(chanscan,1,1) eq 'H') then begin
   f = (5+findgen(210))*144.981*1.00159/320.0
endif else if (strmid(chanscan,1,2) eq 'LL') then begin
   f = (7+findgen(182))*36.245*1.00159/320.0
endif else begin
   f = (5+findgen(43))*144.981*1.00159/320.0
endelse
;
;  Read the JCJ spectra files
;
openr,lunr,file(0),/get_lun,/share
finfo=fstat(lunr)
num_recs=finfo.size/finfo.rec_len
r=assoc(lunr,fel_fsl_dsk_st(num_recs),0)
fsl_dsk = r(0)
close,lunr & free_lun,lunr
if (strmid(chanscan,1,1) eq 'H') then begin
   spec0 = fsl_dsk.spec_data.spec(5:214)
endif else if (strmid(chanscan,1,2) eq 'LL') then begin
   spec0 = fsl_dsk.spec_data.spec(7:188)
endif else begin
   spec0 = fsl_dsk.spec_data.spec(5:47)
endelse
pixel = fsl_dsk.attitude.pixel_no
nsize = size(spec0)
jcj = complexarr(n_files,nsize(1),nsize(2))
jcj(0,*,*) = spec0
;
for j=1,n_files-1 do begin
   openr,lunr,file(j),/get_lun,/share
   finfo=fstat(lunr)
   num_recs=finfo.size/finfo.rec_len
   r=assoc(lunr,fel_fsl_dsk_st(num_recs),0)
   fsl_dsk = r(0)
   close,lunr & free_lun,lunr
   if (strmid(chanscan,1,1) eq 'H') then begin
      spec = fsl_dsk.spec_data.spec(5:214)
   endif else if (strmid(chanscan,1,2) eq 'LL') then begin
      spec = fsl_dsk.spec_data.spec(7:188)
   endif else begin
      spec = fsl_dsk.spec_data.spec(5:47)
   endelse
   jcj(j,*,*) = spec - spec0
endfor

return
end
