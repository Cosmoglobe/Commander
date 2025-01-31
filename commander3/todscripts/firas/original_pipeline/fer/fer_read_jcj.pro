pro fer_read_jcj,smode,pixel,model,f,jcj

;
; Read the JCJ file list
;
file = strarr(20)
n_files = 0
a = ''

openr,1,'csdr$firas_ref:fer_scr_jcj.txt'
while (not eof(1)) do begin
   readf,1,a
   file(n_files) = 'csdr$firas_in:fcf_dsk_' + smode + '.' + a
   n_files = n_files + 1
endwhile
close,1

file = file(0:n_files-1)
file = strupcase(file)

;
; Determine the calibration model solution and Nyquist frequency
;
field = 'spec_data.model_label,coad_spec_data.nyquist_icm'
stat = read_skymap(file(0),field,pixel,model,nyquist)
model = model(0)
nyquist = nyquist(0)

;
; Build the frequency vector
;
if (strupcase(strmid(smode,1,1)) eq 'H') then begin
   bins = '5:171'
   f = (4+findgen(167))*nyquist/256.0
endif else if (strupcase(strmid(smode,1,2)) eq 'LS') then begin
   bins = '5:38'
   f = (4+findgen(34))*nyquist/256.0
endif else if (strupcase(strmid(smode,2,2)) eq 'FL') then begin
   bins = '5:38'
   f = (4+findgen(34))*nyquist/256.0
endif else if (strupcase(strmid(smode,1,2)) eq 'LL') then begin
   bins = '9:156'
   f = (8+findgen(148))*nyquist/256.0
endif

;
;  Read the JCJ spectra files
;
field = 'spec_data.spec(' + bins + ')'
stat = read_skymap(file(0),field,pixel,spec0)
nsize = size(spec0)
jcj = complexarr(n_files,nsize(1),nsize(2))
jcj(0,*,*) = spec0

for j=1,n_files-1 do begin
   stat = read_skymap(file(j),field,pixel,spec)
   jcj(j,*,*) = spec - spec0
endfor


return
end
