pro fer_get_var,chanscan,solution,var_rec
;
; IDL procedure to read the FEX_VAR binary files for FER_ERRORS.
; For the merged short fast and long fast scan modes, the FMS_VAR binary
; files are read.
;
;  Author:  Gene Eplee
;           General Sciences Corp.
;           18 August 1994
;

;
;  Get the correct file name for the merged short fast and long fast datasets.
;
file = 'csdr$firas_var:fex_var_' + chanscan + '.' + solution
if (chanscan eq 'RHLF') then begin
   file = 'csdr$firas_var:fms_var_rhfa' + '.' + solution
endif
if (chanscan eq 'RLSF') then begin
   file = 'csdr$firas_var:fms_var_rlfa' + '.' + solution
endif
if (chanscan eq 'LHLF') then begin
   file = 'csdr$firas_var:fms_var_lhfa' + '.' + solution
endif
if (chanscan eq 'LLSF') then begin
   file = 'csdr$firas_var:fms_var_llfa' + '.' + solution
endif
;
openr,lunr,file,/get_lun,/share
finfo=fstat(lunr)
nrecs=finfo.size/finfo.rec_len
;
r=assoc(lunr,fer_fex_var(nrecs))
var_rec = r(0)
;
close,lunr
free_lun,lunr
;
return
end
