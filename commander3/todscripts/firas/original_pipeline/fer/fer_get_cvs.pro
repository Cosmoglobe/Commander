pro fer_get_cvs,chanscan,solution,cvs_rec
;
; IDL procedure to read the FEX_CVS binary files for FER_ERRORS.
;
;  Author:  Gene Eplee
;           General Sciences Corp.
;           9 August 1994
;
file = 'csdr$firas_cvs:fex_cvs_' + chanscan + '.' + solution
openr,lunr,file,/get_lun,/share
finfo=fstat(lunr)
nrecs=finfo.size/finfo.rec_len
;
r=assoc(lunr,fer_fex_cvs(nrecs))
cvs_rec = r(0)
;
close,lunr
free_lun,lunr
;
return
end
