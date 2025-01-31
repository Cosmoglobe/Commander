;_______________________________________________________________________________
;
;+NAME/ONE-LINE DESCRIPTION:
;    FRD_CAL_WEIGHTS computes calibration model solution relative weights
;
; DESCRIPTION:
;    FRD_CAL_WEIGHTS computes calibration model solution relative weights from
;    the C-Vector and the number of cold-null cal ifgs that went into the
;    calibration model solution.  The C-Vector and the ncal ifgs are read from
;    the FEX_CVS files.
;
; CALLING SEQUENCE:
;    cal_wgt = FRD_CAL_WEIGHTS (infile, smode)
;
; ARGUMENTS (I= input, O=output)
;    infile  I  string  FEX_CVS file name
;    smode   I  string  Scan mode for FRD_MERGE_CSP
;
; WARNINGS:
;    The valid scan modes are LOWF, HIGH, LOW2, HIG2, HRES, and LRES.
;#
; COMMON BLOCKS:  none
;
; PROCEDURE:
;       (1)  Read the FEX_CVS file.
;       (2)  Compute the cal weights
;
; LIBRARY CALLS:
;    Calls FRD_FEX_CVS to define the FEX_CVS record strucutre.
;
; MODIFICATION HISTORY:
;    Written by Gene Eplee, General Sciences Corp., 27 October 1994
;    Modified by Ken Jensen, Hughes STX, 1 Dec 1994, LOW2 and HIG2 modes added.
;-
;______________________________________________________________________________
;
function frd_cal_weights,file,smode
;
; Read the CVS file.
;
openr,lunr,file,/get_lun,/share
finfo=fstat(lunr)
nrecs=finfo.size/finfo.rec_len
;
r=assoc(lunr,frd_fex_cvs(nrecs))
cvs_rec = r(0)
;
close,lunr
free_lun,lunr
;
; Read the C-Vector^2 from the CVS record.
;
if ((smode eq 'LOWF')or(smode eq 'LOW2')or(smode eq 'LRES')) then begin
   jstart = 4
   nfreq = 34
endif else if ((smode eq 'HIGH')or(smode eq 'HIG2')) then begin
   jstart = 4
   nfreq = 167
endif else if (smode eq 'HRES') then begin
   jstart = 8
   nfreq = 148
endif
cvs = dblarr(nfreq)
for  j=0,nfreq-1 do cvs(j) = cvs_rec.cvector(j+jstart)
;
; Compute the cal weight
;
cal_wgt = TOTAL(cvs_rec.ncal_ifgs/cvs)

return,cal_wgt
end
