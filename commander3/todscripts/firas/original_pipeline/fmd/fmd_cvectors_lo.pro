Pro FMD_CVECTORS_LO
;

;
;  FMD_CVECTORS_LO Makes the CVECTORS_LO Reference Data Set
;
;
;  Required Logicals :
;
;     CSDR$FIRAS_OUT  : Directory containing the xxx_CVECTOR.ISS .
;
;     CSDR$FIRAS_REF  : Directory where CVECTORS_LO.ISS will be sent.
;
;
;  Written by : Ken Jensen, Hughes STX, 22-Apr-1997
;
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
outtrans = STRUPCASE(outtrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,' '
;

; Restore the Destriped LLS C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:lls_cvector.iss'
f_lo = freq_lls
;

; Restore the Destriped RLS C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:rls_cvector.iss'
;
IF (MAX(ABS(f_lo-freq_rls)) NE 0.) THEN BEGIN
 PRINT,'FMD_CVECTORS_LO : LLS/RLS Frequency Mismatch !'
 RETURN
ENDIF
;

; Restore the Destriped LSF C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:lsf_cvector.iss'
;
IF (MAX(ABS(f_lo-freq_lsf)) NE 0.) THEN BEGIN
 PRINT,'FMD_CVECTORS_LO : LLS/LSF Frequency Mismatch !'
 RETURN
ENDIF
;

; Restore the Destriped RSF C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:rsf_cvector.iss'
;
IF (MAX(ABS(f_lo-freq_rsf)) NE 0.) THEN BEGIN
 PRINT,'FMD_CVECTORS_LO : LLS/RSF Frequency Mismatch !'
 RETURN
ENDIF
;

; Make the CVECTORS_LO Save Set
; -----------------------------
sname = 'csdr$firas_ref:cvectors_lo.iss'
SAVE,filename=sname,f_lo,cvec_lls,cvec_rls,cvec_lsf,cvec_rsf
;
PRINT,' '
PRINT,'IDL Save Set "' + reftrans(0) + 'CVECTORS_LO.ISS" Created.'
PRINT,' '
;

chanscan=0 & stripes_descrip=0 & freq_dependence=0 & cvec_cut=0
freq_lls=0 & freq_rls=0 & freq_lsf=0 & freq_rsf=0
ndf_cvec_lls=0 & ndf_cvec_rls=0 & ndf_cvec_lsf=0 & ndf_cvec_rsf=0
;

RETURN
END
