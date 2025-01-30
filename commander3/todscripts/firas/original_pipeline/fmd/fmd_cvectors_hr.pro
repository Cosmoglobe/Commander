Pro FMD_CVECTORS_HR
;

;
;  FMD_CVECTORS_HR Makes the CVECTORS_HR Reference Data Set
;
;
;  Required Logicals :
;
;     CSDR$FIRAS_OUT  : Directory containing the xxx_CVECTOR.ISS .
;
;     CSDR$FIRAS_REF  : Directory where CVECTORS_HR.ISS will be sent.
;
;
;  Written by : Ken Jensen, Hughes STX, 07-May-1997
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

; Restore the Destriped LLF C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:llf_cvector.iss'
f_hr = freq_llf
;

; Restore the Destriped RLF C-Vector
; ----------------------------------
RESTORE,'csdr$firas_out:rlf_cvector.iss'
;
IF (MAX(ABS(f_hr-freq_rlf)) NE 0.) THEN BEGIN
 PRINT,'FMD_CVECTORS_HR : LLF/RLF Frequency Mismatch !'
 RETURN
ENDIF
;

; Make the CVECTORS_HR Save Set
; -----------------------------
sname = 'csdr$firas_ref:cvectors_hr.iss'
SAVE,filename=sname,f_hr,cvec_llf,cvec_rlf
;
PRINT,' '
PRINT,'IDL Save Set "' + reftrans(0) + 'CVECTORS_HR.ISS" Created.'
PRINT,' '
;

chanscan=0 & stripes_descrip=0 & freq_dependence=0 & cvec_cut=0
freq_llf=0 & freq_rlf=0 & ndf_cvec_llf=0 & ndf_cvec_rlf=0
;

RETURN
END
