;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_COV_HIGH creates and saves the binary records of HIGH C-Matrix .
;
;  Written by :  Ken Jensen, Hughes STX, 23-May-1997
;-
;______________________________________________________________________________
;
Pro FMD_COV_HIGH,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_COV_HIGH : Called Incorrectly : FMD_COV_HIGH,error'
 print,' '
 return
endif
;

; Logical Translations
; --------------------
PRINT,'Logical Translations :'
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
PRINT,' '
PRINT,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
PRINT,' '
;

; Restore the saveset
; -------------------
PRINT,'Restoring HIGH_COVAR.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:high_covar.iss'
;

; Define the FMD_COV records
; --------------------------
freq_high = [freq_high_2,freq_high_3,freq_high_4]
nrec = N_ELEMENTS(freq_high)
cov_rec = FMD_COV_ST(nrec)
;

; Fill in the FMD_COV records
; ---------------------------
PRINT,' '
PRINT,'Filling in the FMD_COV records.'
PRINT,' '
FOR i=0,nrec-1 DO BEGIN
 cov_rec(i).row_num = i
 cov_rec(i).row_covr(0:nrec-1) = covar_high(i,*)
ENDFOR
;

; Write the FMD_COV file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_COV file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_cov_high.pass4'
rec_len = 730            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,cov_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_COV_HIGH.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & cvec_cut=0 & freq_dependence=0 & max_frac=0 & n_stripes=0
ndf_covar_high=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
