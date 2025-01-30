;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_COV_HIF2 creates and saves the binary records of HIGH_2 C-Matrix .
;
;  Written by :  Ken Jensen, Hughes STX, 23-May-1997
;-
;______________________________________________________________________________
;
Pro FMD_COV_HIF2,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_COV_HIF2 : Called Incorrectly : FMD_COV_HIF2,error'
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
PRINT,'Restoring HIGH_2_COVAR.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:high_2_covar.iss'
;

; Define the FMD_COV records
; --------------------------
nrec = N_ELEMENTS(freq_high_2)
cov_rec = FMD_COV_ST(nrec)
;

; Fill in the FMD_COV records
; ---------------------------
PRINT,' '
PRINT,'Filling in the FMD_COV records.'
PRINT,' '
FOR i=0,nrec-1 DO BEGIN
 cov_rec(i).row_num = i
 cov_rec(i).row_covr(0:nrec-1) = covar_high_2(i,*)
ENDFOR
;

; Write the FMD_COV file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_COV file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_cov_hif2.pass4'
rec_len = 730            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,cov_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_COV_HIF2.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & cvec_cut=0 & freq_dependence=0 & max_frac=0 & n_stripes=0
ndf_covar_high2=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
