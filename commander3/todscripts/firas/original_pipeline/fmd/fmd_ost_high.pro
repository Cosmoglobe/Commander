;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_OST_HIGH creates and saves the binary records of HIGH stripe arrays
;               in the orthogonal basis.
;
;  Written by :  Ken Jensen,  Hughes STX,  09-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_OST_HIGH,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_OST_HIGH : Called Incorrectly : FMD_OST_HIGH,error'
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
PRINT,'Restoring HIGH_ERRORS.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:high_errors.iss'
;

; Define the FMD_OST records
; --------------------------
sz = SIZE(beta)
nrec = sz(2)
ost_rec = FMD_OST_ST(nrec)
;

; Fill in the FMD_OST records
; ---------------------------
PRINT,' '
PRINT,'Filling in the FMD_OST records.'
PRINT,' '
sz = SIZE(gamma_high)
nf = sz(1)
FOR i=0,nrec-1 DO BEGIN
 ost_rec(i).str_gama(0:nf-1) = gamma_high(*,i)
 ost_rec(i).str_beta(cmp_px) = FLOAT(beta(*,i))
ENDFOR
;

; Write the FMD_OST file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_OST file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_ost_high.pass4'
rec_len = 25304             ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,ost_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_OST_HIGH.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & diag=0 & stripe_contrib=0 & c_high=0 & f_hi=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
