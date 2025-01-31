;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_CVC_LOWF creates and saves the binary records of LOWF C-Vector .
;
;-
;______________________________________________________________________________
;
Pro FMD_CVC_LOWF,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_CVC_LOWF : Called Incorrectly : FMD_CVC_LOWF,error'
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
PRINT,'Restoring LOWF_ERRORS.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:lowf_errors.iss'
;

; Define the FMD_CVC record
; -------------------------
cvc_rec = FMD_CVC_ST(1)
;

; Fill in the FMD_CVC record
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_CVC record.'
PRINT,' '
nf = N_ELEMENTS(c_lowf)
cvc_rec.c_vector(0:nf-1) = FLOAT(c_lowf)
IF (nf LT 182) THEN cvc_rec.c_vector(nf:181) = 0.
;

; Write the FMD_CVC file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_CVC file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_cvc_lowf.pass4'
rec_len = 728            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
WRITEU,1,cvc_rec(0)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_CVC_LOWF.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
stripes_descrip=0 & cmp_px=0 & pcvr=0 & rect=0 & diag=0 & beta=0 & omega=0
square=0 & d_inv=0 & stripe_conv=0 & destriper_wgt=0 & stripe_contrib=0
lmat=0 & chanscan=0 & f_lo=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
