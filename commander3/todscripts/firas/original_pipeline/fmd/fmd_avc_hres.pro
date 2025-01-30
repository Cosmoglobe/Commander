;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_AVC_HRES creates and saves the binary records of HRES A-Vector .
;
;-
;______________________________________________________________________________
;
Pro FMD_AVC_HRES,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_AVC_HRES : Called Incorrectly : FMD_AVC_HRES,error'
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
PRINT,'Restoring HRES_AVECTOR.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:hres_avector.iss'
;

; Define the FMD_AVC record
; -------------------------
avc_rec = FMD_AVC_ST(1)
;

; Fill in the FMD_AVC record
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_AVC record.'
PRINT,' '
avc_rec.a_vector(0:31) = FLOAT(avec_hres(0:31))
;

; Write the FMD_AVC file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_AVC file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_avc_hres.pass4'
rec_len = 128            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
WRITEU,1,avc_rec(0)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_AVC_HRES.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & cvec_cut=0 & freq_dependence=0 & max_frac=0
freq_hres=0 & delta_freq_hres=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
