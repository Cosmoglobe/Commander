;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_CSQ_HIF2 creates and saves the binary records of HIGH_2 pixel chi-squared.
;
;  Written by :  Ken Jensen,  Hughes STX,  03-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_CSQ_HIF2,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_CSQ_HIF2 : Called Incorrectly : FMD_CSQ_HIF2,error'
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
PRINT,'Restoring HIGH_2_SKYMAP.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:high_2_skymap.iss'
;

; Define the FMD_CSQ records
; --------------------------
nrec = N_ELEMENTS(cmp_px)
csq_rec = FMD_CSQ_ST(nrec)
;

; Skymap Indices
; --------------
PIX2XY,cmp_px,xp,yp,res=6,/six
;

; Fill in the FMD_CSQ records
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_CSQ records.'
PRINT,' '
nf = N_ELEMENTS(freq_high_2)
FOR i=0,nrec-1 DO BEGIN
 csq_rec(i).pixel = cmp_px(i)
 xpi = xp(i)
 ypi = yp(i)
 csq_rec(i).gal_lon = FLOAT(l_high(xpi,ypi))
 csq_rec(i).gal_lat = FLOAT(b_high(xpi,ypi))
 csq_rec(i).weight = FLOAT(n_high(xpi,ypi))
 csq_rec(i).dof = LONG(nf) * nc_high(xpi,ypi)
 csq_rec(i).comb_csq = FLOAT(chi2_high_2(xpi,ypi))
ENDFOR
;

; Write the FMD_CSQ file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_CSQ file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_csq_hif2.pass4'
rec_len = 24            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,csq_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_CSQ_HIF2.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & stripes_descrip=0 & pcvr=0 & rect=0 & beta=0
c_high_2=0 & s_high_2=0 & z_high_2=0 & stripe_contrib=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
