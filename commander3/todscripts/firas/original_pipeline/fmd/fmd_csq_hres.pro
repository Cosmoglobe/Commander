;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_CSQ_HRES creates and saves the binary records of HRES pixel chi-squared.
;
;  Written by :  Ken Jensen,  Hughes STX,  03-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_CSQ_HRES,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_CSQ_HRES : Called Incorrectly : FMD_CSQ_HRES,error'
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
PRINT,'Restoring HRES_SKYMAP.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:hres_skymap.iss'
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
nf = N_ELEMENTS(freq_hres)
FOR i=0,nrec-1 DO BEGIN
 csq_rec(i).pixel = cmp_px(i)
 xpi = xp(i)
 ypi = yp(i)
 csq_rec(i).gal_lon = FLOAT(l_hres(xpi,ypi))
 csq_rec(i).gal_lat = FLOAT(b_hres(xpi,ypi))
 csq_rec(i).weight = FLOAT(n_hres(xpi,ypi))
 csq_rec(i).dof = LONG(nf) * nc_hres(xpi,ypi)
 csq_rec(i).comb_csq = FLOAT(chi2_hres(xpi,ypi))
ENDFOR
;

; Write the FMD_CSQ file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_CSQ file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_csq_hres.pass4'
rec_len = 24            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,csq_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_CSQ_HRES.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & stripes_descrip=0 & pcvr=0 & rect=0 & beta=0 & c_hres=0 & s_hres=0
stripe_contrib=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
