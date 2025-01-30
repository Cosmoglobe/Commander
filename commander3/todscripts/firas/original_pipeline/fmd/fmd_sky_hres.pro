;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_SKY_HRES creates and saves the binary records of destriped calibrated
;               HRES sky spectra .
;
;  Written by :  Ken Jensen,  Hughes STX,  04-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_Sky_HRES,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_SKY_HRES : Called Incorrectly : FMD_SKY_HRES,error'
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

; Define the FMD_SKY records
; --------------------------
nrec = N_ELEMENTS(cmp_px)
sky_rec = FMD_SKY_ST(nrec)
;

; Skymap Indices
; --------------
PIX2XY,cmp_px,xp,yp,res=6,/six
;

; Fill in the FMD_SKY records
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_SKY records.'
PRINT,' '
nf = N_ELEMENTS(freq_hres)
FOR i=0,nrec-1 DO BEGIN
 sky_rec(i).pixel = cmp_px(i)
 xpi = xp(i)
 ypi = yp(i)
 sky_rec(i).gal_lon = FLOAT(l_hres(xpi,ypi))
 sky_rec(i).gal_lat = FLOAT(b_hres(xpi,ypi))
 sky_rec(i).weight = FLOAT(n_hres(xpi,ypi))
 sky_rec(i).str_used = BYTE(stripe_contrib(i)) 
 sky_rec(i).spectrum(0:nf-1) = FLOAT(s_hres(xpi,ypi,*))
ENDFOR
;

; Write the FMD_SKY file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_SKY file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_sky_hres.pass4'
rec_len = 745            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,sky_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_SKY_HRES.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
chanscan=0 & stripes_descrip=0 & pcvr=0 & rect=0 & beta=0
c_hres=0 & nc_hres=0 & chi2_hres=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
