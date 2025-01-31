Function FMD_Pixel_Wgt,pixel=pixel,weights=weights,sky_idx=sky_idx, $
                       cvecs=cvecs,frac_wgt=frac_wgt,frac_fwgt=frac_fwgt, $
                       px_wgt=px_wgt
;
;
; FMD_PIXEL_WGT computes pixel weights and coadd fraction of pixel weight,
;               (frequency-dependent as well as frequency-independent).
;
;
; REQUIRED LOGICALS : None
;
;
; PROGRAMS Called   : None
;
;
; HISTORY : Written by Ken Jensen, Hughes STX, 12-May-1997
;
;
;

; Initialize Error Status
; -----------------------
error = 1
;

; Initialize Arrays
; -----------------
px_wgt = FLTARR(6144)        ; Pixel Weight
;
dummy = SIZE(cvecs)
nf = dummy(1)                ; Number of Frequencies
nch = dummy(2)               ; Number of CHANSCANs
nc = N_ELEMENTS(pixel)       ; Number of Coadds
;
frac_wgt = 0.*pixel          ; Coadd Fraction of Pixel Weight (Freq-Independent)
wfac = 0.*pixel + 1.         ; Weight Adjustment = 1. if Uncombined Data
;

IF (KEYWORD_SET(sky_idx)) THEN BEGIN
 ;
 ; Frequency Dependent Weight Adjustment
 ; -------------------------------------
 wfac = 0.*cvecs + 1.
 fac0 = cvecs(*,0)^2. * TOTAL(1./(cvecs(*,0)^2.))
 wfac(*,1) = fac0 / ( cvecs(*,1)^2. * TOTAL(1./(cvecs(*,1)^2.)) )
 IF (nch GT 2) THEN $
     wfac(*,2) = fac0 / ( cvecs(*,2)^2. * TOTAL(1./(cvecs(*,2)^2.)) )
 IF (nch GT 3) THEN $
     wfac(*,3) = fac0 / ( cvecs(*,3)^2. * TOTAL(1./(cvecs(*,3)^2.)) )
 ;
 frac_fwgt = FLTARR(nf,nc)    ; Coadd Fraction of Pixel Weight (Freq-Dependent)
 ;
ENDIF
;

; Good Coadds
; -----------
ng = WHERE(weights gt 0.,cg)
IF (cg le 0) THEN BEGIN
 PRINT,'FMD_PIXEL_WGT : No Good Coadds Found !'
 RETURN,error
ENDIF
;

; List of Pixels with Good Coadds
; -------------------------------
hpix = HISTOGRAM(pixel(ng),min=0,max=6143)
ng = WHERE(hpix GE 1,cg) 
IF (cg LT 1)THEN BEGIN
 PRINT,'FMD_PIXEL_WGT : No Pixels With Good Coadds !'
 RETURN,error
ENDIF
;

FOR jj=0,cg-1 DO BEGIN       ; Loop over Pixels with Good Coadds
;
 j = ng(jj)                   ; Pixel Number
 nx = WHERE(pixel eq j,cx)    ; Index of Coadds in Pixel
;
 IF (cx LE 0) THEN $
     PRINT,'FMD_PIXEL_WGT : Pixel' + STRCOMPRESS(STRING(j)) + ' : No Coadds'
;
 IF (cx GT 0) THEN BEGIN
 ;
  wix = weights(nx)            ; Destriper weights
  wp = TOTAL(wix)              ; Total Destriper Weight
  px_wgt(j) = wp               ; = Pixel Weight
  frac_wgt(nx) = wix / wp      ; Coadd Fraction of Weight in Pixel
 ;
  IF (KEYWORD_SET(cvecs)) THEN BEGIN
   ;
   idxx = sky_idx(nx)           ; Sky Index of Coadds in Pixel
   ;
   FOR k=0,nf-1 DO BEGIN        ; Loop over Frequencies
    wif = wix * wfac(k,idxx)           ; Frequency-Dependent Weight
    frac_fwgt(k,nx) = wif / TOTAL(wif) ; Coadd Fraction of Freq-Dependent Weight
   ENDFOR
   ;
  ENDIF
  ;
 ENDIF
 ;
ENDFOR
;

; No error
; --------
error = 0
;

RETURN,error
END
