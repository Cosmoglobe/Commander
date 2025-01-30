FUNCTION FMD_RESID, pixel=pixel,spec=spec,sky_wgts_ds=sky_wgts_ds, $
                    sky_idx=sky_idx,mask=mask,w_frac=w_frac, $
                    wf_frac=wf_frac,px_spec=px_spec,resid=resid
;
;
; FMD_RESID computes weighted mean sky pixel spectra and coadd residuals.
;
;
; KEYWORDS   :
;
;  PIXEL (I)        : Coadd Pixel Number
;
;  SPEC (I)         : Coadd Spectra (MJy/sr)
;
;  SKY_WGTS_DS (I)  : Coadd Destriper Weights
;
;  SKY_IDX (I)      : Channel/Scan-Mode Index for each Coadd (Optional)
;
;  MASK (I)         : = 0 if Coadd is masked.
;
;  W_FRAC (I)       : Coadd Fraction of Weight in Pixel 
;
;  WF_FRAC (I)      : Coadd Fraction of Frequency-Dependent
;                     Weight in Pixel (Optional)
;
;  PX_SPEC (O)      : Pixel Spectra (MJy/sr)
;
;  RESID (O)        : Coadd Residuals (MJy/sr)
;
;
; RETURNS :  ERROR  =  Return Status
;
;
; REQUIRED LOGICALS : None
;
;
; PROGRAMS Called   : None
;
;
; HISTORY : Written by Ken Jensen, Hughes STX, 10-Apr-97
;
;
;

; Initialize Return Status
; ------------------------
error = 1
;

; Initialize Arrays
; -----------------
nc = N_ELEMENTS(pixel)
nf = N_ELEMENTS(spec) / nc
pspec = FLTARR(nf,6144)     ; Pixel Spectrum
resid = FLTARR(nf,nc)       ; Coadd Residual
;

; Good Coadds
; -----------
good = WHERE((sky_wgts_ds GT 0.)and(mask GT 0),cgood)
IF (cgood LE 0) THEN BEGIN
 PRINT,'FMD_RESID : No Good Coadds Found !'
 RETURN,error
ENDIF
;

hpix = HISTOGRAM(pixel(good),min=0,max=6143)
;

ng = WHERE(hpix GE 1,cg)  ;  List of Pixels with Good Coadds
IF (cg LT 1) THEN BEGIN
 PRINT,'FMD_RESID : No Pixels With Good Coadds !'
 RETURN,error
ENDIF
;

spix = STRCOMPRESS(STRING(cg))
PRINT,' '
PRINT,'Spectra will be Computed for' + spix + ' Pixels.' 
PRINT,' '
;

ng2 = WHERE(hpix GT 1,cg2)  ;  List of Pixels with Multiple Good Coadds
IF (cg2 LT 1)THEN BEGIN
 PRINT,'FMD_RESID : No Pixels With Multiple Good Coadds !'
 RETURN,error
ENDIF
;

spix = STRCOMPRESS(STRING(cg2))
scoad = STRCOMPRESS(STRING(LONG(TOTAL(hpix(ng2)))))
PRINT,' '
PRINT,'Residuals will be Computed for'+scoad+' Coadds in'+spix+' Pixels.' 
PRINT,' '
;

; PIXEL SPECTRA  and  COADD RESIDUALS
; -----------------------------------
;
FOR jj=0,cg-1 DO BEGIN       ; Loop over Pixels with Good Coadds
 ;
 j = ng(jj)                        ; Pixel Number
 nx = WHERE(pixel EQ j,cx)         ; Index of Coadds in Pixel
 wix = sky_wgts_ds(nx) * mask(nx) ; Destriper weights * DIRBE mask
 nw = WHERE(wix GT 0.,cw)          ; CW = Number of Good Coadds in Pixel
 ;
 IF (cw LE 0) THEN BEGIN
  sj = STRCOMPRESS(STRING(j))
  PRINT,'FMD_RESID : Good Coadds Not Found for Pixel' + sj + ' !!'
  RETURN,error
 ENDIF
 ;
 IF (cw EQ 1) THEN pspec(*,j) = spec(*,nx(nw(0)))
 ;
 IF (cw GT 1) THEN BEGIN

  ; Pixel Spectrum
  ; --------------
  ;
  ; IF (WF_FRAC), use Frequency-Dependent Weights, ELSE use Freq-Independent
  ; ------------------------------------------------------------------------
  IF (KEYWORD_SET(wf_frac)) THEN BEGIN
   FOR k=0,nf-1 DO pspec(k,j) = TOTAL(spec(k,nx) * wf_frac(k,nx))
  ENDIF ELSE BEGIN
   FOR k=0,nf-1 DO pspec(k,j) = TOTAL(spec(k,nx) * w_frac(nx))
  ENDELSE
  ;

  ; Coadd Residuals
  ; ---------------
  FOR k=0,cx-1 DO  resid(*,nx(k)) = spec(*,nx(k)) - pspec(*,j)
  ;

 ENDIF
 ;
ENDFOR
;

; SKYMAPS
; -------
PIX2XY,INDGEN(6144),xp,yp,data=pspec,raster=px_spec,res=6,/six
;

error = 0
;

RETURN,error
END
