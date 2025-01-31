FUNCTION FMD_CVECTOR, pixel=pixel,spec=spec,sky_wgts_ds=sky_wgts_ds, $
                      sky_idx=sky_idx,dmask=dmask,mask=mask,wfac=wfac, $
                      w_frac=w_frac,wf_frac=wf_frac,n_stripes=n_stripes, $
                      max_frac=max_frac,good_dvec=good_dvec,ndf_cvec=ndf_cvec, $
                      good_cvec=good_cvec
;
;
; FMD_CVECTOR computes the C-Vector .
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
;  DMASK (I)        : = 0 if Coadd is masked by DIRBE_CUT
;
;  MASK (I)         : = 0 if Coadd is masked by CVEC_CUT or LATCUT/LONCUT
;
;  WFAC (I)         : Frequency-Dependent Weight Factor (Optional)
;
;  W_FRAC (I)       : Coadd Fraction of Weight in Pixel 
;
;  WF_FRAC (I)      : Coadd Fraction of Frequency-Dependent
;                     Weight in Pixel (Optional)
;
;  N_STRIPES (I)    : Number of Stripes 
;
;  MAX_FRAC (I)     : Maximum Fraction of Pixel Weight Allowed for
;                     C-Vector calculation
;
;  GOOD_DVEC (I)    : Index of Coadds used in D-Vector calculation
;
;  NDF_CVEC (O)     : Degrees of Freedom in C-Vector calculation
;
;  GOOD_CVEC (I)    : Index of Coadds used in C-Vector calculation
;
;
; RETURNS :  CVEC  =  C-Vector (MJy/sr)
;
;
; REQUIRED LOGICALS : None
;
;
; PROGRAMS Called   : FMD_MATCH
;
;
; HISTORY : Written by Ken Jensen, Hughes STX, 31-May-97
;
;
;

; Initialize Arrays
; -----------------
;
nc = N_ELEMENTS(pixel)
nf = N_ELEMENTS(spec) / nc
;
cvec = FLTARR(nf)           ; C-Vector (new calculation)
pspec = FLTARR(nf,6144)     ; Pixel Spectrum
resid = FLTARR(nf,nc)       ; Coadd Residual
chi2 = FLTARR(nf,nc)        ; Coadd Chi-Squared
;

; Good Coadds
; -----------
good = WHERE((sky_wgts_ds GT 0.)and(dmask GT 0),cgood)
IF (cgood LE 0) THEN BEGIN
 PRINT,'FMD_CVECTOR : No Good Coadds Found !'
 RETURN,cvec
ENDIF
;

hpix = HISTOGRAM(pixel(good),min=0,max=6143)
;

ng = WHERE(hpix GE 1,cg)  ;  List of Pixels with Good Coadds
IF (cg LT 1) THEN BEGIN
 PRINT,'FMD_CVECTOR : No Pixels With Good Coadds !'
 RETURN,cvec
ENDIF
;

spix = STRCOMPRESS(STRING(cg))
PRINT,' '
PRINT,'Spectra will be Computed for' + spix + ' Pixels.' 
PRINT,' '
;

ng2 = WHERE(hpix GT 1,cg2)  ;  List of Pixels with Multiple Good Coadds
IF (cg2 LT 1)THEN BEGIN
 PRINT,'FMD_CVECTOR : No Pixels With Multiple Good Coadds !'
 RETURN,cvec
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
 wix = sky_wgts_ds(nx) * dmask(nx) ; Destriper weights * DIRBE mask
 nw = WHERE(wix GT 0.,cw)          ; CW = Number of Good Coadds in Pixel
 ;
 IF (cw LE 0) THEN BEGIN
  sj = STRCOMPRESS(STRING(j))
  PRINT,'FMD_CVECTOR : Good Coadds Not Found for Pixel' + sj + ' !!'
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

; C-VECTOR (New Calculation)
; --------------------------
;
; Good Coadds for C-Vector
; ------------------------
rt = TOTAL(resid,1)
good2 = $
  WHERE((rt NE 0.)and(w_frac GT 0.)and(w_frac LT max_frac)and(mask GT 0),cgood2)
;
; IF No Coadds, return with error
; -------------------------------
IF (cgood2 le 0) THEN BEGIN
 PRINT,'FMD_CVECTOR : No Coadds for C-Vector Calculation !'
 RETURN,cvec
ENDIF
;

; Match to Good_DVEC Coadds
; -------------------------
FMD_MATCH,good2,good_dvec,good3,dummy,'Y',error
good_cvec = good2(good3)
;

hg = HISTOGRAM(pixel(good_cvec),min=0,max=6143)
ngg = WHERE(hg GT 0,cgg)
sg = STRCOMPRESS(STRING(N_ELEMENTS(good_cvec)))
sgg = STRCOMPRESS(STRING(cgg))
PRINT,' '
PRINT,'C-Vector will be Computed from'+sg+' Coadds in'+sgg+' Pixels.' 
PRINT,' '
;

; Frequency-Independent Weight for C-Vector
; -----------------------------------------
freq_flag = 'N'
;

; IF ((SKY_IDX)and(WFAC)and(WF_FRAC)), Frequency-Dependent Weight for C-Vector
; ----------------------------------------------------------------------------
IF ((KEYWORD_SET(sky_idx))and(KEYWORD_SET(wfac))and(KEYWORD_SET(wf_frac))) $
     THEN freq_flag = 'Y'
;

; Frequency-Independent Weights
; -----------------------------
norm = sky_wgts_ds(good_cvec) / (1. - w_frac(good_cvec))
;

FOR i=0,nf-1 DO BEGIN  ; Loop over Frequencies
 ;
 ; IF (FREQ_FLAG), Frequency-Dependent Weights
 ; -------------------------------------------
 IF (freq_flag EQ 'Y') THEN BEGIN
  sky_wgts_n = sky_wgts_ds * REFORM(wfac(i,sky_idx))
  norm(*) = sky_wgts_n(good_cvec) / (1. - REFORM(wf_frac(i,good_cvec)))
 ENDIF
 ;
 ; C-Variance
 ; ----------
 cvec(i) = TOTAL( (resid(i,good_cvec)^2.) * norm )
 ;
ENDFOR
;

ndf_cvec = N_ELEMENTS(good_cvec) - n_stripes   ;  Degrees of Freedom for C-Vector
;
cvec = SQRT(cvec/ndf_cvec)      ;  C-Vector
;

RETURN,cvec
END
