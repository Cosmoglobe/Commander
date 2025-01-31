FUNCTION FMD_COVAR,nfreq=nfreq,resid=resid,weight=weight,mask=mask,px=px,$
                   cvec_cut=cvec_cut,n_stripes=n_stripes,max_frac=max_frac,$
                   wfac=wfac,w_frac=w_frac,wf_frac=wf_frac,sky_idx=sky_idx,$
                   covar=covar,ndf_covar=ndf_covar
;

;
;  FMD_COVAR computes the covariance of destriped coadd spectra.
;
;
;  PROGRAMS Called       :  COORCONV
;
;
;  HISTORY : Written by Ken Jensen, Hughes STX, 11-Apr-1997.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;; BEGIN
;
; Initialize Return Error
; -----------------------
error = 1
;

; Number of Frequencies
; ---------------------
n_freq = N_ELEMENTS(resid) / N_ELEMENTS(weight)
;

; Sanity Check
; ------------
IF (n_freq NE nfreq) THEN BEGIN
 PRINT,'FMD_COVAR : Error in Array Dimensions !'
 RETURN,error
ENDIF
;

; C-Matrix Mask (from CVEC_CUT)
; -----------------------------
cmask = mask
good_cut = 'Y'
IF ((cvec_cut(0) LE 0)or(cvec_cut(0) GT 90)) THEN good_cut = 'N'
IF ((cvec_cut(1) LE 0)or(cvec_cut(1) GT 180)) THEN good_cut = 'N'
;
IF (good_cut eq 'Y') THEN BEGIN
 ;
 ll = COORCONV(px,infmt='p',inco='f',outfmt='l',outco='g')
 pix_glat = ll(*,1)
 pix_glon = ll(*,0)
 nq = WHERE(pix_glon gt 180.,cq)
 IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) - 360.
 nq = WHERE(pix_glon lt -180.,cq)
 IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) + 360.
  bad_cvec = $
    WHERE((ABS(pix_glat) LE cvec_cut(0))and(ABS(pix_glon) LE cvec_cut(1)),cbad)
 IF (cbad GT 0) THEN cmask(bad_cvec) = 0
 ;
ENDIF
;

; Good Coadds for C-Matrix
; ------------------------
rt = TOTAL(resid,1)  ; Total Residual per Coadd
;
ng = WHERE((rt NE 0.)and(w_frac GT 0.)and(w_frac LT max_frac)and(cmask EQ 1),cg)
IF (cg LE 0) THEN BEGIN
 PRINT,'FMD_COVAR : No Good Coadds'
 RETURN,error
ENDIF
;

; Good Pixels for C-Matrix
; ------------------------
pxg = px(ng)
hpix = HISTOGRAM(pxg,min=0)
npix = WHERE(hpix GT 1,cpix)
IF (cpix LE 0) THEN BEGIN
 PRINT,'FMD_COVAR : No Good Pixels with Multiple Coadds'
 RETURN,error
ENDIF
;

; C-MATRIX
; --------
;
sg = STRCOMPRESS(STRING(cg))
spix = STRCOMPRESS(STRING(cpix))
PRINT,' '
PRINT,'C-Matrix will be Computed from'+sg+' Coadds in'+spix+' Pixels.' 
PRINT,' '
;

rg = FLTARR(n_freq,cg)  ;  Initialize Residuals matrix
;

; Frequency-Independent Weighting
; -------------------------------
IF (NOT KEYWORD_SET(wfac)) THEN FOR i=0,n_freq-1 DO $
    rg(i,*) = resid(i,ng(*)) * SQRT(weight(ng)/(1. - w_frac(ng)))
;

; Frequency-Dependent Weighting
; -----------------------------
IF (KEYWORD_SET(wfac)) THEN FOR i=0,n_freq-1 DO rg(i,*) = $
    resid(i,ng) * SQRT( weight(ng) * wfac(i,sky_idx(ng)) / (1.-wf_frac(i,ng)) )
;

; C-Matrix
; --------
covar = rg # TRANSPOSE(rg)
;

; Degrees of Freedom for C-Matrix
; -------------------------------
ndf_covar = cg - n_stripes
;

; Covariance
; ----------
covar = covar / ndf_covar
;

error = 0
;

RETURN,error
END
