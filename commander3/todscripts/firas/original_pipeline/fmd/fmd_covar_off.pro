FUNCTION FMD_COVAR_OFF,nfreq=nfreq,resid1=resid1,resid2=resid2,weight=weight,$
                     mask=mask,px=px,cvec_cut=cvec_cut,n_stripes=n_stripes,$
                     max_frac=max_frac,wfac=wfac,w_frac=w_frac,wf_frac=wf_frac,$
                     sky_idx=sky_idx,covar=covar,ndf_covar=ndf_covar
;

;
;  FMD_COVAR_OFF computes the covariance of destriped coadd spectra.
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

; Sanity Checks
; -------------
IF (N_ELEMENTS(nfreq) NE N_ELEMENTS(n_stripes)) THEN BEGIN
 PRINT,'FMD_COVAR : Error in N_STRIPES !'
 help,nfreq,n_stripes
 READ,ncont
 RETURN,error
ENDIF
;
n_freq1 = N_ELEMENTS(resid1) / N_ELEMENTS(weight)
IF (n_freq1 NE nfreq(0)) THEN BEGIN
 PRINT,'FMD_COVAR_OFF : Error in Array Dimensions !'
 RETURN,error
ENDIF
;
n_freq2 = N_ELEMENTS(resid2) / N_ELEMENTS(weight)
IF (n_freq2 NE nfreq(1)) THEN BEGIN
 PRINT,'FMD_COVAR_OFF : Error in Array Dimensions !'
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
rt = TOTAL(resid1,1) + TOTAL(resid2,1)  ; Total Residual per Coadd
;
ng = WHERE((rt NE 0.)and(w_frac GT 0.)and(w_frac LT max_frac)and(cmask EQ 1),cg)
IF (cg LE 0) THEN BEGIN
 PRINT,'FMD_COVAR_OFF : No Good Coadds'
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


; Residuals for Good Coadds
; -------------------------
rg1 = resid1(*,ng)
rg2 = resid2(*,ng)
;

; Frequency-Independent Weighting
; -------------------------------
IF (NOT KEYWORD_SET(wfac)) THEN FOR i=0,n_freq1-1 DO $
    rg1(i,*) = rg1(i,*) * SQRT(weight(ng)/(1. - w_frac(ng)))
;
IF (NOT KEYWORD_SET(wfac)) THEN FOR i=0,n_freq2-1 DO $
    rg2(i,*) = rg2(i,*) * SQRT(weight(ng)/(1. - w_frac(ng)))
;

; Frequency-Dependent Weighting
; -----------------------------
IF (KEYWORD_SET(wfac)) THEN FOR i=0,n_freq1-1 DO rg1(i,*) = $
    rg1(i,*) * SQRT( weight(ng) * wfac(i,sky_idx(ng)) / (1.-wf_frac(i,ng)) )
;
IF (KEYWORD_SET(wfac)) THEN FOR i=0,n_freq2-1 DO rg2(i,*) = $
    rg2(i,*) * SQRT( weight(ng) * wfac(i,sky_idx(ng)) / (1.-wf_frac(i,ng)) )
;

; C-Matrix
; --------
covar = rg1 # TRANSPOSE(rg2)
;

; Degrees of Freedom for C-Matrix
; -------------------------------
ndf_covar = cg - 0.5*TOTAL(n_stripes)
;

; Covariance
; ----------
covar = covar / ndf_covar
;

error = 0
;

RETURN,error
END
