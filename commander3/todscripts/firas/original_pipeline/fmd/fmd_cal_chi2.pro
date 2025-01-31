FUNCTION FMD_CAL_CHI2, xcal=xcal,del_temp=del_temp,freq=freq,cal_idx=cal_idx,$
                       cal_wgts=cal_wgts,cal_spec=cal_spec,wfac=wfac,$
                       cvec=cvec,resid=resid,chi2=chi2
;

;
; FMD_CAL_CHI2F computes calibration coadd residuals and chi-squared, using
;               frequency-dependent coadd weights.
;
;
; KEYWORDS   :
;
;  XCAL (I)         : XCAL Temperature (K)
;
;  DEL_TEMP (I)     : XCAL Offset (K)
;
;  FREQ (I)         : Frequency Array (icm)
;
;  CAL_IDX (I)      : Channel/Scan-Mode Index for each Coadd
;
;  CAL_WGTS (I)     : Coadd Weights
;
;  WFAC (I)         : Frequency-Dependent Weight Factor (Optional)
;
;  CAL_SPEC (I)     : Calibration Coadd Spectra (MJy/sr)
;
;  CVEC (I)         : Combined C-Vector (MJy/sr)
;
;  RESID (O)        : Coadd Residuals (MJy/sr)
;
;  CHI2 (O)         : Coadd Chi-Squared
;
;
;
; RETURNS : ERROR = Error Status
;
;
; REQUIRED LOGICALS : None
;
;
; PROGRAMS Called   : PLANCK
;
;
; HISTORY : Written by Ken Jensen, Hughes STX, 11-Dec-96
;
;
;

; Initialize Error Status
; -----------------------
error = 1
;

; Initialize Arrays
; -----------------
nf = N_ELEMENTS(freq)
nc = N_ELEMENTS(xcal)
resid = FLTARR(nf,nc)  ; Coadd Residual
chi2 = FLTARR(nf,nc)   ; Coadd Chi-Squared
;

; COADD RESIDUALS and CHI-SQUARED
; -------------------------------
;
FOR j=0,nc-1 DO BEGIN

 ; Planck Spectrum
 ; ---------------
 ref = PLANCK(xcal(j)-del_temp,freq,units='i',/mjy)
 ;

 ; Coadd Residuals
 ; ---------------
 resid(*,j) = cal_spec(*,j) - ref
 ;

 ; Coadd Chi-Squared
 ; -----------------
 ;

 ; IF (WFAC and CAL_IDX), use Frequency-Dependent Weights
 ; ------------------------------------------------------
 IF ((KEYWORD_SET(wfac))and(KEYWORD_SET(cal_idx))) THEN BEGIN
  idxx = cal_idx(j)  ; CHANSCAN Index
  FOR k=0,nf-1 DO chi2(k,j) = $
    (resid(k,j)^2.) * cal_wgts(j) * wfac(k,idxx) / cvec(k)^2.
 ENDIF ELSE BEGIN
 ;
 ; IF NOT (WFAC and CAL_IDX), use Frequency-Independent Weights
 ; ------------------------------------------------------------
   FOR k=0,nf-1 DO $
    chi2(k,j) = (resid(k,j)^2.) * cal_wgts(j) / cvec(k)^2.
 ;
 ENDELSE
 ;

ENDFOR
;

; "Good" Cal Coadds
; -----------------
good = WHERE(cal_wgts GT 0.,cgood)
;

; Chi-Squared per DOF Summed over Frequencies (Good Coadds)
; ---------------------------------------------------------
ctg = TOTAL(chi2(*,good),1) / nf
;

; No error
; --------
error = 0
;

RETURN,error
END
