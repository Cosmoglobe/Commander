FUNCTION FMD_DVECTOR, var=var,mask=mask,nifgs=nifgs,st_sub=st_sub,$
                      sky_wgts=sky_wgts,good=good,ndf_dvec=ndf_dvec
;
;
; FMD_DVECTOR computes the D-Vector .
;
;
; KEYWORDS   :
;
;  VAR (I)         : Coadd Variance ((MJy/sr)^2.)
;
;  MASK (I)        : Coadd Mask (if = 0, coadd is excluded)
;
;  NIFGS (I)       : Coadd Size
;
;  ST_SUB (I)      : Secondary_Template Subtraction (0=NO,1=YES)
;
;  SKY_WGTS (I)    : Coadd Weights
;
;  GOOD (O)        : Index of Coadds used in D-Vector calculation
;
;  NDF_CVEC (O)    : Degrees of Freedom in D-Vector calculation
;
;
; RETURNS :  DVEC  =  D-Vector (MJy/sr)
;
;
; REQUIRED LOGICALS : None
;
;
; PROGRAMS Called   : None
;
;
; HISTORY : Written by Ken Jensen, Hughes STX, 23-Apr-97
;
;
;

; Initialize Arrays
; -----------------
;
nc = N_ELEMENTS(nifgs)
nf = N_ELEMENTS(var) / nc
;
dvec = FLTARR(nf)           ; D-Vector
;

; Good Coadds for D-Vector
; ------------------------
good = WHERE((sky_wgts GT 0.)and(mask GT 0),cgood)
;

; IF No Coadds, return with error
; -------------------------------
IF (cgood LE 0) THEN BEGIN
 PRINT,'FMD_DVECTOR : No Coadds for D-Vector Calculation !'
 RETURN,dvec
ENDIF
;
sg = STRCOMPRESS(STRING(cgood))
PRINT,' '
PRINT,'D-Vector will be Computed from'+sg+' Coadds'
PRINT,' '
;

; Compute D-Vector
; ----------------
num = var(*,good) # ((nifgs(good)-1-st_sub(good))*sky_wgts(good))
den = TOTAL(nifgs(good)-1-st_sub(good))
dvec = FLOAT(SQRT(num/den))
;

; Degrees of Freedom in D-Vector Calculation
; ------------------------------------------
ndf_dvec = LONG(den)


RETURN,dvec
END
