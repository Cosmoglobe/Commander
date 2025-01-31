Pro FMD_QUALS,bselect
;

;
;  FMD_QUALS drives the FMD_SAVE_QUALS procedure, which makes an IDL
;  save set of qualifiers for running FMD (FIRAS Mega-Destriper).
;
;
;  ARGUMENTS :
;
;   BSELECT (I) : = 0 for Low Channel LRES
;                   1 for High Channel Band_1
;                   2 for High Channel Band_2
;                   3 for High Channel Band_3
;                   4 for High Channel Band_4
;                   5 for Low Channel HRES
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_DEFAULT.ISS, and
;                        where FMD_QUALS_xxx.ISS will be sent.
;
;
;  Programs Called    :  FMD_SAVE_QUALS
;
;
;  EXAMPLE : $ UIDL
;            UIDL> FMD_QUALS,2
;
;                  (will call FMD_SAVE_QUALS to store FMD qualifiers for
;                   high channel Band_2 destriping in IDL save set
;                   FMD_QUALS_HI_2.ISS .)
;
;
;  Written by :  Ken Jensen,  Hughes STX,  23-Apr-97
;
;
;  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;  $$                                               $$
;  $$  This file is designed to be edited whenever  $$
;  $$  a modification of the qualifiers is desired. $$
;  $$                                               $$
;  $$  All changes to FMD qualifiers should be      $$
;  $$  made inside this program !                   $$
;  $$                                               $$
;  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;
;
;  Modification History : K.Jensen, HSTX, 04-May-97,
;                         Re-define BADCOADD, BADCAL parameters
;
;  


if N_Params() ne 1 then begin
 print,'FMD_QUALS,bselect'
 return
endif
;

; Qualifiers Common to All Bands
; ------------------------------
x_tm = [1,-1]              ; Range of X(time) for Legendre Polynomials Pn(X)
refd = 2.3                 ; Reference Dihedral Temperature
db_cut = [1.e6,1.e6,1.e6]  ; Maximum Allowed DIRBE Gradients [G10,G9,G8]
sky_dihd = [1.5,5.5]       ; Allowed SKY Dihedral Temps
tminx = [2.6,2.6,2.6,2.6]  ; Minimum Allowed CAL Temps [XCAL,ICAL,SKYH,REFH]
tmaxx = [2.8,2.8,2.8,2.8]  ; Maximum Allowed CAL Temps [XCAL,ICAL,SKYH,REFH]
cal_dihd = [1.5,5.5]       ; Allowed CAL Dihedral Temps
hot = [4.,5.,1.5,1.1]      ; Parameters for De-Weighting Hot CAL Coadds
convg = 0.01               ; Gain Convergence Threshold
iter = 20                  ; Maximum Allowed Gain Iterations
;

; LO Qualifers
; ------------
IF (bselect EQ 0) THEN BEGIN
;
 fb = [2.,21.5]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,0,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = INDGEN(7) + 1        ; Orders of Legendre polynomial time kernals.
 horn = [1,1]              ; Horn Stripe Flags [6K,4K]
 pow = [2.]                ; Powers of Dihedral Temperature kernals.
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-1.974,-1.6261,-1.974,-1.6261] ; Reference Bolometer Responsivity
 rmsb = 0.0462                          ; RMS Bolometer Responsivity
;
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [14,29]        ; Number of Stripes ([Individual,Combined])
;
 lat = 5  &  lon = 30      ; LAT/LON Mask for Destriper
 cv_cut = [15,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; HI_1 Qualifers
; --------------
IF (bselect EQ 1) THEN BEGIN
;
 fb = [2.,21.5]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,0,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = INDGEN(7) + 1        ; Orders of Legendre polynomial time kernals.
 horn = [1,1]              ; Horn Stripe Flags [6K,4K]
 pow = [2.]                ; Powers of Dihedral Temperature kernals.
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-.5235,-.5607,-.4919,-.5312]  ; Reference Bolometer Responsivity
 rmsb = 0.014                          ; RMS Bolometer Responsivity
;
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [14,29]        ; Number of Stripes ([Individual,Combined])
;
 lat = 5  &  lon = 30      ; LAT/LON Mask for Destriper
 cv_cut = [15,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; HI_2 Qualifers
; --------------
IF (bselect EQ 2) THEN BEGIN
;
 fb = [20.,45.]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,0,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = INDGEN(3) + 1        ; Orders of Legendre polynomial time kernals.
 horn = [1,0]              ; Horn Stripe Flags [6K,4K]
 pow = [2.]                ; Powers of Dihedral Temperature kernals.
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-.5235,-.5607,-.4919,-.5312]  ; Reference Bolometer Responsivity
;
;rmsb = 0.014                          ; RMS Bolometer Responsivity
 rmsb = 0.                 ; RMS Bolometer Responsivity (Disabled if EQ 0.)
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [8,17]         ; Number of Stripes ([Individual,Combined])
;
 lat = 8  &  lon = 100     ; LAT/LON Mask for Destriper
 cv_cut = [20,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; HI_3 Qualifers
; --------------
IF (bselect EQ 3) THEN BEGIN
;
 fb = [45.,70.]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,1,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = [1]                  ; Orders of Legendre polynomial time kernals.
 horn = [1,0]              ; Horn Stripe Flags [6K,4K]
 pow = [2.]                ; Powers of Dihedral Temperature kernals
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-.5235,-.5607,-.4919,-.5312]  ; Reference Bolometer Responsivity
;
 rmsb = 0.014                          ; RMS Bolometer Responsivity
; rmsb = 0.                 ; RMS Bolometer Responsivity (Disabled if EQ 0.)
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [7,19]         ; Number of Stripes ([Individual,Combined])
;
 lat = 8  &  lon = 100     ; LAT/LON Mask for Destriper
 cv_cut = [20,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; HI_4 Qualifers
; --------------
IF (bselect EQ 4) THEN BEGIN
;
 fb = [70.,98.]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,1,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = [-1]                 ; Orders of Legendre polynomial time kernals.
 horn = [0,0]              ; Horn Stripe Flags [6K,4K]
 pow = [0.]                ; Powers of Dihedral Temperature kernals.
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-.5235,-.5607,-.4919,-.5312]  ; Reference Bolometer Responsivity
;
; rmsb = 0.014                          ; RMS Bolometer Responsivity
 rmsb = 0.                 ; RMS Bolometer Responsivity (Disabled if EQ 0.)
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [3,6]          ; Number of Stripes ([Individual,Combined])
;
 lat = 8  &  lon = 100     ; LAT/LON Mask for Destriper
 cv_cut = [20,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; HRES Qualifers
; --------------
IF (bselect EQ 5) THEN BEGIN
;
 fb = [0.5,21.5]            ; Frequency Range (icm)
 dirb = [1,1,0]            ; DIRBE Stripes [G10,G9,G8]
 zod = [0,0,0]             ; ZODI Stripes [Z10,Z9,Z8]
 lp = INDGEN(7) + 1        ; Orders of Legendre polynomial time kernals.
 horn = [1,0]              ; Horn Stripe Flags [6K,4K]
 pow = [2.]                ; Powers of Dihedral Temperature kernals.
 cut = 0.                  ; Dihedral Temperature Cut
;
 refb = [-1.974,-1.6261]   ; Reference Bolometer Responsivity
 rmsb = 0.0462             ; RMS Bolometer Responsivity
;
 cmn = ['Y','N','N']       ; Common Stripe Flags (TM,DIHD,BOL)
 nstripes = [9,17]         ; Number of Stripes ([Individual,Combined])
;
 lat = 5  &  lon = 30      ; LAT/LON Mask for Destriper
 cv_cut = [15,180]         ; LAT/LON Mask for C-Vector
 maxfrac = 0.9             ; Maximum Coadd Fraction of Pixel-Weight for COVAR
;
ENDIF
;

; Make the Qualifiers Save Set
; ----------------------------
err = FMD_SAVE_QUALS(bselect,fb=fb,dirb=dirb,db_cut=db_cut,zod=zod,$
      lp=lp,x_tm=x_tm,horn=horn,pow=pow,cut=cut,refd=refd,refb=refb,rmsb=rmsb,$
      sky_dihd=sky_dihd,tminx=tminx,tmaxx=tmaxx,cal_dihd=cal_dihd,$
      hot=hot,cv_cut=cv_cut,lat=lat,lon=lon,convg=convg,iter=iter,cmn=cmn,$
      nstripes=nstripes,maxfrac=maxfrac)
;

; Re-Define Parameters
; --------------------
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_lsf=0 & badcoadd_rsf=0
badcoadd_lhs=0 & badcoadd_rhs=0 & badcoadd_lhf=0 & badcoadd_rhf=0
badcoadd_llf=0 & badcoadd_rlf=0
badcal_lls=0 & badcal_rls=0 & badcal_lsf=0 & badcal_rsf=0
badcal_lhs=0 & badcal_rhs=0 & badcal_lhf=0 & badcal_rhf=0
badcal_llf=0 & badcal_rlf=0
;

RETURN
END
