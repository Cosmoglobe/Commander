Pro FMD_ZODI_SKYMAP,band,error
;

;  FMD_ZODI_SKYMAP reads and concatenates Band_n zodi spectra
;  for the four high chanscans, reads the combined HIGH coadd weights,
;  computes coadd weighted zodi pixel spectra, and stores them in
;  IDL save set HIGH_n_ZODI_SKYMAP.ISS .
;
;
;  ARGUMENTS (I/O)    :
;
;   BAND (I)          :  Frequency Band ( 2, 3, or 4 ).
;
;   ERROR (O)         :  Return Error Status
;
;
;   Required Logicals :
;
;    CSDR$FIRAS_REF    =  Directory containing FMD_QUALS_HI_n.ISS .
;
;    CSDR$FIRAS_IN     =  Directory containing xxx_n_CALSPEC.ISS,
;                         HIGH_WEIGHTS.ISS, and xxx_n_ZODISPEC.ISS,
;                         and where output HIGH_n_SKYMAP.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_ZODI_SKYMAP,2,error
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 19-Jun-1997.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;; BEGIN
;
; Initialize Return Error
; -----------------------
error = 1
;

; Procedure Invoked Correctly ?  If not, signal and RETURN
; --------------------------------------------------------
IF N_Params() ne 2 THEN BEGIN 
 PRINT,'FMD_ZODI_SKYMAP : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_ZODI_SKYMAP,band,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_ZODI_SKYMAP,band,error'
 PRINT,' '
 RETURN
ENDIF
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
intrans = STRUPCASE(intrans)
;
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

sband = STRCOMPRESS(STRING(band),/remove_all)
;

RESTORE,'csdr$firas_in:lhs_' + sband + '_zodispec.iss'
zspec = TEMPORARY(zodi_spec)
;
RESTORE,'csdr$firas_in:rhs_' + sband + '_zodispec.iss'
zspec = [[zspec],[TEMPORARY(zodi_spec)]]
;
RESTORE,'csdr$firas_in:lhf_' + sband + '_zodispec.iss'
zspec = [[zspec],[TEMPORARY(zodi_spec)]]
;
RESTORE,'csdr$firas_in:rhf_' + sband + '_zodispec.iss'
zspec = [[zspec],[TEMPORARY(zodi_spec)]]
;
RESTORE,'csdr$firas_in:high_weights.iss'
;

cf = N_ELEMENTS(zspec) / N_ELEMENTS(px)
;
good = WHERE(sky_wgts_ds GT 0.)
pxg = px(good)
wg = sky_wgts_ds(good)
;
dummy = FMD_PIX_SUM(pixel=pxg,data=zspec(0,good),wgt=wg,cmp_px=cmp_px)
dummy = FLTARR(cf,N_ELEMENTS(cmp_px))
;
FOR i=0,cf-1 DO $
 dummy(i,*) = FMD_PIX_SUM(pixel=pxg,data=zspec(i,good),wgt=wg,cmp_px=cmp_px)
;

IF (band EQ 2) THEN BEGIN
 ;
 PIX2XY,cmp_px,data=dummy,res=6,/six,ras=z_high_2
 ;
 fz_hi2 = f_hi2
 sname = 'csdr$firas_in:high_2_zodi_skymap.iss'
 SAVE,filename=sname,fz_hi2,z_high_2
 ;
ENDIF
;

IF (band EQ 3) THEN BEGIN
 ;
 PIX2XY,cmp_px,data=dummy,res=6,/six,ras=z_high_3
 ;
 fz_hi3 = f_hi3
 sname = 'csdr$firas_in:high_3_zodi_skymap.iss'
 SAVE,filename=sname,fz_hi3,z_high_3
 ;
ENDIF
;

IF (band EQ 4) THEN BEGIN
 ;
 PIX2XY,cmp_px,data=dummy,res=6,/six,ras=z_high_4
 ;
 fz_hi4 = f_hi4
 sname = 'csdr$firas_in:high_4_zodi_skymap.iss'
 SAVE,filename=sname,fz_hi4,z_high_4
 ;
ENDIF
;

PRINT,' '
PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_'+sband+'_ZODI_SKYMAP.ISS" Created.'
PRINT,' '
;

chanscan=0 & freq_high=0 & pixel_wgt=0 & frac_wgt=0 & cal_wgts_ds=0
sky_idx=0 & cal_idx=0 & scale_factor=0
;

error=0
;

RETURN
END
