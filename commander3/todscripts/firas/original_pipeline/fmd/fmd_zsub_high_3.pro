Pro FMD_ZSUB_HIGH_3,error
;

;
;  FMD_ZSUB_HIGH_3 reads calibrated and zodi Band_3 spectra for the
;  four high chanscans, subtracts the zodi spectra, concatenates the
;  zodi-subtracted spectra, and stores them in IDL save set
;  HIGH_3_CALSPEC_ZSUB.ISS .
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;   Required Logicals :
;
;    CSDR$FIRAS_REF    =  Directory containing FMD_QUALS_HI_3.ISS .
;
;    CSDR$FIRAS_IN     =  Directory containing xxx_3_CALSPEC.ISS and
;                         xxx_3_ZODISPEC.ISS, and where output IDL save
;                         set HIGH_3_CALSPEC_ZSUB.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_ZSUB_HIGH_3,error
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 26-Mar-1997.
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
IF N_Params() ne 1 THEN BEGIN 
 PRINT,'FMD_ZSUB_HIGH_3 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_ZSUB_HIGH_3,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_ZSUB_HIGH_3,error'
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

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hi_3.iss'
;

; Frequency Band
; --------------
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHS_3_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lhs_3_calspec.iss'
PRINT,' '
;
fx = f_hi3
nf = WHERE((fx ge freq_band(0))and(fx le freq_band(1)),cf)
IF (cf le 0) THEN BEGIN
 PRINT,''
 PRINT,'No Data Found in FREQ_BAND Range !'
 PRINT,''
 RETURN
ENDIF
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHS_3_ZODISPEC.ISS'
RESTORE,'csdr$firas_in:lhs_3_zodispec.iss'
PRINT,' '
;
nf = WHERE(fx NE f_hi3,cf)
IF (cf GT 0) THEN BEGIN
 PRINT,''
 PRINT,'LHS_3 : CALSPEC/ZODISPEC Frequency Mismatch !'
 PRINT,''
 RETURN
ENDIF
;
spec0 = TEMPORARY(spec) - TEMPORARY(zodi_spec)
cal_spec0 = TEMPORARY(cal_spec)
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS_3_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rhs_3_calspec.iss'
PRINT,' '
;
fx = f_hi3
nf = WHERE((fx ge freq_band(0))and(fx le freq_band(1)),cf)
IF (cf le 0) THEN BEGIN
 PRINT,''
 PRINT,'No Data Found in FREQ_BAND Range !'
 PRINT,''
 RETURN
ENDIF
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS_3_ZODISPEC.ISS'
RESTORE,'csdr$firas_in:rhs_3_zodispec.iss'
PRINT,' '
;
nf = WHERE(fx NE f_hi3,cf)
IF (cf GT 0) THEN BEGIN
 PRINT,''
 PRINT,'RHS_3 : CALSPEC/ZODISPEC Frequency Mismatch !'
 PRINT,''
 RETURN
ENDIF
;
spec0 = [[spec0],[TEMPORARY(spec) - TEMPORARY(zodi_spec)]]
cal_spec0 = [[cal_spec0],[TEMPORARY(cal_spec)]]
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHF_3_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lhf_3_calspec.iss'
PRINT,' '
;
fx = f_hi3
nf = WHERE((fx ge freq_band(0))and(fx le freq_band(1)),cf)
IF (cf le 0) THEN BEGIN
 PRINT,''
 PRINT,'No Data Found in FREQ_BAND Range !'
 PRINT,''
 RETURN
ENDIF
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHF_3_ZODISPEC.ISS'
RESTORE,'csdr$firas_in:lhf_3_zodispec.iss'
PRINT,' '
;
nf = WHERE(fx NE f_hi3,cf)
IF (cf GT 0) THEN BEGIN
 PRINT,''
 PRINT,'LHF_3 : CALSPEC/ZODISPEC Frequency Mismatch !'
 PRINT,''
 RETURN
ENDIF
;
spec0 = [[spec0],[TEMPORARY(spec) - TEMPORARY(zodi_spec)]]
cal_spec0 = [[cal_spec0],[TEMPORARY(cal_spec)]]
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHF_3_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rhf_3_calspec.iss'
PRINT,' '
;
fx = f_hi3
nf = WHERE((fx ge freq_band(0))and(fx le freq_band(1)),cf)
IF (cf le 0) THEN BEGIN
 PRINT,''
 PRINT,'No Data Found in FREQ_BAND Range !'
 PRINT,''
 RETURN
ENDIF
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHF_3_ZODISPEC.ISS'
RESTORE,'csdr$firas_in:rhf_3_zodispec.iss'
PRINT,' '
;
nf = WHERE(fx NE f_hi3,cf)
IF (cf GT 0) THEN BEGIN
 PRINT,''
 PRINT,'RHF_3 : CALSPEC/ZODISPEC Frequency Mismatch !'
 PRINT,''
 RETURN
ENDIF
;
spec = [[spec0],[TEMPORARY(spec) - TEMPORARY(zodi_spec)]]
cal_spec = [[cal_spec0],[TEMPORARY(cal_spec)]]
;

sname = 'csdr$firas_in:high_3_calspec_zsub.iss'
SAVE,filename=sname,f_hi3,spec,cal_spec
PRINT,' '
PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_3_CALSPEC_ZSUB.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
;

; Return with No Error
; --------------------
error = 0
;

del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0
xtm=0 & cmn_tm=0 & step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0
zodi_array=0 & good_sky_dihd=0 & tmin=0 & tmax=0 & good_cal_dihd=0
hot_cal=0 & cvec_cut=0 & latcut=0 & loncut=0 & n_stripes=0 & max_frac=0
gain_convg=0 & gain_iter=0
badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_lhf=0
badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
;

RETURN
END
