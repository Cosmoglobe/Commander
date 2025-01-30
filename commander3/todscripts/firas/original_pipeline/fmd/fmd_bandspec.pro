Pro FMD_BANDSPEC,band,error
;

;
;  FMD_BANDSPEC reads spectra for the four high chanscans, applies a
;  frequency cut, and stores them in IDL save sets.
;
;
;  ARGUMENTS (I/O)    :
;
;   BAND (I)          :  Frequency Band (1, 2, 3, or 4)
;
;   ERROR (O)         :  Return Error Status
;
;
;   Required Logicals :
;
;    CSDR$FIRAS_REF    =  Directory containing FMD_QUALS_HI_n.ISS .
;
;    CSDR$FIRAS_IN     =  Directory containing xxx.ISS and xxx_CALSPEC.ISS,
;                         and where output xxx_n_CALSPEC.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_BANDSPEC,2,error
;
;                    (will read the xxx and xxx_CALSPEC data for each of
;                     LHS, LHF, RHS, and RHF, obtain a frequency cut from
;                     FMD_QUALS_HI_2.ISS, apply the frequency cut to the
;                     real part of the calibrated spectra, and save each
;                     array of spectra in the xxx_2_CALSPEC.ISS.)
;
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
IF N_Params() ne 2 THEN BEGIN 
 PRINT,'FMD_BANDSPEC : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_BANDSPEC,band,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_BANDSPEC,2,error'
 PRINT,' '
 RETURN
ENDIF
;

; Is BAND Valid ?
; ---------------
IF ((band LT 1)or(band GT 4)) THEN BEGIN
 PRINT,'FMD_BANDSPEC : BAND Must Be 1, 2, 3, or 4 !'
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
sband = STRTRIM(STRCOMPRESS(STRING(band)),2)
RESTORE,'csdr$firas_ref:fmd_quals_hi_'+sband+'.iss'
;

; Frequency Band
; --------------
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
;

; LHS
; ---
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHS_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lhs_calspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf)
 IF (cf le 0) THEN BEGIN
  PRINT,''
  PRINT,'No Data Found in FREQ_BAND Range !'
  PRINT,''
  RETURN
 ENDIF
;
cal_spec = FLOAT(cal_sp(nf,*))
spec = FLOAT(sp(nf,*))
;

IF (band eq 1) THEN BEGIN
 f_hi1 = f_hi(nf)
 sname = 'csdr$firas_in:lhs_1_calspec.iss'
 SAVE,filename=sname,f_hi1,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_1_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 2) THEN BEGIN
 f_hi2 = f_hi(nf)
 sname = 'csdr$firas_in:lhs_2_calspec.iss'
 SAVE,filename=sname,f_hi2,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_2_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 3) THEN BEGIN
 f_hi3 = f_hi(nf)
 sname = 'csdr$firas_in:lhs_3_calspec.iss'
 SAVE,filename=sname,f_hi3,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_3_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 4) THEN BEGIN
 f_hi4 = f_hi(nf)
 sname = 'csdr$firas_in:lhs_4_calspec.iss'
 SAVE,filename=sname,f_hi4,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHS_4_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;

; RHS
; ---
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rhs_calspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf)
 IF (cf le 0) THEN BEGIN
  PRINT,''
  PRINT,'No Data Found in FREQ_BAND Range !'
  PRINT,''
  RETURN
 ENDIF
;
cal_spec = FLOAT(cal_sp(nf,*))
spec = FLOAT(sp(nf,*))
;

IF (band eq 1) THEN BEGIN
 f_hi1 = f_hi(nf)
 sname = 'csdr$firas_in:rhs_1_calspec.iss'
 SAVE,filename=sname,f_hi1,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_1_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 2) THEN BEGIN
 f_hi2 = f_hi(nf)
 sname = 'csdr$firas_in:rhs_2_calspec.iss'
 SAVE,filename=sname,f_hi2,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_2_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 3) THEN BEGIN
 f_hi3 = f_hi(nf)
 sname = 'csdr$firas_in:rhs_3_calspec.iss'
 SAVE,filename=sname,f_hi3,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_3_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 4) THEN BEGIN
 f_hi4 = f_hi(nf)
 sname = 'csdr$firas_in:rhs_4_calspec.iss'
 SAVE,filename=sname,f_hi4,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHS_4_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;

; LHF
; ---
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lhf_calspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf)
 IF (cf le 0) THEN BEGIN
  PRINT,''
  PRINT,'No Data Found in FREQ_BAND Range !'
  PRINT,''
  RETURN
 ENDIF
;
cal_spec = FLOAT(cal_sp(nf,*))
spec = FLOAT(sp(nf,*))
;

IF (band eq 1) THEN BEGIN
 f_hi1 = f_hi(nf)
 sname = 'csdr$firas_in:lhf_1_calspec.iss'
 SAVE,filename=sname,f_hi1,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_1_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 2) THEN BEGIN
 f_hi2 = f_hi(nf)
 sname = 'csdr$firas_in:lhf_2_calspec.iss'
 SAVE,filename=sname,f_hi2,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_2_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 3) THEN BEGIN
 f_hi3 = f_hi(nf)
 sname = 'csdr$firas_in:lhf_3_calspec.iss'
 SAVE,filename=sname,f_hi3,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_3_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 4) THEN BEGIN
 f_hi4 = f_hi(nf)
 sname = 'csdr$firas_in:lhf_4_calspec.iss'
 SAVE,filename=sname,f_hi4,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LHF_4_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;

; RHF
; ---
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rhf_calspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf)
 IF (cf le 0) THEN BEGIN
  PRINT,''
  PRINT,'No Data Found in FREQ_BAND Range !'
  PRINT,''
  RETURN
 ENDIF
;
cal_spec = FLOAT(cal_sp(nf,*))
spec = FLOAT(sp(nf,*))
;

IF (band eq 1) THEN BEGIN
 f_hi1 = f_hi(nf)
 sname = 'csdr$firas_in:rhf_1_calspec.iss'
 SAVE,filename=sname,f_hi1,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_1_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 2) THEN BEGIN
 f_hi2 = f_hi(nf)
 sname = 'csdr$firas_in:rhf_2_calspec.iss'
 SAVE,filename=sname,f_hi2,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_2_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 3) THEN BEGIN
 f_hi3 = f_hi(nf)
 sname = 'csdr$firas_in:rhf_3_calspec.iss'
 SAVE,filename=sname,f_hi3,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_3_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 4) THEN BEGIN
 f_hi4 = f_hi(nf)
 sname = 'csdr$firas_in:rhf_4_calspec.iss'
 SAVE,filename=sname,f_hi4,spec,cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'RHF_4_CALSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0 & sky_glitch=0 & cal_glitch=0 & sky_s0=0 & cal_s0=0
glon=0 & glat=0 & sky_lbl=0 & cal_lbl=0 & sky_bol_volt=0 & cal_bol_volt=0
solution=0 & good_cal_dihd=0 & good_sky_dihd=0 & step_up=0 & step_dn=0
del_temp=0 & tau=0 & dihed_cut_min=0 & dihed_pow=0 & tmin=0 & tmax=0
dirbe_array=0 & dirbe_cut=0 & pow_tm=0 & tau0=0 & ref_tm=0 & cmn_tm=0
ref_dihd=0 & dihed_cut=0 & hot_cal=0 & cvec_cut=0 & latcut=0 & loncut=0
gain_convg=0 & gain_iter=0 & diff_px=0 & sl_weight_rat=0 & galcut=0
lp_order=0 & xtm=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & cmn_dihd=0
zodi_array=0 & n_stripes=0 & max_frac=0 & px=0 & tm=0 & cal_tm=0
badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
badcoadd_lhs=0 & badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhf=0
d_lhs=0 & n_lhs=0 & b_lhs=0 & l_lhs=0
d_rhs=0 & n_rhs=0 & b_rhs=0 & l_rhs=0
d_lhf=0 & n_lhf=0 & b_lhf=0 & l_lhf=0
d_rhf=0 & n_rhf=0 & b_rhf=0 & l_rhf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
