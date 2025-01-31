Pro FMD_CONCAT_HIGH,file_select,band,error
;

;
;  FMD_CONCAT_HIGH reads data for the four high chanscans, applies a
;  frequency cut, concatenates the arrays, and stores them in an IDL
;  save set for later processing by FMD_DESTRIPER.
;
;
;  ARGUMENTS (I/O)    :
;
;   FILE_SELECT (I)   :  3-element array controls which data sets
;                        will be concatenated
;                        ( [xxx,xxx_CALSPEC,xxx_DIFFSPEC] ).
;
;   BAND (I)          :  Frequency Band (1, 2, 3, or 4)
;
;   ERROR (O)         :  Return Error Status
;
;
;   Required Logicals :
;
;    CSDR$FIRAS_REF    =  Directory containing FMD_QUALS_HI_x.ISS.
;                         FMD_DIRBE_FUNC_xxx_n.ISS .
;
;    CSDR$FIRAS_IN     =  Directory containing xxx.ISS and xxx_CALSPEC.ISS,
;                         and where output HIGH.ISS, HIGH_n_CALSPEC.ISS,
;                         and HIGH_n_DIFFSPEC.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CONCAT_HIGH,[1,1,0],2,error
;
;                    (will read the xxx and xxx_CALSPEC data for each of
;                     LHS, LHF, RHS, and RHF, obtain a frequency cut from
;                     FMD_QUALS_HI_2.ISS, apply the frequency cut to the
;                     real part of the calibrated spectra, concatenate the
;                     real part of the arrays, and save the concatenated 
;                     arrays in HIGH.ISS and HIGH_2_CALSPEC.ISS.)
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 13-Jun-1997.
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
IF N_Params() ne 3 THEN BEGIN 
 PRINT,'FMD_CONCAT_HIGH : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CONCAT_HIGH,file_select,band,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CONCAT_HIGH,[1,1,0],2,error'
 PRINT,' '
 RETURN
ENDIF
;

; Is FILE_SELECT Valid ?
; ----------------------
IF (TOTAL(file_select(0:2)) LE 0.) THEN BEGIN
 PRINT,'FMD_CONCAT_HIGH : FILE_SELECT Error !'
 RETURN
ENDIF
;

; Is BAND Valid ?
; ---------------
IF ((band LT 1)or(band GT 4)) THEN BEGIN
 PRINT,'FMD_CONCAT_HIGH : BAND Must Be 1, 2, 3, or 4 !'
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

; IF (FILE_SELECT(0)) THEN Make HIGH.ISS
; --------------------------------------
IF (file_select(0) GT 0) THEN BEGIN
;
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'LHS.ISS'
 RESTORE,'csdr$firas_in:lhs.iss'
 PRINT,' '
;
 sky_idx = BYTE(0*px)
 px0 = TEMPORARY(px)
 tm0 = TEMPORARY(tm)
 glon0 = TEMPORARY(glon)
 glat0 = TEMPORARY(glat)
 sky_dihd0 = TEMPORARY(sky_dihd)
 sky_glitch0 = TEMPORARY(sky_glitch)
 sky_s00 = TEMPORARY(sky_s0)
 scan0 = TEMPORARY(scan)
 sky_wgts0 = TEMPORARY(sky_wgts)
 nifgs0 = TEMPORARY(nifgs)
 sub0 = TEMPORARY(st_sub)
 idx0 = TEMPORARY(fsl_idx)
 lbl0 = TEMPORARY(sky_lbl)
 ;
 cal_idx = BYTE(0*xcal)
 xcal0 = TEMPORARY(xcal)
 cal_tm0 = TEMPORARY(cal_tm)
 ical0 = TEMPORARY(ical)
 refh0 = TEMPORARY(refh)
 skyh0 = TEMPORARY(skyh)
 dihd0 = TEMPORARY(dihd)
 cal_glitch0 = TEMPORARY(cal_glitch)
 cal_s00 = TEMPORARY(cal_s0)
 cal_wgts0 = TEMPORARY(cal_wgts)
 cal_nifgs0 = TEMPORARY(cal_nifgs)
 clbl0 = TEMPORARY(cal_lbl)
 ;

 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'RHS.ISS'
 RESTORE,'csdr$firas_in:rhs.iss'
 PRINT,' '
 ;
 sky_idx = [[sky_idx],BYTE(0*px)+1B]
 px0 = [TEMPORARY(px0),px]
 tm0 = [TEMPORARY(tm0),tm]
 glon0 = [TEMPORARY(glon0),TEMPORARY(glon)]
 glat0 = [TEMPORARY(glat0),TEMPORARY(glat)]
 sky_dihd0 = [TEMPORARY(sky_dihd0),TEMPORARY(sky_dihd)]
 sky_glitch0 = [TEMPORARY(sky_glitch0),TEMPORARY(sky_glitch)]
 sky_s00 = [TEMPORARY(sky_s00),TEMPORARY(sky_s0)]
 scan0 = [TEMPORARY(scan0),TEMPORARY(scan)]
 sky_wgts0 = [TEMPORARY(sky_wgts0),TEMPORARY(sky_wgts)]
 nifgs0 = [TEMPORARY(nifgs0),TEMPORARY(nifgs)]
 sub0 = [sub0,TEMPORARY(st_sub)]
 idx0 = [idx0,TEMPORARY(fsl_idx)]
 lbl0 = [lbl0,TEMPORARY(sky_lbl)]
 ;
 cal_idx = [cal_idx,BYTE(0*xcal)+1B]
 xcal0 = [TEMPORARY(xcal0),xcal]
 cal_tm0 = [TEMPORARY(cal_tm0),cal_tm]
 ical0 = [TEMPORARY(ical0),ical]
 refh0 = [TEMPORARY(refh0),refh]
 skyh0 = [TEMPORARY(skyh0),skyh]
 dihd0 = [TEMPORARY(dihd0),dihd]
 cal_glitch0 = [TEMPORARY(cal_glitch0),TEMPORARY(cal_glitch)]
 cal_s00 = [TEMPORARY(cal_s00),TEMPORARY(cal_s0)]
 cal_wgts0 = [TEMPORARY(cal_wgts0),TEMPORARY(cal_wgts)]
 cal_nifgs0 = [TEMPORARY(cal_nifgs0),TEMPORARY(cal_nifgs)]
 clbl0 = [clbl0,TEMPORARY(cal_lbl)]
 ;

 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'LHF.ISS'
 RESTORE,'csdr$firas_in:lhf.iss'
 PRINT,' '
 ;
 sky_idx = [[sky_idx],BYTE(0*px)+2B]
 px0 = [TEMPORARY(px0),px]
 tm0 = [TEMPORARY(tm0),tm]
 glon0 = [TEMPORARY(glon0),TEMPORARY(glon)]
 glat0 = [TEMPORARY(glat0),TEMPORARY(glat)]
 sky_dihd0 = [TEMPORARY(sky_dihd0),TEMPORARY(sky_dihd)]
 sky_glitch0 = [TEMPORARY(sky_glitch0),TEMPORARY(sky_glitch)]
 sky_s00 = [TEMPORARY(sky_s00),TEMPORARY(sky_s0)]
 scan0 = [TEMPORARY(scan0),TEMPORARY(scan)]
 sky_wgts0 = [TEMPORARY(sky_wgts0),TEMPORARY(sky_wgts)]
 nifgs0 = [TEMPORARY(nifgs0),TEMPORARY(nifgs)]
 sub0 = [sub0,TEMPORARY(st_sub)]
 idx0 = [idx0,TEMPORARY(fsl_idx)]
 lbl0 = [lbl0,TEMPORARY(sky_lbl)]
 ;
 cal_idx = [cal_idx,BYTE(0*xcal)+2B]
 xcal0 = [TEMPORARY(xcal0),xcal]
 cal_tm0 = [TEMPORARY(cal_tm0),cal_tm]
 ical0 = [TEMPORARY(ical0),ical]
 refh0 = [TEMPORARY(refh0),refh]
 skyh0 = [TEMPORARY(skyh0),skyh]
 dihd0 = [TEMPORARY(dihd0),dihd]
 cal_glitch0 = [TEMPORARY(cal_glitch0),TEMPORARY(cal_glitch)]
 cal_s00 = [TEMPORARY(cal_s00),TEMPORARY(cal_s0)]
 cal_wgts0 = [TEMPORARY(cal_wgts0),TEMPORARY(cal_wgts)]
 cal_nifgs0 = [TEMPORARY(cal_nifgs0),TEMPORARY(cal_nifgs)]
 clbl0 = [clbl0,TEMPORARY(cal_lbl)]
 ;

 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'RHF.ISS'
 RESTORE,'csdr$firas_in:rhf.iss'
 PRINT,' '
 ;
 sky_idx = [[sky_idx],BYTE(0*px)+3B]
 px = [TEMPORARY(px0),px]
 tm = [TEMPORARY(tm0),tm]
 glon = [TEMPORARY(glon0),TEMPORARY(glon)]
 glat = [TEMPORARY(glat0),TEMPORARY(glat)]
 sky_dihd = [TEMPORARY(sky_dihd0),TEMPORARY(sky_dihd)]
 sky_glitch = [TEMPORARY(sky_glitch0),TEMPORARY(sky_glitch)]
 sky_s0 = [TEMPORARY(sky_s00),TEMPORARY(sky_s0)]
 scan = [TEMPORARY(scan0),TEMPORARY(scan)]
 sky_wgts = [TEMPORARY(sky_wgts0),TEMPORARY(sky_wgts)]
 nifgs = [TEMPORARY(nifgs0),TEMPORARY(nifgs)]
 st_sub = [sub0,TEMPORARY(st_sub)]
 fsl_idx = [idx0,TEMPORARY(fsl_idx)]
 sky_lbl = [lbl0,TEMPORARY(sky_lbl)]
 ;
 cal_idx = [cal_idx,BYTE(0*xcal)+3B]
 xcal = [TEMPORARY(xcal0),xcal]
 cal_tm = [TEMPORARY(cal_tm0),cal_tm]
 ical = [TEMPORARY(ical0),ical]
 refh = [TEMPORARY(refh0),refh]
 skyh = [TEMPORARY(skyh0),skyh]
 dihd = [TEMPORARY(dihd0),dihd]
 cal_glitch = [TEMPORARY(cal_glitch0),TEMPORARY(cal_glitch)]
 cal_s0 = [TEMPORARY(cal_s00),TEMPORARY(cal_s0)]
 cal_wgts = [TEMPORARY(cal_wgts0),TEMPORARY(cal_wgts)]
 cal_nifgs = [TEMPORARY(cal_nifgs0),TEMPORARY(cal_nifgs)]
 cal_lbl = [clbl0,TEMPORARY(cal_lbl)]
 ;
 chan_label = 'LHS_RHS_LHF_RHF'
 sname = 'csdr$firas_in:high.iss'
 SAVE,filename=sname,chan_label,sky_idx,px,tm,glat,glon,sky_dihd,$
      sky_glitch,sky_s0,scan,sky_wgts,nifgs,cal_idx,xcal,cal_tm,$
      ical,refh,skyh,dihd,cal_glitch,cal_s0,cal_wgts,cal_nifgs,$
      st_sub,fsl_idx,sky_lbl,cal_lbl
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'HIGH.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (TOTAL(file_select(1:2)) GT 0.) THEN BEGIN
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
ENDIF
;

; IF (FILE_SELECT(1)) THEN Make HIGH_CALSPEC_x.ISS
; ------------------------------------------------
IF (file_select(1) GT 0) THEN BEGIN
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
 cal_spec0 = FLOAT(cal_sp(nf,*))
 spec0 = FLOAT(sp(nf,*))
;
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'RHS_CALSPEC.ISS'
 RESTORE,'csdr$firas_in:rhs_calspec.iss'
 PRINT,' '
;
 f_hi = f * freq_corr  ; Frequency Correction
 nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf2)
 IF (cf2 ne cf) THEN BEGIN
  PRINT,''
  PRINT,'LHS/RHS Frequency Arrays are Not Compatible !'
  PRINT,''
  RETURN
 ENDIF
;
 cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
 spec0 = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
;
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'LHF_CALSPEC.ISS'
 RESTORE,'csdr$firas_in:lhf_calspec.iss'
 PRINT,' '
;
 f_hi = f * freq_corr  ; Frequency Correction
 nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf3)
 IF (cf3 ne cf) THEN BEGIN
  PRINT,''
  PRINT,'LHS/LHF Frequency Arrays are Not Compatible !'
  PRINT,''
  RETURN
 ENDIF
;
 cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
 spec0 = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
;
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'RHF_CALSPEC.ISS'
 RESTORE,'csdr$firas_in:rhf_calspec.iss'
 PRINT,' '
;
 f_hi = f * freq_corr  ; Frequency Correction
 nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf4)
 IF (cf4 ne cf) THEN BEGIN
  PRINT,''
  PRINT,'LHS/RHF Frequency Arrays are Not Compatible !'
  PRINT,''
  RETURN
 ENDIF
;
 cal_spec = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
 spec = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
;

 IF (band eq 1) THEN BEGIN
  f_hi1 = f_hi(nf)
  sname = 'csdr$firas_in:high_1_calspec.iss'
  SAVE,filename=sname,f_hi1,spec,cal_spec
  PRINT,' '
  PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_1_CALSPEC.ISS" Created.'
  PRINT,' '
 ENDIF
;
 IF (band eq 2) THEN BEGIN
  f_hi2 = f_hi(nf)
  sname = 'csdr$firas_in:high_2_calspec.iss'
  SAVE,filename=sname,f_hi2,spec,cal_spec
  PRINT,' '
  PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_2_CALSPEC.ISS" Created.'
  PRINT,' '
 ENDIF
;
 IF (band eq 3) THEN BEGIN
  f_hi3 = f_hi(nf)
  sname = 'csdr$firas_in:high_3_calspec.iss'
  SAVE,filename=sname,f_hi3,spec,cal_spec
  PRINT,' '
  PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_3_CALSPEC.ISS" Created.'
  PRINT,' '
 ENDIF
;
 IF (band eq 4) THEN BEGIN
  f_hi4 = f_hi(nf)
  sname = 'csdr$firas_in:high_4_calspec.iss'
  SAVE,filename=sname,f_hi4,spec,cal_spec
  PRINT,' '
  PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_4_CALSPEC.ISS" Created.'
  PRINT,' '
 ENDIF
;
ENDIF
;

; IF (FILE_SELECT(2)) THEN Make HIGH_n_DIFFSPEC.ISS
; -------------------------------------------------
;
IF (file_select(2) LE 0) THEN BEGIN
 error = 0
 RETURN
ENDIF
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHS_DIFFSPEC.ISS'
RESTORE,'csdr$firas_in:lhs_diffspec.iss'
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
cal_spec0 = FLOAT(diff_cal_sp(nf,*))
spec0 = FLOAT(diff_sp(nf,*))
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS_DIFFSPEC.ISS'
RESTORE,'csdr$firas_in:rhs_diffspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf2)
IF (cf2 ne cf) THEN BEGIN
 PRINT,''
 PRINT,'LHS/RHS Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(diff_cal_sp(nf,*))]]
spec0 = [[TEMPORARY(spec0)],[FLOAT(diff_sp(nf,*))]]
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHF_DIFFSPEC.ISS'
RESTORE,'csdr$firas_in:lhf_diffspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf3)
IF (cf3 ne cf) THEN BEGIN
 PRINT,''
 PRINT,'LHS/LHF Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(diff_cal_sp(nf,*))]]
spec0 = [[TEMPORARY(spec0)],[FLOAT(diff_sp(nf,*))]]
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHF_DIFFSPEC.ISS'
RESTORE,'csdr$firas_in:rhf_diffspec.iss'
PRINT,' '
;
f_hi = f * freq_corr  ; Frequency Correction
nf = WHERE((f_hi ge freq_band(0))and(f_hi le freq_band(1)),cf4)
IF (cf4 ne cf) THEN BEGIN
 PRINT,''
 PRINT,'LHS/RHF Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
diff_cal_spec = [[TEMPORARY(cal_spec0)],[FLOAT(diff_cal_sp(nf,*))]]
diff_spec = [[TEMPORARY(spec0)],[FLOAT(diff_sp(nf,*))]]
;

IF (band eq 1) THEN BEGIN
 f_hi1 = f_hi(nf)
 sname = 'csdr$firas_in:high_1_diffspec.iss'
 SAVE,filename=sname,f_hi1,diff_spec,diff_cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_1_DIFFSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 2) THEN BEGIN
 f_hi2 = f_hi(nf)
 sname = 'csdr$firas_in:high_2_diffspec.iss'
 SAVE,filename=sname,f_hi2,diff_spec,diff_cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_2_DIFFSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 3) THEN BEGIN
 f_hi3 = f_hi(nf)
 sname = 'csdr$firas_in:high_3_diffspec.iss'
 SAVE,filename=sname,f_hi3,diff_spec,diff_cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_3_DIFFSPEC.ISS" Created.'
 PRINT,' '
ENDIF
;
IF (band eq 4) THEN BEGIN
 f_hi4 = f_hi(nf)
 sname = 'csdr$firas_in:high_4_diffspec.iss'
 SAVE,filename=sname,f_hi4,diff_spec,diff_cal_spec
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'HIGH_4_DIFFSPEC.ISS" Created.'
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
zodi_array=0 & n_stripes=0 & max_frac=0 & time=0
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
