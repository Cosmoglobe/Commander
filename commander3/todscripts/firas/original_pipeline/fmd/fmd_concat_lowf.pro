Pro FMD_CONCAT_LOWF,file_select,error
;

;
;  FMD_CONCAT_LOWF reads data for the four low chanscans, applies a
;  frequency cut, concatenates the arrays, and stores them in an IDL
;  save set for later processing by FMD_DESTRIPER.
;
;
;  ARGUMENTS (I/O)    :
;
;   FILE_SELECT (I)   :  2-element array controls which data sets
;                        will be concatenated
;                        ( [xxx,xxx_CALSPEC] ).
;
;   ERROR (O)         :  Return Error Status
;
;
;   Required Logicals :
;
;    CSDR$FIRAS_REF    =  Directory containing FMD_QUALS_LO.ISS.
;
;    CSDR$FIRAS_IN     =  Directory containing xxx.ISS and xxx_CALSPEC.ISS,
;                         and where output LOWF.ISS and LOWF_CALSPEC.ISS
;                         will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CONCAT_LOWF,[1,1],error
;
;                    (will read the xxx and xxx_CALSPEC data for each of
;                     LLS, RLS, LSF, and RSF, obtain a frequency cut from
;                     FMD_QUALS_LO.ISS, apply the frequency cut to the
;                     real part of the calibrated spectra, concatenate the
;                     real part of the arrays, and save the concatenated 
;                     arrays in LOWF.ISS and LOWF_CALSPEC.ISS.)
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 27-May-1997.
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
 PRINT,'FMD_CONCAT_LOWF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CONCAT_LOWF,file_select,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CONCAT_LOWF,[1,1],error'
 PRINT,' '
 RETURN
ENDIF
;

; Is FILE_SELECT Valid ?
; ----------------------
IF (TOTAL(file_select(0:1)) LE 0.) THEN BEGIN
 PRINT,'FMD_CONCAT_LOWF : FILE_SELECT Error !'
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

; IF (FILE_SELECT(0)) THEN Make LOWF.ISS
; --------------------------------------
IF (file_select(0) GT 0) THEN BEGIN
;
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'LLS.ISS'
 RESTORE,'csdr$firas_in:lls.iss'
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
 PRINT,'Restoring ' + intrans(0) + 'RLS.ISS'
 RESTORE,'csdr$firas_in:rls.iss'
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
 PRINT,'Restoring ' + intrans(0) + 'LSF.ISS'
 RESTORE,'csdr$firas_in:lsf.iss'
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
 PRINT,'Restoring ' + intrans(0) + 'RSF.ISS'
 RESTORE,'csdr$firas_in:rsf.iss'
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
 chan_label = 'LLS_RLS_LSF_RSF'
 sname = 'csdr$firas_in:lowf.iss'
 SAVE,filename=sname,chan_label,sky_idx,px,tm,glat,glon,sky_dihd,$
      sky_glitch,sky_s0,scan,sky_wgts,nifgs,st_sub,cal_idx,xcal,cal_tm,$
      ical,refh,skyh,dihd,cal_glitch,cal_s0,cal_wgts,cal_nifgs,$
      fsl_idx,sky_lbl,cal_lbl
 PRINT,' '
 PRINT,'IDL Save Set "' + intrans(0) + 'LOWF.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (file_select(1) EQ 0) THEN RETURN
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_lo.iss'
;
; Frequency Band
; --------------
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
;

PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LLS_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lls_calspec.iss'
PRINT,' '
;
f_lo = f * freq_corr  ; Frequency Correction
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
IF (cf LE 0) THEN BEGIN
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
PRINT,'Restoring ' + intrans(0) + 'RLS_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rls_calspec.iss'
PRINT,' '
;
f_lo = f * freq_corr  ; Frequency Correction
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf2)
IF (cf2 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'LLS/RLS Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
spec0 = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LSF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lsf_calspec.iss'
PRINT,' '
;
f_lo = f * freq_corr  ; Frequency Correction
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf3)
IF (cf3 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'LLS/LSF Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
cal_spec0 = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
spec0 = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
;
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RSF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rsf_calspec.iss'
PRINT,' '
;
f_lo = f * freq_corr  ; Frequency Correction
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf4)
IF (cf4 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'LLS/RSF Frequency Arrays are Not Compatible !'
 PRINT,''
 RETURN
ENDIF
;
cal_spec = [[TEMPORARY(cal_spec0)],[FLOAT(cal_sp(nf,*))]]
spec = [[TEMPORARY(spec0)],[FLOAT(sp(nf,*))]]
f_lo = f_lo(nf)
;

sname = 'csdr$firas_in:lowf_calspec.iss'
SAVE,filename=sname,f_lo,spec,cal_spec
PRINT,' '
PRINT,'IDL Save Set "' + intrans(0) + 'LOWF_CALSPEC.ISS" Created.'
PRINT,' '
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
badcal_rls=0 & badcal_rsf=0 & badcal_lls=0 & badcal_lsf=0
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lsf=0
d_lls=0 & n_lls=0 & b_lls=0 & l_lls=0
d_rls=0 & n_rls=0 & b_rls=0 & l_rls=0
d_lsf=0 & n_lsf=0 & b_lsf=0 & l_lsf=0
d_rsf=0 & n_rsf=0 & b_rsf=0 & l_rsf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
