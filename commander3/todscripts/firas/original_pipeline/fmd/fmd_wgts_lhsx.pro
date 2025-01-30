Pro FMD_WGTS_LHSX,error
;
;
;  FMD_WGTS_LHSX drives the FMD_PIXEL_WGT procedures to create IDL save sets
;  of LHSS weights, using the default bad coadd indices.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_PIXEL_WGT
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS,
;                        FMD_QUALS_HI_4.ISS, and FMD_QUALS_DEFAULT.ISS,
;                        and where LHSX_WEIGHTS.ISS will be sent.
;
;    CSDR$FIRAS_IN    =  Directory containing LHS.ISS .
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_WGTS_LHSX,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-1997.
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
 PRINT,'FMD_WGTS_LHSX : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_WGTS_LHSX,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 RETURN
ENDIF
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
intrans = STRUPCASE(intrans)
;
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
outtrans = STRUPCASE(outtrans)
;
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hi_2.iss'
freq1 = freq_band(0)
;
RESTORE,'csdr$firas_ref:fmd_quals_hi_4.iss'
freq2 = freq_band(1)
;
freq_band = [freq1,freq2]
;

; Restore Default Bad Coadd Indices
; ---------------------------------
RESTORE,'csdr$firas_ref:fmd_quals_default.iss'
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHS.ISS'
RESTORE,'csdr$firas_in:lhs.iss'
;

; List of FMD Qualifiers
; ----------------------
PRINT,' '
PRINT,'FMD Qualifiers are :'
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
PRINT,'TMIN =',tmin                     ; Minimum controllable temps
PRINT,'TMAX =',tmax                     ; Maximum controllable temps
PRINT,'GOOD_CAL_DIHD =',good_cal_dihd   ; Good Cal DIHED temps
PRINT,'GOOD_SKY_DIHD =',good_sky_dihd   ; Good Sky DIHED temps
PRINT,' '
PRINT,'BADCOADD_LHS =',badcoadd_lhs     ; Bad LHS SKY Coadds
PRINT,' '
PRINT,'BADCAL_LHS =',badcal_lhs         ; Bad LHS CAL Coadds
PRINT,' '
PRINT,' '
;

n_sky = N_ELEMENTS(px)     ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)   ; Number of CAL Coadds
;

; Build BADCOADD Index
; --------------------
badcoadd = BYTARR(n_sky)
IF (badcoadd_lhs(0) ge 0) THEN badcoadd(badcoadd_lhs) = 1B
;

; Build BADCAL Index
; --------------------
badcal = BYTARR(n_cal)
IF (badcal_lhs(0) ge 0) THEN badcal(badcal_lhs) = 1B
;

; Initialize Destriper Coadd Weights
; ----------------------------------
sky_wgts_ds = sky_wgts
cal_wgts_ds = cal_wgts
;

; De-weight Bad Sky Coadds
; ------------------------
bad = WHERE(badcoadd eq 1,cbad)
IF (cbad GT 0) THEN sky_wgts_ds(bad) = 0.
;

; Deweight Sky Spectra with Temps out-of-Range
; --------------------------------------------
bad = WHERE( sky_dihd LT good_sky_dihd(0) OR sky_dihd GT good_sky_dihd(1) )
IF (bad(0) ge 0) THEN sky_wgts_ds(bad) = 0.
;

; De-weight Bad Cal Coadds
; ------------------------
bad = WHERE(badcal eq 1,cbad)
IF (cbad GT 0) THEN cal_wgts_ds(bad) = 0.
;

; De-weight Calibration Spectra Early in Mission
; ----------------------------------------------
bad = WHERE(cal_tm lt 28.35,cbad)
IF (cbad GT 0) THEN cal_wgts_ds(bad) = 0.
;

; Deweight Calibration Spectra with Temps out-of-Range
; ----------------------------------------------------
bad = WHERE(xcal LE tmin(0) OR xcal GE tmax(0) OR $
           ical LE tmin(1)  OR ical GE tmax(1) OR $
           skyh LE tmin(2)  OR skyh GE tmax(2) OR $
           refh LE tmin(3)  OR refh GE tmax(3) OR $
           dihd LT good_cal_dihd(0) OR dihd GT good_cal_dihd(1) )
IF (bad(0) ge 0) THEN cal_wgts_ds(bad) = 0.
;


;  Pixel Weights and Coadd Fractional Weights
;  ------------------------------------------
error = FMD_PIXEL_WGT(pixel=px,weights=sky_wgts_ds, $
                      frac_wgt=frac_wgt,px_wgt=pixel_wgt)
;

IF (error NE 0) THEN BEGIN
 PRINT,''
 PRINT,'FMD_PIXEL_WGT Returned With Error !'
 PRINT,''
 RETURN
ENDIF
;

; Make LHSX_WEIGHTS.ISS Save Set
; ---------------------------------
chanscan = 'LHS'
sname = 'csdr$firas_ref:lhsx_weights.iss'
SAVE,filename=sname,chanscan,px,sky_wgts,sky_wgts_ds,pixel_wgt,$
                    frac_wgt,cal_wgts,cal_wgts_ds
PRINT,' '
PRINT,'IDL Save Set "' + reftrans(0) + 'LHSX_WEIGHTS.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
hot_cal=0 & gain_convg=0 & gain_iter=0 & sky_glitch=0 & tm=0 & glat=0 & glon=0
sky_s0=0 & sky_bol_volt=0 & cal_glitch=0 & cal_bol_volt=0 & chanscan_id=0
cal_s0=0 & del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & maxfrac=0
xtm=0 & cmn_tm=0 & step_up=0 & step_dn=0 & chan_label=0 & f_hi3=0 & st_sub=0
ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0 & ref_dihd=0 & n_stripes=0
dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & cvec_cut=0 & latcut=0 & loncut=0
zodi_array=0 & max_frac=0 & solution=0 & f=0 & galcut=0 & scan=0 & cal_lbl=0
sky_lbl=0 & time=0 & lp_order=0 & nifgs=0 & cal_nifgs=0 & sl_weight_rat=0
fsl_idx=0 & n_lhs=0 & b_lhs=0 & l_lhs=0 & d_lhs=0
badcoadd_lhf=0 & badcoadd_rhs=0 & badcoadd_rhf=0
badcal_lhf=0 & badcal_rhs=0 & badcal_rhf=0
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_lsf=0 & badcoadd_rsf=0
badcal_lls=0 & badcal_rls=0 & badcal_lsf=0 & badcal_rsf=0
badcoadd_llf=0 & badcoadd_rlf=0 & badcal_llf=0 & badcal_rlf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
