Pro FMD_WGTS_HRES,error
;

;  FMD_WGTS_HRES drives the FMD_PIXEL_WGT procedure to create
;  IDL save sets of combined low channel weights.
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
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HR.ISS,
;                        CVECTORS_HR.ISS, and xxx_WEIGHTS.ISS.
;
;    CSDR$FIRAS_IN    =  Directory containing HRES.ISS, and where output
;                        HRES_WEIGHTS.ISS and HRES_WEIGHTS_F.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_WGTS_HRES,error
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
 PRINT,'FMD_WGTS_HRES : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_WGTS_HRES,error'
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
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
intrans = STRUPCASE(intrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,' '
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hr.iss'
;

PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HRES.ISS'
RESTORE,'csdr$firas_in:hres.iss'
px0 = TEMPORARY(px)
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'LLF_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:llf_weights.iss'
sky_wgts0 = sky_wgts_ds
cal_wgts0 = cal_wgts_ds
;

PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'RLF_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:rlf_weights.iss'
sky_wgts1 = sky_wgts_ds
cal_wgts1 = cal_wgts_ds
;

; Check CHANSCAN Index
; --------------------
dummy = [0*sky_wgts0,0*sky_wgts1+1]
nx = WHERE(dummy NE sky_idx,cx)
IF (cx GT 0) THEN BEGIN
 PRINT,'FMD_WGTS_HRES : Error in SKY_IDX !'
 READ,ncont
 RETURN
ENDIF
;
dummy = [0*cal_wgts0,0*cal_wgts1+1]
nx = WHERE(dummy NE cal_idx,cx)
IF (cx GT 0) THEN BEGIN
 PRINT,'FMD_WGTS_HRES : Error in CAL_IDX !'
 READ,ncont
 RETURN
ENDIF
;

; Restore C-Vectors for Individual Low CHANSCANs
; ----------------------------------------------
RESTORE,'csdr$firas_ref:cvectors_hr.iss'
;

; Frequency Band
; --------------
nf = WHERE((f_hr ge freq_band(0))and(f_hr le freq_band(1)),cf)
freq_hres = f_hr(nf)
;

; Frequency-Independent Coadd Weight Adjustment
; ---------------------------------------------
scale_factor = [1.,TOTAL(1/cvec_rlf(nf)^2) / TOTAL(1/cvec_llf(nf)^2)]
;

sky_wgts_ds =  [sky_wgts0,sky_wgts1] * scale_factor(sky_idx)

cal_wgts_ds =  [cal_wgts0,cal_wgts1] * scale_factor(cal_idx)
;

; Number of Pixels with Good Data
; -------------------------------
good = WHERE(sky_wgts_ds GT 0.)
hpix = HISTOGRAM(px0(good),min=0)
npix = WHERE(hpix GT 0,cpix)
;

; Normalization Defined so that total sky weight = Number of Pixels
; ------------------------------------------------------------------
norm = FLOAT(cpix) / TOTAL(sky_wgts_ds)
;

; Normalized Weights
; ------------------
sky_wgts_ds = norm * sky_wgts_ds
cal_wgts_ds = norm * cal_wgts_ds
;

cvecs = [[cvec_llf(nf)],[cvec_rlf(nf)]]
;

;  Pixel Weights and Coadd Fractional Weights
;  ------------------------------------------
px = TEMPORARY(px0)
error = FMD_PIXEL_WGT(pixel=px,weights=sky_wgts_ds,sky_idx=sky_idx, $
                      cvecs=cvecs,frac_wgt=frac_wgt,frac_fwgt=frac_fwgt, $
                      px_wgt=pixel_wgt)
;

IF (error NE 0) THEN BEGIN
 PRINT,''
 PRINT,'FMD_PIXEL_WGT Returned With Error !'
 PRINT,''
 RETURN
ENDIF
;

; Make HRES_WEIGHTS.ISS Save Set
; --------------------------------
chanscan = 'HRES'
sname = 'csdr$firas_in:hres_weights.iss'
SAVE,filename=sname,chanscan,freq_hres,px,sky_wgts_ds,pixel_wgt,frac_wgt,$
                    cal_wgts_ds,sky_idx,cal_idx,scale_factor
PRINT,' '
PRINT,'IDL Save Set "' + intrans(0) + 'HRES_WEIGHTS.ISS" Created.'
PRINT,' '
;

; Make HRES_WEIGHTS_F.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_in:hres_weights_f.iss'
SAVE,filename=sname,chanscan,px,frac_fwgt,scale_factor
PRINT,' '
PRINT,'IDL Save Set "'+ intrans(0) + 'HRES_WEIGHTS_F.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
hot_cal=0 & gain_convg=0 & gain_iter=0
del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & maxfrac=0
cmn_tm=0 & step_up=0 & step_dn=0 & chan_label=0
ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0 & ref_dihd=0 & n_stripes=0
dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & cvec_cut=0 & latcut=0 & loncut=0
lp_order=0 & xtm=0 & zodi_array=0 & good_sky_dihd=0 & tmin=0 & tmax=0
good_cal_dihd=0 & max_frac=0 & tm=0 & glat=0 & glon=0 & sky_dihd=0
sky_glitch=0 & sky_s0=0 & scan=0 & sky_wgts=0 & nifgs=0 & st_sub=0
xcal=0 & cal_tm=0 & ical=0 & refh=0 & skyh=0 & dihd=0 & cal_s0=0
st_sub=0 & cal_glitch=0 & cal_wgts=0 & cal_nifgs=0
sky_wgts=0 & cal_wgts=0 & model_var=0 & tm_up=0 & tm_dn=0 & gpow=0 & bpow=0
nscan=0 & delt=0 & tmcut=0 & lp_orderv=0 & fsl_idx=0 & sky_lbl=0 & cal_lbl=0
badcoadd_rlf=0 & badcoadd_llf=0 & badcal_rlf=0 & badcal_llf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
