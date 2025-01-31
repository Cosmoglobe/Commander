Pro FMD_COVAR_HI_23,error
;

;
;  FMD_COVAR_HI_23 drives the FMD_COVAR procedure to create an IDL
;  save set of combined high-channel Band_2/Band_3 covariance.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS
;                        and FMD_QUALS_HI_3.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing HIGH_WEIGHTS.ISS.
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_2_SKY_RESID.ISS,
;                        and HIGH_3_SKY_RESID.ISS, and where output
;                        HIGH_23_COVAR.ISS will be sent.
;
;
;    PROGRAMS Called  :  FMD_COVAR_OFF
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_COVAR_HI_23,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 12-May-97.
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
 PRINT,'FMD_COVAR_HI_23 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_COVAR_HI_23,error'
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

; Restore FMD Qualifiers
; ----------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'FMD_QUALS_HI_2.ISS'
RESTORE,'csdr$firas_ref:fmd_quals_hi_2.iss'
nst = [n_stripes(1)]
;
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'FMD_QUALS_HI_3.ISS'
RESTORE,'csdr$firas_ref:fmd_quals_hi_3.iss'
nst = [nst,n_stripes(1)]
;

; Restore Coadd Residuals
; -----------------------
PRINT,'Restoring ' + outtrans(0) + 'HIGH_2_SKY_RESID.ISS'
RESTORE,'csdr$firas_out:high_2_sky_resid.iss'
resid1 = TEMPORARY(sky_resid_high2)
;
PRINT,'Restoring ' + outtrans(0) + 'HIGH_3_SKY_RESID.ISS'
RESTORE,'csdr$firas_out:high_3_sky_resid.iss'
resid2 = TEMPORARY(sky_resid_high3)
;

nfreq = [N_ELEMENTS(freq_high_2),N_ELEMENTS(freq_high_3)]
;

; Restore Coadd Weights
; ---------------------
PRINT,'Restoring HIGH_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:high_weights.iss'
;

; Compute C-Matrices
; ------------------
;
PRINT,'Calling FMD_COVAR_OFF'

error =  FMD_COVAR_OFF(nfreq=nfreq,resid1=resid1,resid2=resid2,$
                       weight=sky_wgts_ds,mask=sky_mask,cvec_cut=cvec_cut,$
                       px=px,n_stripes=nst,max_frac=max_frac,w_frac=frac_wgt,$
                       covar=covar_high_23,ndf_covar=ndf_covar_hi_23)
;

IF (error NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HI_23 : Error Returned from FMD_COVAR_OFF !'
 RETURN
ENDIF
;

; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HIGH'
freq_dependence = 'N'
n_stripes = nst
sname = 'csdr$firas_out:high_23_covar.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,n_stripes,max_frac,$
                    freq_high_2,freq_high_3,covar_high_23,ndf_covar_hi_23
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_23_COVAR.ISS" Created.'
PRINT,' '
;

del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0
cmn_tm=0 & step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0
zodi_array=0 & good_sky_dihd=0 & tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0
latcut=0 & loncut=0 & gain_convg=0 & gain_iter=0 & chan_label=0 & pixel_wgt=0
stripes_descrip=0 & cal_wgts_ds=0 & cal_idx=0 & freq_band=0 & sky_idx=0
badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_lhf=0
scale_factor=0 & badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
freq_high=0 & sp_high_2=0 & z_high_2=0 & sp_high_3=0 & z_high_3=0
;

RETURN
END
