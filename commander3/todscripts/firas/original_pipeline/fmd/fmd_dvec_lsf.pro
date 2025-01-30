Pro FMD_DVEC_LSF,error
;

;
;  FMD_DVEC_LSF drives the FMD_DVECTOR procedure to compute
;  the LLFA D-Vector .
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_DVECTOR
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_LO.ISS and
;                        LSF_WEIGHTS.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing LSF.ISS and LSF_VAR.ISS
;
;    CSDR$FIRAS_OUT   =  Directory containing LSF_CVECTOR.ISS and LSF_EJG.ISS,
;                        where output LSF_DVECTOR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_DVEC_LSF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-97.
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
 PRINT,'FMD_DVEC_LSF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_DVEC_LSF,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_DVEC_LSF,error'
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
RESTORE,'csdr$firas_ref:fmd_quals_lo.iss'
;

; Restore C-Vector
; ----------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LSF_CVECTOR.ISS'
RESTORE,'csdr$firas_out:lsf_cvector.iss'
;

; Restore IDL Save Set of Coadd Variances
; ---------------------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LSF_VAR.ISS'
RESTORE,'csdr$firas_in:lsf_var.iss'
;

; Frequency Cut
; -------------
freq = f_lsf * freq_corr
nf = WHERE((freq GE freq_band(0))and(freq LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(freq_lsf)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_DVEC_LSF : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;
var = var_lsf(nf,*)
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LSF.ISS'
RESTORE,'csdr$firas_in:lsf.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'LSF_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:lsf_weights.iss'
;

; Restore LSF_EJG Save Set
; ------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LSF_EJG.ISS'
RESTORE,'csdr$firas_out:lsf_ejg.iss'
;

n_sky = N_ELEMENTS(px)    ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)  ; Number of CAL Coadds
;

; Coadd Mask
; ----------
mask = cvec_mask
;
badf = WHERE(frac_wgt GE max_frac,cbad)   ; Bad Fraction of Pixel Weight
IF (cbad GT 0) THEN mask(badf) = 0
;
vt = TOTAL(var_lsf,1)
badv = WHERE(vt LE 0.,cbad)               ; Bad Coadd Variance
IF (cbad GT 0) THEN mask(badv) = 0
;
bad_cvec = WHERE(mask NE 1,cbad)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cbad)) + ' Coadds Masked from D-VECTOR'
PRINT,' '
;

; Compute D-Vector
; ----------------
dvec_lsf = $
 FMD_DVECTOR(var=var,mask=mask,nifgs=nifgs,st_sub=st_sub,sky_wgts=sky_wgts_ds,$
             good=good_dvec_lsf,ndf_dvec=ndf_dvec_lsf)
;

; Make LSF_DVECTOR.ISS Save Set
; -----------------------------
chanscan = 'LSF'
freq_dependence = 'N'
dvec_cut = cvec_cut
sname = 'csdr$firas_out:lsf_dvector.iss'
SAVE,filename=sname,chanscan,freq_dependence,dvec_cut,max_frac,freq_lsf,$
                    dvec_lsf,good_dvec_lsf,ndf_dvec_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSF_DVECTOR.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
del_temp=0 & n_stripes=0 & stripes_descrip=0 & var_tm=0 & cal_var_tm=0
tm=0 & cal_tm=0 & ical=0 & refh=0 & skyh=0 & dihd=0 & sky_glitch=0 & sky_s0=0
sky_dihd=0 & cal_glitch=0 & cal_wgts=0 & cal_s0=0 & cal_wgts_ds=0
chan_label=0 & sky_mask=0 & dirbe_mask=0 & time=0 & model_var=0 & lp_orderv=0
cal_nifgs=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0 & latcut=0 & loncut=0
sky_wgts=0 & gain_convg=0 & gain_iter=0 & glat=0 & glon=0 & fsl_idx=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lls=0 & badcoadd_lsf=0
badcal_rls=0 & badcal_rsf=0 & badcal_lls=0 & badcal_lsf=0
sl_weight_rat=0 & var_lbl=0 & cal_var_lbl=0
galcut_llfl=0 & f_llfl=0 & d_llfl=0 & var_llfl=0 & cal_var_llfl=0
galcut_llfs=0 & f_llfs=0 & d_llfs=0 & var_llfs=0 & cal_var_llfs=0
ejg_lsf=0 & c_lsf=0 & ndf_ejg_lsf=0 & chi2_lsf=0
n_lsf=0 & b_lsf=0 & l_lsf=0 & d_lsf=0 & f_lo=0
cvec_lsf=0 & ndf_cvec_lsf=0 & galcut_lsf=0 & cal_var_lsf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
