Pro FMD_DSP_RSF,error
;

;
;  FMD_DSP_RSF drives the FMD_APPLY procedure to apply destriper
;  corrections to RLFA data, using default coadd weights.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_APPLY
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_IN    =  Directory containing RSF.ISS and RSF_CALSPEC.ISS.
;
;    CSDR$FIRAS_OUT   =  Directory containing RSF_EJG.ISS and RSF_FUNC.ISS,
;                        and where RSF_DSP.ISS will be sent.
;
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_DSP_RSF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 17-Jun-97.
;
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
 PRINT,'FMD_DSP_RSF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_DSP_RSF,error'
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
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_IN     == ' + intrans
PRINT,'CSDR$FIRAS_OUT    == ' + outtrans
PRINT,' '
;

; Restore Stripes and Kernels
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RSF_EJG.ISS'
RESTORE,'csdr$firas_out:rsf_ejg.iss'
;
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RSF_FUNC.ISS'
RESTORE,'csdr$firas_out:rsf_func.iss'
;

; Restore RSF Data
; ----------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RSF.ISS'
RESTORE,'csdr$firas_in:rsf.iss'
;

; Restore IDL Save Set of RSF Undestriped Spectra
; -------------------------------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RSF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:rsf_calspec.iss'
;

; Frequency Correction
; --------------------
fx = f * freq_corr
;

; Frequency Band
; --------------
nf = WHERE((fx GE freq_band(0))and(fx LE freq_band(1)),cf)
IF (cf ne N_ELEMENTS(f_lo)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_DSP_RSF : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Frequency Cut
; -------------
spec = FLOAT(sp(nf,*))
cal_spec = FLOAT(cal_sp(nf,*))
;

; Concatenate SKY/CAL Kernels
; ---------------------------
func = [[sky_func_rsf],[cal_func_rsf]]
;

; Apply Stripes to Compute Destriped Coadd Spectra
; ------------------------------------------------
;
dummy = FMD_APPLY(freq=f_lo,func=func,n_krnl=n_krnl,n_cmn=n_cmn, $
                  spec=TRANSPOSE(spec),cal_spec=TRANSPOSE(cal_spec),$
   		  ejg=ejg_rsf,sky_idx=sky_idx,cal_idx=cal_idx,$
                  dsp=dsp,cal_dsp=cal_dsp,/pure)
;

; Make IDL Save Set
; -----------------
chanscan = 'RSF'
sname = 'csdr$firas_out:rsf_dsp.iss'
SAVE,filename=sname,chanscan,stripes_descrip,f_lo,dsp,cal_dsp
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RSF_DSP.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0 & sky_glitch=0 & cal_glitch=0 & sky_s0=0 & cal_s0=0
glon=0 & glat=0 & sky_lbl=0 & cal_lbl=0 & st_sub=0 & sl_weight_rat=0
chan_label=0 & sky_wgts=0 & sky_wgts_ds=0 & sky_idx=0 & sky_mask=0
dirbe_mask=0 & cal_wgts=0 & cal_wgts_ds=0 & fsl_idx=0
px=0 & tm=0 & sky_dihd=0 & cal_idx=0 & xcal=0 & ical=0 & refh=0
galcut=0 & skyh=0 & dihd=0 & cal_tm=0 & solution=0 & scan=0 & time=0
stripe_order=0 & stripe_id=0 & dirbe_array=0 & lp_order=0 & xtm=0 & cmn_tm=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & cmn_dihd=0 & ref_s0=0 & rms_s0=0
step_up=0 & step_dn=0 & dirbe_cut=0 & max_frac=0 & n_stripes=0 & cmn_bol=0
latcut=0 & loncut=0 & del_temp=0 & freq_corr=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & dihed_cut_min=0 & cvec_mask=0 & cvec_cut=0
d_rsf=0 & n_rsf=0 & b_rsf=0 & l_rsf=0
c_rsf=0 & ndf_ejg_rsf=0 & chi2_rsf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
