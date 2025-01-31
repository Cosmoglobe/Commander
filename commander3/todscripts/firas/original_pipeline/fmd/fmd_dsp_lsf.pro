Pro FMD_DSP_LSF,error
;

;
;  FMD_DSP_LSF drives the FMD_APPLY procedure to apply destriper
;  corrections to LLFA data, using default coadd weights.
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
;    CSDR$FIRAS_IN    =  Directory containing LSF.ISS and LSF_CALSPEC.ISS.
;
;    CSDR$FIRAS_OUT   =  Directory containing LSF_EJG.ISS and LSF_FUNC.ISS,
;                        and where LSF_DSP.ISS will be sent.
;
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_DSP_LSF,error
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
 PRINT,'FMD_DSP_LSF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_DSP_LSF,error'
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
PRINT,'Restoring ' + outtrans(0) + 'LSF_EJG.ISS'
RESTORE,'csdr$firas_out:lsf_ejg.iss'
;
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LSF_FUNC.ISS'
RESTORE,'csdr$firas_out:lsf_func.iss'
;

; Restore LSF Data
; ----------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LSF.ISS'
RESTORE,'csdr$firas_in:lsf.iss'
;

; Restore IDL Save Set of LSF Undestriped Spectra
; -------------------------------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LSF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lsf_calspec.iss'
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
 PRINT,'FMD_DSP_LSF : Error in Frequency Array !'
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
func = [[sky_func_lsf],[cal_func_lsf]]
;

; Apply Stripes to Compute Destriped Coadd Spectra
; ------------------------------------------------
;
dummy = FMD_APPLY(freq=f_lo,func=func,n_krnl=n_krnl,n_cmn=n_cmn, $
                  spec=TRANSPOSE(spec),cal_spec=TRANSPOSE(cal_spec),$
   		  ejg=ejg_lsf,sky_idx=sky_idx,cal_idx=cal_idx,$
                  dsp=dsp,cal_dsp=cal_dsp,/pure)
;

; Make IDL Save Set
; -----------------
chanscan = 'LSF'
sname = 'csdr$firas_out:lsf_dsp.iss'
SAVE,filename=sname,chanscan,stripes_descrip,f_lo,dsp,cal_dsp
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSF_DSP.ISS" Created.'
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
d_lsf=0 & n_lsf=0 & b_lsf=0 & l_lsf=0
c_lsf=0 & ndf_ejg_lsf=0 & chi2_lsf=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
