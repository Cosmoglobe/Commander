Pro FMD_DSP_LOWF,error
;

;
;  FMD_DSP_LOWF drives the FMD_APPLY procedure to apply destriper
;  corrections to combined low channel data.
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
;    CSDR$FIRAS_IN   =  Directory containing LOWF.ISS and LOWF_CALSPEC.ISS .
;
;    CSDR$FIRAS_OUT  =  Directory containing LOWF_EJG.ISS and LOWF_FUNC.ISS,
;                       and where LOWF_DSP.ISS will be sent.
;
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_DSP_LOWF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-1997.
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
 PRINT,'FMD_DSP_LOWF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_DSP_LOWF,error'
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
PRINT,'Restoring ' + outtrans(0) + 'LOWF_EJG.ISS'
RESTORE,'csdr$firas_out:lowf_ejg.iss'
;
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LOWF_FUNC.ISS'
RESTORE,'csdr$firas_out:lowf_func.iss'
;

; Restore LOWF Data
; -----------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LOWF.ISS'
RESTORE,'csdr$firas_in:lowf.iss'
;


; Restore Undestriped Spectra
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LOWF_CALSPEC.ISS'
RESTORE,'csdr$firas_in:lowf_calspec.iss'
;

; Frequency Band
; --------------
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
IF (cf ne N_ELEMENTS(f_lo)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_DSP_LOWF : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Concatenate SKY/CAL Kernels
; ---------------------------
func = [[sky_func_lowf],[cal_func_lowf]]
;

; Apply Stripes to Compute Destriped Coadd Spectra
; ------------------------------------------------
;
dummy = FMD_APPLY(freq=f_lo,func=func,n_krnl=n_krnl,n_cmn=n_cmn, $
                  spec=TRANSPOSE(spec),cal_spec=TRANSPOSE(cal_spec),$
   		  ejg=ejg_lowf,sky_idx=sky_idx,cal_idx=cal_idx,$
                  dsp=dsp,cal_dsp=cal_dsp,/pure)
;

; Make IDL Save Set
; -----------------
chanscan = 'LOWF'
sname = 'csdr$firas_out:lowf_dsp.iss'
SAVE,filename=sname,chanscan,stripes_descrip,f_lo,dsp,cal_dsp
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LOWF_DSP.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0 & sky_glitch=0 & cal_glitch=0 & sky_s0=0 & cal_s0=0
glon=0 & glat=0 & sky_lbl=0 & cal_lbl=0 & sky_bol_volt=0 & cal_bol_volt=0
chan_label=0 & badcoadd=0 & sky_wgts=0 & sky_wgts_ds=0 & sky_idx=0 & sky_mask=0
dirbe_mask=0 & badcal=0 & cal_wgts=0 & cal_wgts_ds=0 & pcvr=0 & rect=0 & diag=0
px=0 & cmp_px=0 & tm=0 & sky_dihd=0 & cal_idx=0 & xcal=0 & ical=0 & refh=0
sl_weight_rat=0 & galcut=0 & chanscan_id=0 & skyh=0 & dihd=0 & st_sub=0
gain_convg=0 & gain_iter=0 & chi2_loch=0 & cal_tm=0 & fsl_idx=0
scan=0 & sky_wgts_fsl=0 & cal_wgts_fsl=0 & cvec_mask=0 & cvec_cut=0
c_lowf=0 & ndf_ejf_lowf=0 & chi2_lowf=0
n_lowf=0 & l_lowf=0 & b_lowf=0 & gamma_lowf=0
stripe_order=0 & stripe_id=0 & dirbe_array=0 & lp_order=0 & xtm=0 & cmn_tm=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & cmn_dihd=0 & ref_s0=0 & rms_s0=0
step_up=0 & step_dn=0 & dirbe_cut=0 & max_frac=0 & n_stripes=0 & cmn_bol=0
latcut=0 & loncut=0 & del_temp=0 & freq_corr=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & ndf_ejg_lowf=0 & dihed_cut_min=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
