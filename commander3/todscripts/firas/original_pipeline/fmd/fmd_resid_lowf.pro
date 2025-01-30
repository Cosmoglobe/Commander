Pro FMD_RESID_LOWF,error
;
;
;  FMD_RESID_LOWF drives the FMD_RESID procedure to compute pixel spectra
;  and coadd residuals for combined low channel destriped data,
;  using frequency-independent coadd weights.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_RESID
;
;
 ;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_LO.ISS and
;                        CVECTORS_LO.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing LOWF.ISS, LOWF_WEIGHTS.ISS,
;                        and LOWF_WEIGHTS_F.ISS .
;
;    CSDR$FIRAS_OUT   =  Directory containing LOWF_EJG.ISS and
;                        LOWF_DSP.ISS, and where IDL save set
;                        LOWF_SKY_RESID.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_RESID_LOWF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-1997
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
 PRINT,'FMD_RESID_LOWF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_RESID_LOWF,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_RESID_LOWF,error'
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

; Frequency-Dependent Weight Flag
; -------------------------------
freq_dependence = 'N'
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_lo.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LOWF_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:lowf_weights.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'Restoring ' + intrans(0) + 'LOWF_WEIGHTS_F.ISS'
 RESTORE,'csdr$firas_in:lowf_weights_f.iss'
ENDIF
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LOWF_DSP.ISS'
RESTORE,'csdr$firas_out:lowf_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_lo)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_RESID_LOWF : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;
nf2 = WHERE((freq_lowf GE freq_band(0))and(freq_lowf LE freq_band(1)),cf2)
IF (cf2 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'FMD_RESID_LOWF : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Frequency Cut
; -------------
freq = f_lo(nf)
;

IF (freq_dependence EQ 'Y') THEN BEGIN

 ; FRAC_FWGT Frequency Cut
 ; -----------------------
 frac_fwgt = frac_fwgt(nf2,*)
 ;

 ; Restore Individual C-Vectors
 ; ----------------------------
 RESTORE,'csdr$firas_ref:cvectors_lo.iss'
 ;

 ; Are Frequencies Compatible ?
 ; ----------------------------
 nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
 IF (cf NE N_ELEMENTS(freq)) THEN BEGIN
  PRINT,'FMD_RESID_LOWF : C-Vector Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 aa = MAX(ABS(freq-f_lo(nf)))
 IF (aa NE 0.) THEN BEGIN
  PRINT,'FMD_RESID_LOWF : C-Vector Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;

 ; Apply Frequency Cut to C-Vectors
 ; --------------------------------
 cvec_lls = cvec_lls(nf)
 cvec_rls = cvec_rls(nf)
 cvec_lsf = cvec_lsf(nf)
 cvec_rsf = cvec_rsf(nf)
 ;

 ; Concatenate C-Vectors
 ; ---------------------
 cvecs = [[cvec_lls],[cvec_rls],[cvec_lsf],[cvec_rsf]]
 ;

 ; Frequency Dependent Weight Adjustment
 ; -------------------------------------
 ;
 ; Normalization defined so that LLSS weights are unaffected
 ; ---------------------------------------------------------
 wfac = 0.*cvecs + 1.
 fac0 = cvecs(*,0)^2. * scale_factor(0)
 wfac(*,1) = fac0 / cvecs(*,1)^2. / scale_factor(1)
 wfac(*,2) = fac0 / cvecs(*,2)^2. / scale_factor(2)
 wfac(*,3) = fac0 / cvecs(*,3)^2. / scale_factor(3)
 ;

ENDIF
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LOWF.ISS'
RESTORE,'csdr$firas_in:lowf.iss'
;

; Restore LOWF_EJG Save Set
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LOWF_EJG.ISS'
RESTORE,'csdr$firas_out:lowf_ejg.iss'
;

n_sky = N_ELEMENTS(px)    ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)  ; Number of CAL Coadds
;

nmask = WHERE(dirbe_mask EQ 0,cmask)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cmask)) + ' Coadds Masked by DIRBE_CUT'
PRINT,' '
;

; Compute Pixel Spectra and Coadd Residuals
; -----------------------------------------
;
IF (freq_dependence EQ 'Y') THEN error = $
    FMD_RESID(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,sky_idx=sky_idx,$
              mask=dirbe_mask,w_frac=frac_wgt,wf_frac=frac_fwgt,$
              px_spec=sp_lowf,resid=sky_resid_lowf)
;

IF (freq_dependence NE 'Y') THEN error = $
    FMD_RESID(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,sky_idx=sky_idx,$
              mask=dirbe_mask,w_frac=frac_wgt,$
              px_spec=sp_lowf,resid=sky_resid_lowf)
;

; Make LOWF_SKY_RESID.ISS Save Set
; ----------------------------------
chanscan = 'LOWF'
chan_label = 'LLS_RLS_LSF_RSF'
freq_lowf = f_lo
;
sname = 'csdr$firas_out:lowf_sky_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,$
                    sky_idx,sky_wgts_ds,frac_wgt,sky_mask,freq_lowf,$
                    sky_resid_lowf,sp_lowf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LOWF_SKY_RESID.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
fsl_idx=0 & sky_lbl=0 & cal_lbl=0
nifgs=0 & cal_nifgs=0 & freq_lowf=0 & gamma_lowf=0
freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0 & latcut=0 & loncut=0 & st_sub=0
gain_convg=0 & gain_iter=0 & glat=0 & glon=0 & sky_wgts=0 & cal_wgts=0
zodi_array=0 & tm=0 & sky_dihd=0 & sky_glitch=0 & sky_s0=0 & scan=0
badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lls=0 & badcoadd_lsf=0
badcal_rls=0 & badcal_rsf=0 & badcal_lls=0 & badcal_lsf=0 & cvec_mask=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0 & del_temp=0
cvec_cut=0 & n_stripes=0 & max_frac=0 & cal_wgts_ds=0 & cal_idx=0 & cal_dsp=0
cal_tm=0 & ical=0 & refh=0 & skyh=0 & dihd=0 & cal_glitch=0 & cal_s0=0
ejg_lowf=0 & c_lowf=0 & ndf_ejg_lowf=0 & chi2_lowf=0 & chi2_loch=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END

;

