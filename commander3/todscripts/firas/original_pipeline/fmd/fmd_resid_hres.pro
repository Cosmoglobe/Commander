Pro FMD_RESID_HRES,error
;
;
;  FMD_RESID_HRES drives the FMD_RESID procedure to compute pixel spectra
;  and coadd residuals for combined HRES destriped data,
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
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HR.ISS and
;                        CVECTORS_HR.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing HRES.ISS, HRES_WEIGHTS.ISS,
;                        and HRES_WEIGHTS_F.ISS .
;
;    CSDR$FIRAS_OUT   =  Directory containing HRES_EJG.ISS and
;                        HRES_DSP.ISS, and where IDL save set
;                        HRES_SKY_RESID.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_RESID_HRES,error
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
 PRINT,'FMD_RESID_HRES : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_RESID_HRES,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_RESID_HRES,error'
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
RESTORE,'csdr$firas_ref:fmd_quals_hr.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HRES_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:hres_weights.iss'
;
IF (freq_dependence EQ 'Y') THEN BEGIN
 PRINT,'Restoring ' + intrans(0) + 'HRES_WEIGHTS_F.ISS'
 RESTORE,'csdr$firas_in:hres_weights_f.iss'
ENDIF
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'HRES_DSP.ISS'
RESTORE,'csdr$firas_out:hres_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_hr GE freq_band(0))and(f_hr LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hr)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_RESID_HRES : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;
nf2 = WHERE((freq_hres GE freq_band(0))and(freq_hres LE freq_band(1)),cf2)
IF (cf2 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'FMD_RESID_HRES : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Frequency Cut
; -------------
freq = f_hr(nf)
;

IF (freq_dependence EQ 'Y') THEN BEGIN

 ; FRAC_FWGT Frequency Cut
 ; -----------------------
 frac_fwgt = frac_fwgt(nf2,*)
 ;

 ; Restore Individual C-Vectors
 ; ----------------------------
 RESTORE,'csdr$firas_ref:cvectors_hr.iss'
 ;

 ; Are Frequencies Compatible ?
 ; ----------------------------
 nf = WHERE((f_hr GE freq_band(0))and(f_hr LE freq_band(1)),cf)
 IF (cf NE N_ELEMENTS(freq)) THEN BEGIN
  PRINT,'FMD_RESID_HRES : C-Vector Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 aa = MAX(ABS(freq-f_hr(nf)))
 IF (aa NE 0.) THEN BEGIN
  PRINT,'FMD_RESID_HRES : C-Vector Frequencies are NOT Compatible !'
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
PRINT,'Restoring ' + intrans(0) + 'HRES.ISS'
RESTORE,'csdr$firas_in:hres.iss'
;

; Restore HRES_EJG Save Set
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'HRES_EJG.ISS'
RESTORE,'csdr$firas_out:hres_ejg.iss'
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
              px_spec=sp_hres,resid=sky_resid_hres)
;

IF (freq_dependence NE 'Y') THEN error = $
    FMD_RESID(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,sky_idx=sky_idx,$
              mask=dirbe_mask,w_frac=frac_wgt,$
              px_spec=sp_hres,resid=sky_resid_hres)
;

; Make HRES_SKY_RESID.ISS Save Set
; ----------------------------------
chanscan = 'HRES'
chan_label = 'LLF_RLF'
freq_hres = f_hr
;
sname = 'csdr$firas_out:hres_sky_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,$
                    sky_idx,sky_wgts_ds,frac_wgt,sky_mask,freq_hres,$
                    sky_resid_hres,sp_hres
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_SKY_RESID.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0 & freq_hres=0 & gamma_hres=0
fsl_idx=0 & sky_lbl=0 & cal_lbl=0
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
ejg_hres=0 & c_hres=0 & ndf_ejg_hres=0 & chi2_hres=0 & chi2_hrch=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END

;

