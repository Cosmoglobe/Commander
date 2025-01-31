Pro FMD_COVAR_LOWF,error
;

;
;  FMD_COVAR_LOWF drives the FMD_COVAR procedure to create IDL
;  save sets of combined low-channel covariance matrices.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_LO.ISS and
;                        CVECTORS_LO.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing LOWF_WEIGHTS.ISS.
;
;    CSDR$FIRAS_OUT   =  Directory containing LOWF_SKY_RESID.ISS, and
;                        where output LOWF_COVAR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_COVAR_LOWF,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 29-Apr-97.
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
 PRINT,'FMD_COVAR_LOWF : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_COVAR_LOWF,error'
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
PRINT,'Restoring ' + reftrans(0) + 'FMD_QUALS_LO.ISS'
RESTORE,'csdr$firas_ref:fmd_quals_lo.iss'
nst = n_stripes(1)
;

; Restore Coadd Residuals
; -----------------------
PRINT,'Restoring LOWF_SKY_RESID.ISS'
RESTORE,'csdr$firas_out:lowf_sky_resid.iss'
freq = freq_lowf
;
nfreq = N_ELEMENTS(freq)
;

; Restore Coadd Weights
; ---------------------
PRINT,'Restoring LOWF_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:lowf_weights.iss'
;

; Frequency-Dependent Weight Factor
; ---------------------------------
IF (freq_dependence EQ 'Y') THEN BEGIN
 ;

 ; Restore Individual C-Vectors
 ; ----------------------------
 RESTORE,'csdr$firas_ref:cvectors_lo.iss'
 ;

 ; Are Frequencies Compatible ?
 ; ----------------------------
 nf2 = WHERE((freq_lowf GE freq_band(0))and(freq_lowf LE freq_band(1)),cf2)
 IF (cf2 NE N_ELEMENTS(freq)) THEN BEGIN
  PRINT,'FMD_COVAR_LOWF : C-Matrix Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
 IF (cf NE N_ELEMENTS(freq)) THEN BEGIN
  PRINT,'FMD_COVAR_LOWF : C-Matrix Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 aa = MAX(ABS(freq-f_lo(nf)))
 IF (aa NE 0.) THEN BEGIN
  PRINT,'FMD_COVAR_LOWF : C-Matrix Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;

 ; Apply Frequency Cut to C-Vectors
 ; --------------------------------
 cvec_lhs = cvec_lhs(nf)
 cvec_rhs = cvec_rhs(nf)
 cvec_lhf = cvec_lhf(nf)
 cvec_rhf = cvec_rhf(nf)
 ;

 ; Concatenate C-Vectors
 ; ---------------------
 cvecs = [[cvec_lhs],[cvec_rhs],[cvec_lhf],[cvec_rhf]]
 ;

 ; Restore Frequency-Dependent Pixel Fractional Weights
 ; ----------------------------------------------------
 PRINT,'Restoring LOWF_WEIGHTS_F.ISS'
 RESTORE,'csdr$firas_out:lowf_weights_f.iss'
 ;
 ; FRAC_FWGT Frequency Cut
 ; -----------------------
 frac_fwgt = frac_fwgt(nf2,*)
 ;

 ; Frequency Dependent Weight Adjustment
 ; -------------------------------------
 ;
 ; Normalization defined so that RHSS weights are unaffected
 ; ---------------------------------------------------------
 wfac = 0.*cvecs + 1.
 fac0 = cvecs(*,1)^2. * scale_factor(1)
 wfac(*,0) = fac0 / cvecs(*,0)^2. / scale_factor(0)
 wfac(*,2) = fac0 / cvecs(*,2)^2. / scale_factor(2)
 wfac(*,3) = fac0 / cvecs(*,3)^2. / scale_factor(3)
 ;

ENDIF
;

; Compute C-Matrices
; ------------------
;
PRINT,'Calling FMD_COVAR'

IF (freq_dependence NE 'Y') THEN error = $
    FMD_COVAR(nfreq=nfreq,resid=sky_resid_lowf,weight=sky_wgts_ds,mask=sky_mask,$
              cvec_cut=cvec_cut,px=px,n_stripes=nst,max_frac=max_frac,$
              w_frac=frac_wgt,covar=covar_lowf,ndf_covar=ndf_covar_lowf)
;

IF (freq_dependence EQ 'Y') THEN error = $
    FMD_COVAR(nfreq=nfreq,resid=sky_resid_lowf,weight=sky_wgts_ds,mask=sky_mask,$
              cvec_cut=cvec_cut,px=px,n_stripes=nst,max_frac=max_frac,$
              w_frac=frac_wgt,sky_idx=sky_idx,wfac=wfac,wf_frac=frac_fwgt,$
              covar=covar_lowf,ndf_covar=ndf_covar_lowf)
;

IF (error NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_LOWF : Error Returned from FMD_COVAR !'
 RETURN
ENDIF
;

; Create LOWF_COVAR Save Set
; --------------------------
chanscan = 'LOWF'
sname = 'csdr$firas_out:lowf_covar.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,n_stripes,max_frac,$
                    freq_lowf,covar_lowf,ndf_covar_lowf
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LOWF_COVAR.ISS" Created.'
PRINT,' '
;

del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0
cmn_tm=0 & step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0
zodi_array=0 & good_sky_dihd=0 & tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0
latcut=0 & loncut=0 & gain_convg=0 & gain_iter=0 & chan_label=0 & pixel_wgt=0
stripes_descrip=0 & cal_wgts_ds=0 & cal_idx=0
badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_lhf=0
badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
sp_lowf=0 & z_lowf=0
;

RETURN
END
