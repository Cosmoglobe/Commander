Pro FMD_COVAR_HRES2,error
;

;
;  FMD_COVAR_HRES2 drives the FMD_COVAR procedure to create IDL
;  save sets of combined HRES "Band_2" covariance matrices.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         : Return Error Status
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HR.ISS and
;                        CVECTORS_HR.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing HRES_WEIGHTS.ISS.
;
;    CSDR$FIRAS_OUT   =  Directory containing HRES_SKY_RESID.ISS, and
;                        where output HRES_2_COVAR.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_COVAR_HRES2,error
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
 PRINT,'FMD_COVAR_HRES2 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_COVAR_HRES2,error'
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
PRINT,'Restoring ' + reftrans(0) + 'FMD_QUALS_HR.ISS'
RESTORE,'csdr$firas_ref:fmd_quals_hr.iss'
nst = n_stripes(1)
;

; Restore Coadd Residuals
; -----------------------
PRINT,'Restoring HRES_SKY_RESID.ISS'
RESTORE,'csdr$firas_out:hres_sky_resid.iss'
;
nf = WHERE((freq_hres GE freq_band(0))and(freq_hres LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(freq_hres)) THEN BEGIN
 PRINT,'FMD_COVAR_HRES1 : Error in Frequency Array !'
 error = 1
 RETURN
ENDIF
;

; Frequency Cut  (Band_2)
; -----------------------
fidx = INDGEN(cf/3) + cf/3
;
freq = freq_hres(fidx)
nfreq = N_ELEMENTS(freq)
sky_resid = TEMPORARY(sky_resid_hres)
sky_resid = sky_resid(fidx,*)
;

; Restore Coadd Weights
; ---------------------
PRINT,'Restoring HRES_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:hres_weights.iss'
;

; Frequency-Dependent Weight Factor
; ---------------------------------
IF (freq_dependence EQ 'Y') THEN BEGIN
 ;

 ; Restore Individual C-Vectors
 ; ----------------------------
 RESTORE,'csdr$firas_ref:cvectors_hr.iss'
 ;

 ; Are Frequencies Compatible ?
 ; ----------------------------
 nf2 = WHERE((f_hr GE freq_band(0))and(f_hr LE freq_band(1)),cf2)
 IF (cf2 NE N_ELEMENTS(freq_hres)) THEN BEGIN
  PRINT,'FMD_COVAR_HRES2 : C-Matrix Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 aa = MAX(ABS(freq-f_hr(nf2(fidx))))
 IF (aa NE 0.) THEN BEGIN
  PRINT,'FMD_COVAR_HRES2 : C-Matrix Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;

 ; Apply Frequency Cut to C-Vectors
 ; --------------------------------
 cvec_llf = cvec_llf(nf2(fidx))
 cvec_rlf = cvec_rlf(nf2(fidx))
 ;

 ; Concatenate C-Vectors
 ; ---------------------
 cvecs = [[cvec_llf],[cvec_rlf]]
 ;

 ; Restore Frequency-Dependent Pixel Fractional Weights
 ; ----------------------------------------------------
 PRINT,'Restoring HRES_WEIGHTS_F.ISS'
 RESTORE,'csdr$firas_out:hres_weights_f.iss'
 ;
 ; FRAC_FWGT Frequency Cut
 ; -----------------------
 frac_fwgt = frac_fwgt(nf2(fidx),*)
 ;

 ; Frequency Dependent Weight Adjustment
 ; -------------------------------------
 ;
 wfac = 0.*cvecs + 1.
 fac0 = cvecs(*,0)^2. * scale_factor(0)
 wfac(*,1) = fac0 / cvecs(*,1)^2. / scale_factor(1)
 ;

ENDIF
;

; Compute C-Matrices
; ------------------
;
PRINT,'Calling FMD_COVAR'

IF (freq_dependence NE 'Y') THEN error = $
    FMD_COVAR(nfreq=nfreq,resid=sky_resid,weight=sky_wgts_ds,mask=sky_mask,$
              cvec_cut=cvec_cut,px=px,n_stripes=nst,max_frac=max_frac,$
              w_frac=frac_wgt,covar=covar_hres_2,ndf_covar=ndf_covar_hres2)
;

IF (freq_dependence EQ 'Y') THEN error = $
    FMD_COVAR(nfreq=nfreq,resid=sky_resid,weight=sky_wgts_ds,mask=sky_mask,$
              cvec_cut=cvec_cut,px=px,n_stripes=nst,max_frac=max_frac,$
              w_frac=frac_wgt,sky_idx=sky_idx,wfac=wfac,wf_frac=frac_fwgt,$
              covar=covar_hres_2,ndf_covar=ndf_covar_hres2)
;

IF (error NE 0) THEN BEGIN
 PRINT,'FMD_COVAR_HRES2 : Error Returned from FMD_COVAR !'
 RETURN
ENDIF
;

; Create HIGH_COVAR Save Set
; --------------------------
chanscan = 'HRES_2'
sname = 'csdr$firas_out:hres_2_covar.iss'
SAVE,filename=sname,chanscan,freq_dependence,cvec_cut,n_stripes,max_frac,$
                    freq_hres,covar_hres_2,ndf_covar_hres2
;
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HRES_2_COVAR.ISS" Created.'
PRINT,' '
;

del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0
cmn_tm=0 & step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0
dihed_pow=0 & ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0
zodi_array=0 & good_sky_dihd=0 & tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0
latcut=0 & loncut=0 & gain_convg=0 & gain_iter=0 & chan_label=0 & pixel_wgt=0
stripes_descrip=0 & cal_wgts_ds=0 & cal_idx=0 & sp_hres=0
badcoadd_rlf=0 & badcoadd_llf=0 & badcal_rlf=0 & badcal_llf=0
;

RETURN
END
