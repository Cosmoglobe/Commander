Pro FMD_CHI2_HIGH_4,error
;
;  FMD_CHI2_HIGH_4 drives the FMD_CHI2 and FMD_CAL_CHI2 procedures to
;  compute C-Vector, pixel spectra, coadd residuals, and coadd chi-squared
;  for combined high channel Band_4 destriped data, using frequency
;  dependent coadd weights.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_CHI2
;                        FMD_CAL_CHI2
;                        COORCONV
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_4.ISS and
;                        CVECTORS_HI.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing HIGH_4_ZODI_SKYMAP.ISS, HIGH.ISS,
;                        HIGH_WEIGHTS.ISS, and HIGH_WEIGHTS_F.ISS .
;
;    CSDR$FIRAS_OUT   =  Directory containing HIGH_4_EJG.ISS and
;                        HIGH_4_DSP.ISS, and where IDL save sets
;                        HIGH_4_CVECTOR.ISS, HIGH_4_SKYMAP_F.ISS
;                        or HIGH_4_SKYMAP_2.ISS, HIGH_4_SKY_CHI2.ISS,
;                        HIGH_4_CAL_RESID.ISS, and HIGH_4_CAL_CHI2.ISS
;                        will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CHI2_HIGH_4,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 19-Jun-1997
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
 PRINT,'FMD_CHI2_HIGH_4 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CHI2_HIGH_4,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CHI2_HIGH_4,error'
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
freq_flag = 'Y'
;

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hi_4.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:high_weights.iss'
;
IF (freq_flag EQ 'Y') THEN BEGIN
 PRINT,'Restoring ' + intrans(0) + 'HIGH_WEIGHTS_F.ISS'
 RESTORE,'csdr$firas_in:high_weights_f.iss'
ENDIF
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'HIGH_4_DSP.ISS'
RESTORE,'csdr$firas_out:high_4_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_hi4 GE freq_band(0))and(f_hi4 LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hi4)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CHI2_HIGH_4 : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;
nf2 = WHERE((freq_high GE freq_band(0))and(freq_high LE freq_band(1)),cf2)
IF (cf2 NE cf) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CHI2_HIGH_4 : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Frequency Cut
; -------------
freq = f_hi4(nf)
;

IF (freq_flag EQ 'Y') THEN BEGIN

 ; FRAC_FWGT Frequency Cut
 ; -----------------------
 frac_fwgt = frac_fwgt(nf2,*)
 ;

 ; Restore Individual C-Vectors
 ; ----------------------------
 RESTORE,'csdr$firas_ref:cvectors_hi.iss'
 ;

 ; Are Frequencies Compatible ?
 ; ----------------------------
 nf = WHERE((f_hi GE freq_band(0))and(f_hi LE freq_band(1)),cf)
 IF (cf NE N_ELEMENTS(freq)) THEN BEGIN
  PRINT,'FMD_CHI2_HIGH_4 : C-Vector Frequencies are NOT Compatible !'
  error = 1
  RETURN
 ENDIF
 ;
 aa = MAX(ABS(freq-f_hi(nf)))
 IF (aa NE 0.) THEN BEGIN
  PRINT,'FMD_CHI2_HIGH_4 : C-Vector Frequencies are NOT Compatible !'
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

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH.ISS'
RESTORE,'csdr$firas_in:high.iss'
;

; Restore ZODI Skymap Spectra
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH_4_ZODI_SKYMAP.ISS'
RESTORE,'csdr$firas_in:high_4_zodi_skymap.iss'
;

; Restore HIGH_4_EJG Save Set
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'HIGH_4_EJG.ISS'
RESTORE,'csdr$firas_out:high_4_ejg.iss'
;

n_sky = N_ELEMENTS(px)    ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)  ; Number of CAL Coadds
;

nmask = WHERE(dirbe_mask EQ 0,cmask)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cmask)) + ' Coadds Masked by DIRBE_CUT'
PRINT,' '
;

nmask = WHERE(sky_mask EQ 0,cmask)
PRINT,STRCOMPRESS(STRING(cmask)) + ' Total Coadds Masked from DESTRIPER'
PRINT,''
;

; C-Vector Mask (CVEC_CUT)
; ------------------------
cvec_mask = sky_mask
good_cut = 'Y'
IF ((cvec_cut(0) LE 0)or(cvec_cut(0) GT 90)) THEN good_cut = 'N'
IF ((cvec_cut(1) LE 0)or(cvec_cut(1) GT 180)) THEN good_cut = 'N'
;
IF (good_cut EQ 'Y') THEN BEGIN
;
 ll = COORCONV(px,infmt='p',inco='f',outfmt='l',outco='g')
 pix_glat = ll(*,1)
 pix_glon = ll(*,0)
 nq = WHERE(pix_glon GT 180.,cq)
 IF (cq GT 0) THEN pix_glon(nq) = pix_glon(nq) - 360.
 nq = WHERE(pix_glon LT -180.,cq)
 IF (cq GT 0) THEN pix_glon(nq) = pix_glon(nq) + 360.
;
 bad_cvec = $
   WHERE((ABS(pix_glat) LE cvec_cut(0))and(ABS(pix_glon) LE cvec_cut(1)),cbad)
 IF (cbad GT 0) THEN cvec_mask(bad_cvec) = 0
;
ENDIF
;
bad_cvec = WHERE(cvec_mask NE 1,cbad)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cbad)) + ' Coadds Masked from C-VECTOR'
PRINT,' '
;

; Compute Pixel Spectra, Coadd Residuals, C-Vector, and Coadd Chi-Squared
; -----------------------------------------------------------------------
;
IF (freq_flag EQ 'Y') THEN BEGIN
 ;
 cvec_high_4 = $
  FMD_CHI2(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,sky_idx=sky_idx,$
           dmask=dirbe_mask,mask=cvec_mask,wfac=wfac,w_frac=frac_wgt,$
           wf_frac=frac_fwgt,n_stripes=n_stripes,max_frac=max_frac,$
           px_spec=sp_high_4,chi2=sky_chi2_high_4,$
           ndf_cvec=ndf_cvec_high_4)
 ;

 error = $
  FMD_CAL_CHI2(xcal=xcal,del_temp=del_temp,freq=f_hi4,cal_idx=cal_idx,$
               cal_wgts=cal_wgts_ds,cal_spec=cal_dsp,wfac=wfac,$
               cvec=cvec_high_4,resid=cal_resid_high4,chi2=cal_chi2_high_4)
 ;

ENDIF
;

IF (freq_flag NE 'Y') THEN BEGIN
 ;
 cvec_high_4 = $
  FMD_CHI2(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,$
           dmask=dirbe_mask,mask=cvec_mask,w_frac=frac_wgt,$
           n_stripes=n_stripes,max_frac=max_frac,$
           px_spec=sp_high_4,chi2=sky_chi2_high_4,$
           ndf_cvec=ndf_cvec_high_4)
 ;

 error = $
  FMD_CAL_CHI2(xcal=xcal,del_temp=del_temp,freq=f_hi4,cal_spec=cal_dsp,$
               cal_wgts=cal_wgts_ds,resid=cal_resid_high4,$
               chi2=cal_chi2_high_4,cvec=cvec_high_4)
 ;

ENDIF
;

; Add ZODI to Skymap Spectra
; --------------------------
sp_high_4 = sp_high_4 + z_high_4
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'HIGH_4'
chan_label = 'LHS_RHS_LHF_RHF'
freq_dependence = freq_flag
freq_high_4 = f_hi4
;

IF (freq_flag EQ 'Y') THEN BEGIN
 ;
 ; Make HIGH_4_SKYMAP_F.ISS Save Set
 ; ---------------------------------
 sname = 'csdr$firas_out:high_4_skymap_f.iss'
 SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                     freq_high_4,sp_high_4,z_high_4,cvec_high_4
 PRINT,' '
 PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_SKYMAP_F.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (freq_flag NE 'Y') THEN BEGIN
 ;
 ; Make HIGH_4_SKYMAP_2.ISS Save Set
 ; ---------------------------------
 sname = 'csdr$firas_out:high_4_skymap_2.iss'
 SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                     freq_high_4,sp_high_4,z_high_4,cvec_high_4
 PRINT,' '
 PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_SKYMAP_2.ISS" Created.'
 PRINT,' '
;
ENDIF
;

; Make HIGH_4_CVECTOR.ISS Save Set
; --------------------------------
sname = 'csdr$firas_out:high_4_cvector.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_high_4,cvec_high_4,ndf_cvec_high_4
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_CVECTOR.ISS" Created.'
PRINT,' '
;

; Make HIGH_4_SKY_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:high_4_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,$
                    sky_idx,sky_wgts_ds,frac_wgt,sky_mask,freq_high_4,$
                    cvec_high_4,sky_chi2_high_4
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make HIGH_4_CAL_RESID.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_out:high_4_cal_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,cal_idx,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_high_4,cvec_high_4,$
                    cal_resid_high4
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_CAL_RESID.ISS" Created.'
PRINT,' '
;

; Make HIGH_4_CAL_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:high_4_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,cal_idx,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_high_4,cvec_high_4,$
                    cal_chi2_high_4
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_4_CAL_CHI2.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0 & freq_high=0 & gamma_high_4=0
st_sub=0 & fsl_idx=0 & sky_lbl=0 & cal_lbl=0
freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0 & latcut=0 & loncut=0
gain_convg=0 & gain_iter=0 & glat=0 & glon=0 & sky_wgts=0
zodi_array=0 & tm=0 & sky_dihd=0 & sky_glitch=0 & sky_s0=0 & scan=0 & fz_hi4=0
badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_lhf=0
badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0
ejg_high_4=0 & c_high_4=0 & ndf_ejg_high_4=0 & chi2_high_4=0 & chi2_hi4ch=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
