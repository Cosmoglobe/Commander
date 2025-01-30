Pro FMD_CHI2_LSFX,error
;

;
;  FMD_CHI2_LSFX drives the FMD_CHI2 and FMD_CAL_CHI2 procedures to
;  compute C-Vector, pixel spectra, coadd residuals, and coadd chi-squared
;  for LLFA offset destriped data, using default coadd weights.
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
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_LO.ISS
;
;    CSDR$FIRAS_IN    =  Directory containing LSF.ISS
;
;    CSDR$FIRAS_OUT   =  Directory containing LSFX_DSP.ISS, LSFX_WEIGHTS.ISS,
;                        and LSFX_EJG.ISS, and where output save sets 
;                        LSFX_CVECTOR.ISS, LSFX_SKY_RESID.ISS,
;                        LSFX_SKY_CHI2.ISS, LSFX_CAL_RESID.ISS,
;                        and LSFX_CAL_CHI2.ISS, will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CHI2_LSFX,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 17-Jun-97.
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
 PRINT,'FMD_CHI2_LSFX : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CHI2_LSFX,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CHI2_LSFX,error'
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

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LSFX_DSP.ISS'
RESTORE,'csdr$firas_out:lsfx_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_lo GE freq_band(0))and(f_lo LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_lo)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CHI2_LSFX : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
freq = f_lo(nf)
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
PRINT,'Restoring ' + reftrans(0) + 'LSFX_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:lsfx_weights.iss'
;

; Restore LSFX_EJG Save Set
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LSFX_EJG.ISS'
RESTORE,'csdr$firas_out:lsfx_ejg.iss'
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
cvec_lsf = $
 FMD_CHI2(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,w_frac=frac_wgt,$
          dmask=dirbe_mask,mask=cvec_mask,n_stripes=n_stripes,$
          max_frac=max_frac,px_spec=sp_lsf,resid=sky_resid_lsf,$
          chi2=sky_chi2_lsf,ndf_cvec=ndf_cvec_lsf)
;
error = $
 FMD_CAL_CHI2(xcal=xcal,del_temp=del_temp,freq=f_lo,cal_wgts=cal_wgts,$
              cal_spec=cal_dsp,cvec=cvec_lsf,resid=cal_resid_lsf,$
              chi2=cal_chi2_lsf)
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'LSF'
freq_dependence = 'N'
freq_lsf = f_lo
;

; Make LSFX_SKY_RESID.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_out:lsfx_sky_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,freq_lsf,sky_resid_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSFX_SKY_RESID.ISS" Created.'
PRINT,' '
;

; Make LSFX_CVECTOR.ISS Save Set
; --------------------------------
sname = 'csdr$firas_out:lsfx_cvector.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_lsf,cvec_lsf,ndf_cvec_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSFX_CVECTOR.ISS" Created.'
PRINT,' '
;

; Make LSFX_SKY_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:lsfx_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,cvec_mask,freq_lsf,cvec_lsf,sky_chi2_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSFX_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make LSFX_CAL_RESID.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_out:lsfx_cal_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_lsf,cal_resid_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSFX_CAL_RESID.ISS" Created.'
PRINT,' '
;

; Make LSFX_CAL_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:lsfx_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_lsf,cvec_lsf,cal_chi2_lsf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LSFX_CAL_CHI2.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0
freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0 & latcut=0 & loncut=0
sky_wgts=0 & gain_convg=0 & gain_iter=0 & glat=0 & glon=0
badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lls=0 & badcoadd_lsf=0
badcal_rls=0 & badcal_rsf=0 & badcal_lls=0 & badcal_lsf=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0
ejg_lsf=0 & c_lsf=0 & ndf_ejg_lsf=0 & chi2_lsf=0
n_lsf=0 & b_lsf=0 & l_lsf=0 & d_lsf=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
time=0 & model_var=0 & lp_orderv=0 & sl_weight_rat=0 & fsl_idx=0 & st_sub=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
