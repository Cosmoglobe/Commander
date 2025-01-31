Pro FMD_CHI2_RLFX,error
;

;
;  FMD_CHI2_RLFX drives the FMD_CHI2 and FMD_CAL_CHI2 procedures to
;  compute C-Vector, pixel spectra, coadd residuals, and coadd chi-squared
;  for RLLF offset destriped data, using default coadd weights.
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
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HR.ISS
;
;    CSDR$FIRAS_IN    =  Directory containing RLF.ISS
;
;    CSDR$FIRAS_OUT   =  Directory containing RLFX_DSP.ISS, RLFX_WEIGHTS.ISS,
;                        and RLFX_EJG.ISS, and where output save sets 
;                        RLFX_CVECTOR.ISS, RLFX_SKY_RESID.ISS,
;                        RLFX_SKY_CHI2.ISS, RLFX_CAL_RESID.ISS,
;                        and RLFX_CAL_CHI2.ISS, will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CHI2_RLFX,error
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
 PRINT,'FMD_CHI2_RLFX : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CHI2_RLFX,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CHI2_RLFX,error'
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
RESTORE,'csdr$firas_ref:fmd_quals_hr.iss'
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RLFX_DSP.ISS'
RESTORE,'csdr$firas_out:rlfx_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_hr GE freq_band(0))and(f_hr LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hr)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CHI2_RLFX : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
freq = f_hr(nf)
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RLF.ISS'
RESTORE,'csdr$firas_in:rlf.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'RLFX_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:rlfx_weights.iss'
;

; Restore RLFX_EJG Save Set
; ---------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RLFX_EJG.ISS'
RESTORE,'csdr$firas_out:rlfx_ejg.iss'
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
cvec_rlf = $
 FMD_CHI2(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,w_frac=frac_wgt,$
          dmask=dirbe_mask,mask=cvec_mask,n_stripes=n_stripes,$
          max_frac=max_frac,px_spec=sp_rlf,resid=sky_resid_rlf,$
          chi2=sky_chi2_rlf,ndf_cvec=ndf_cvec_rlf)
;
error = $
 FMD_CAL_CHI2(xcal=xcal,del_temp=del_temp,freq=f_hr,cal_wgts=cal_wgts,$
              cal_spec=cal_dsp,cvec=cvec_rlf,resid=cal_resid_rlf,$
              chi2=cal_chi2_rlf)
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'RLF'
freq_dependence = 'N'
freq_rlf = f_hr
;

; Make RLFX_SKY_RESID.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_out:rlfx_sky_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,freq_rlf,sky_resid_rlf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RLFX_SKY_RESID.ISS" Created.'
PRINT,' '
;

; Make RLFX_CVECTOR.ISS Save Set
; --------------------------------
sname = 'csdr$firas_out:rlfx_cvector.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_rlf,cvec_rlf,ndf_cvec_rlf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RLFX_CVECTOR.ISS" Created.'
PRINT,' '
;

; Make RLFX_SKY_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:rlfx_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,cvec_mask,freq_rlf,cvec_rlf,sky_chi2_rlf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RLFX_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make RLFX_CAL_RESID.ISS Save Set
; ----------------------------------
sname = 'csdr$firas_out:rlfx_cal_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_rlf,cal_resid_rlf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RLFX_CAL_RESID.ISS" Created.'
PRINT,' '
;

; Make RLFX_CAL_CHI2.ISS Save Set
; ---------------------------------
sname = 'csdr$firas_out:rlfx_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_rlf,cvec_rlf,$
                    cal_chi2_rlf
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RLFX_CAL_CHI2.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
nifgs=0 & cal_nifgs=0
freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & lp_order=0 & xtm=0 & cmn_tm=0
step_up=0 & step_dn=0 & ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0
ref_dihd=0 & dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & good_sky_dihd=0
tmin=0 & tmax=0 & good_cal_dihd=0 & hot_cal=0 & latcut=0 & loncut=0
sky_wgts=0 & gain_convg=0 & gain_iter=0 & glat=0 & glon=0 & st_sub=0
badcoadd_rlf=0 & badcoadd_llf=0 & badcal_rlf=0 & badcal_llf=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0
ejg_rlf=0 & c_rlf=0 & ndf_ejg_rlf=0 & chi2_rlf=0
n_rlf=0 & b_rlf=0 & l_rlf=0 & d_rlf=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
time=0 & model_var=0 & lp_orderv=0 & fsl_idx=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
