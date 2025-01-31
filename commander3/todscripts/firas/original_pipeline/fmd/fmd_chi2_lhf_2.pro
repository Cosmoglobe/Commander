Pro FMD_CHI2_LHF_2,error
;

;
;  FMD_CHI2_LHF_2 drives the FMD_CHI2 and FMD_CAL_CHI2 procedures to
;  compute C-Vector, pixel spectra, coadd residuals, and coadd chi-squared
;  for LHFA Band_2 offset destriped data, using frequency-independent
;  coadd weights.
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
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS
;
;    CSDR$FIRAS_IN    =  Directory containing LHF.ISS
;
;    CSDR$FIRAS_OUT   =  Directory containing LHF_2_DSP.ISS, LHF_WEIGHTS.ISS,
;                        and LHF_2_EJG.ISS, and where output save sets 
;                        LHF_2_CVECTOR.ISS, LHF_2_SKYMAP_2.ISS, LHF_2_SKY_CHI2.ISS
;                        and LHF_2_CAL_CHI2.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CHI2_LHF_2,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-97.
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
 PRINT,'FMD_CHI2_LHF_2 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CHI2_LHF_2,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CHI2_LHF_2,error'
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
RESTORE,'csdr$firas_ref:fmd_quals_hi_2.iss'
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LHF_2_DSP.ISS'
RESTORE,'csdr$firas_out:lhf_2_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_hi2 GE freq_band(0))and(f_hi2 LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hi2)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CHI2_LHF_2 : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
freq = f_hi2(nf)
;

; Restore LHF_2 ZODI Spectra
; --------------------------
RESTORE,'csdr$firas_in:lhf_2_zodispec.iss'
nq = WHERE(f_hi2 NE freq,cq)
IF (cq GT 0) THEN BEGIN
 PRINT,' '
 PRINT,'FMD_CHI2_LHF_2 : Zodi Frequency Array is Mismatched !'
 PRINT,' '
 RETURN
ENDIF
;

; Subtract ZODI from Destriped Sky Spectra
; ----------------------------------------
dsp = dsp - zodi_spec
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'LHF.ISS'
RESTORE,'csdr$firas_in:lhf.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'LHF_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:lhf_weights.iss'
;

; Restore LHF_2_EJG Save Set
; -------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'LHF_2_EJG.ISS'
RESTORE,'csdr$firas_out:lhf_2_ejg.iss'
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
cvec_lhf_2 = $
 FMD_CHI2(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,w_frac=frac_wgt,$
          dmask=dirbe_mask,mask=cvec_mask,n_stripes=n_stripes,$
          max_frac=max_frac,px_spec=sp_lhf_2,resid=sky_resid_lhf_2,$
          chi2=sky_chi2_lhf_2,ndf_cvec=ndf_cvec_lhf_2)
;
error = $
 FMD_CAL_CHI2(xcal=xcal,del_temp=del_temp,freq=f_hi2,cal_wgts=cal_wgts,$
              cal_spec=cal_dsp,cvec=cvec_lhf_2,resid=cal_resid_lhf_2,$
              chi2=cal_chi2_lhf_2)
;

; ZODI Pixel Spectra
; ------------------
good = WHERE(sky_wgts_ds GT 0.)
wg = sky_wgts_ds(good)
pxg = px(good)
;
dummy = FMD_PIX_SUM(pixel=pxg,data=zodi_spec(0,good),wgt=wg,cmp_px=cmp_px)
dummy = FLTARR(cf,N_ELEMENTS(cmp_px))
;
FOR i=0,cf-1 DO $
 dummy(i,*) = FMD_PIX_SUM(pixel=pxg,data=zodi_spec(i,good),wgt=wg,cmp_px=cmp_px)
;
PIX2XY,cmp_px,data=dummy,res=6,/six,ras=z_lhf_2
;

; Restore ZODI to skymap
; ----------------------
sp_lhf_2 = sp_lhf_2 + z_lhf_2
;

; Make IDL Save Sets
; ------------------
;
chanscan = 'LHF_2'
freq_dependence = 'N'
freq_lhf_2 = f_hi2
;

; Make LHF_2_SKYMAP_2.ISS Save Set
; ------------------------------
sname = 'csdr$firas_out:lhf_2_skymap_2.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,$
                    freq_lhf_2,sp_lhf_2,cvec_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_SKYMAP_2.ISS" Created.'
PRINT,' '
;

; Make LHF_2_SKY_RESID.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:lhf_2_sky_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,freq_lhf_2,sky_resid_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_SKY_RESID.ISS" Created.'
PRINT,' '
;

; Make LHF_2_CVECTOR.ISS Save Set
; -----------------------------
sname = 'csdr$firas_out:lhf_2_cvector.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_lhf_2,cvec_lhf_2,ndf_cvec_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_CVECTOR.ISS" Created.'
PRINT,' '
;

; Make LHF_2_SKY_CHI2.ISS Save Set
; ------------------------------
sname = 'csdr$firas_out:lhf_2_sky_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,px,tm,$
                    sky_s0,sky_dihd,sky_glitch,sky_wgts_ds,frac_wgt,$
                    sky_mask,cvec_mask,freq_lhf_2,cvec_lhf_2,sky_chi2_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_SKY_CHI2.ISS" Created.'
PRINT,' '
;

; Make LHF_2_CAL_RESID.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:lhf_2_cal_resid.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_lhf_2,cal_resid_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_CAL_RESID.ISS" Created.'
PRINT,' '
;

; Make LHF_2_CAL_CHI2.ISS Save Set
; ------------------------------
sname = 'csdr$firas_out:lhf_2_cal_chi2.iss'
SAVE,filename=sname,chanscan,chan_label,stripes_descrip,freq_dependence,$
                    xcal,cal_tm,ical,skyh,refh,dihd,cal_glitch,$
                    cal_s0,cal_wgts_ds,cal_wgts,freq_lhf_2,cvec_lhf_2,$
                    cal_chi2_lhf_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'LHF_2_CAL_CHI2.ISS" Created.'
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
badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_lhf=0
badcal_rhs=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_lhf=0
pixel_wgt=0 & stripe_order=0 & stripe_id=0 & n_krnl=0 & n_cmn=0
n_lhf=0 & b_lhf=0 & l_lhf=0 & d_lhf=0
ejg_lhf_2=0 & c_lhf_2=0 & ndf_ejg_lhf_2=0 & chi2_lhf_2=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
time=0 & model_var=0 & lp_orderv=0 & fsl_idx=0 & st_sub=0 & sl_weight_rat=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
