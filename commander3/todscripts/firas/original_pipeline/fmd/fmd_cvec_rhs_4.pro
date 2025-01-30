Pro FMD_CVEC_RHS_4,error
;

;
;  FMD_CVEC_RHS_4 drives the FMD_CVECTOR procedure to compute
;  the Band_4 RHSS C-Vector .
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_CVECTOR
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_4.ISS
;
;    CSDR$FIRAS_IN    =  Directory containing RHS.ISS
;
;    CSDR$FIRAS_OUT   =  Directory containing RHS_4_DSP.ISS, RHS_WEIGHTS.ISS,
;                        RHS_4_DVECTOR.ISS, and RHS_4_EJG.ISS, and where
;                        RHS_4_CVECTOR_2.ISS will be sent .
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_CVEC_RHS_4,error
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
 PRINT,'FMD_CVEC_RHS_4 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_CVEC_RHS_4,error'
 PRINT,' '
 PRINT,'Returning with Error.'
 PRINT,' '
 PRINT,'Try again with valid invocation.'
 PRINT,' '
 PRINT,'Example : IDL> FMD_CVEC_RHS_4,error'
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
RESTORE,'csdr$firas_ref:fmd_quals_hi_4.iss'
;

; Restore IDL Save Set of Destriped Spectra
; -----------------------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RHS_4_DSP.ISS'
RESTORE,'csdr$firas_out:rhs_4_dsp.iss'
;

; Frequency Array
; ---------------
nf = WHERE((f_hi4 GE freq_band(0))and(f_hi4 LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hi4)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_CVEC_RHS_4 : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
freq = f_hi4(nf)
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS.ISS'
RESTORE,'csdr$firas_in:rhs.iss'
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'RHS_WEIGHTS.ISS'
RESTORE,'csdr$firas_ref:rhs_weights.iss'
;

; Restore RHS_4_EJG Save Set
; --------------------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RHS_4_EJG.ISS'
RESTORE,'csdr$firas_out:rhs_4_ejg.iss'
;

; Restore D-Vector
; ----------------
PRINT,' '
PRINT,'Restoring ' + outtrans(0) + 'RHS_4_DVECTOR.ISS'
RESTORE,'csdr$firas_out:rhs_4_dvector.iss'
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

; Compute C-Vector
; ----------------
cvec_rhs_4 = $
 FMD_CVECTOR(pixel=px,spec=dsp,sky_wgts_ds=sky_wgts_ds,w_frac=frac_wgt,$
             dmask=dirbe_mask,mask=cvec_mask,n_stripes=n_stripes,$
             max_frac=max_frac,good_dvec=good_dvec_rhs_4,$
             ndf_cvec=ndf_cvec_rhs_4,good_cvec=good_cvec_rhs_4)
;

; Make RHS_4_CVECTOR_2.ISS Save Set
; ---------------------------------
chanscan = 'RHS_4'
freq_dependence = 'N'
freq_rhs_4 = f_hi4
sname = 'csdr$firas_out:rhs_4_cvector_2.iss'
SAVE,filename=sname,chanscan,stripes_descrip,freq_dependence,cvec_cut,$
                    freq_rhs_4,cvec_rhs_4,ndf_cvec_rhs_4,good_cvec_rhs_4
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'RHS_4_CVECTOR_2.ISS" Created.'
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
del_temp=0 & cal_dsp=0 & tm=0 & st_sub=0 & cal_tm=0 & ical=0 & refh=0 & skyh=0
dihd=0 & sky_glitch=0 & cal_glitch=0 & fsl_idx=0 & cal_wgts=0 & sky_dihd=0
sky_s0=0 & cal_s0=0 & cal_wgts=0 & cal_wgts_ds=0 & chan_label=0 & dvec_cut=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
time=0 & model_var=0 & lp_orderv=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
ejg_rhs_4=0 & c_rhs_4=0 & ndf_ejg_rhs_4=0 & chi2_rhs_4=0
n_rhs=0 & b_rhs=0 & l_rhs=0 & d_rhs=0 & dvec_rhs_4=0 & ndf_dvec_rhs_4=0
zodi_array=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0 & sky_lbl=0 & scan=0
time=0 & model_var=0 & lp_orderv=0
tm_up=0 & tm_dn=0 & gpow=0 & bpow=0 & nscan=0 & delt=0 & tmcut=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
