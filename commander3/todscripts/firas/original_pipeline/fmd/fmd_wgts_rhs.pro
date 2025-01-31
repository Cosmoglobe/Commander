Pro FMD_WGTS_RHS,model_var,error
;
;
;  FMD_WGTS_RHS drives the FMD_MODEL_WGT and FMD_PIXEL_WGT procedures
;  to create IDL save sets of RHSS weights.
;
;
;  ARGUMENTS (I/O)    :
;
;   MODEL_VAR (I)     :  "Y(es)" to compute and apply a variance model for
;                        coadd weights.
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_MODEL_WGT
;                        FMD_PIXEL_WGT
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS,
;                        FMD_QUALS_HI_4.ISS, and FMD_BAD_COADD_HI.ISS,
;                        and where RHS_WEIGHTS.ISS will be sent.
;
;    CSDR$FIRAS_IN    =  Directory containing RHS.ISS and RHS_VAR.ISS
;
;    CSDR$FIRAS_OUT   =  Directory where RHS_VAR_MODEL.ISS will be sent.
;
;
;    EXAMPLE : $ uidl
;              UIDL> model_var = 'Y'
;              UIDL> FMD_WGTS_RHS,model_var,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 17-Jun-1997.
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
IF N_Params() ne 2 THEN BEGIN 
 PRINT,'FMD_WGTS_RHS : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_WGTS_RHS,model_var,error'
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

; FMD Qualifiers
; --------------
RESTORE,'csdr$firas_ref:fmd_quals_hi_2.iss'
freq1 = freq_band(0)
;
RESTORE,'csdr$firas_ref:fmd_quals_hi_4.iss'
freq2 = freq_band(1)
;
freq_band = [freq1,freq2]
;

; Restore Bad Coadd Indices
; -------------------------
RESTORE,'csdr$firas_ref:fmd_bad_coadd_hi.iss'
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'RHS.ISS'
RESTORE,'csdr$firas_in:rhs.iss'
;

; Qualifiers for Variance Model
; -----------------------------
nr = [6,100]  ;  Range of "Good" Coadd Sizes
;
gpow = [1.]   ;  Desired Powers of Glitch_Rate Functions,
              ;  (If LE 0, not enabled) .
;
bpow = [1.]   ;  Desired Powers of Bolometer Functions,
              ;  (If LE 0, not enabled) .
;
nscan=5       ;  Desired Number of Scan-Angle COS,SIN Functions,
              ;  (If LE 0, not enabled) .
;
delt = 0.     ;  Time Range (days) for Short Timescale Variance Model,
              ;  (IF LE 0., no short timescale model is applied) .
;
lp_order=[1]   ;  Order of Legendre Polynomials of Time
               ;  (IF LE 0, no Legendre Terms )
;
tm_up = [0.]   ;  Start and Stop Times of Special Time Periods
tm_dn = [0.]   ;  (IF LE 0, disabled)
;
glmax = 6.     ;  Maximum Glitch_Rate for "Good" Coadds
;
glcut = 9999.  ;  Minimum Glitch_Rate for Special High_Glitch Terms
;
tmcut = [-1.,-1.]  ;  Time Range for Special High_Glitch Terms
;
skut = [0,0,0,0,0]     ;  Controls Special Scan-Angle Correction, terms are :
                       ;  [TIME_UP,TIME_DN,SCAN_UP,SCAN_DN,Number_Smoothed]
;

; List of FMD Qualifiers
; ----------------------
PRINT,' '
PRINT,'FMD Qualifiers are :'
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,' '
PRINT,'TMIN =',tmin                     ; Minimum controllable temps
PRINT,'TMAX =',tmax                     ; Maximum controllable temps
PRINT,'GOOD_CAL_DIHD =',good_cal_dihd   ; Good Cal DIHED temps
PRINT,'GOOD_SKY_DIHD =',good_sky_dihd   ; Good Sky DIHED temps
PRINT,' '
PRINT,'BADCOADD_RHS =',badcoadd_rhs     ; Bad RHS SKY Coadds
PRINT,' '
PRINT,'BADCAL_RHS =',badcal_rhs         ; Bad RHS CAL Coadds
PRINT,' '
PRINT,' '
;

n_sky = N_ELEMENTS(px)     ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)   ; Number of CAL Coadds
;

; Build BADCOADD Index
; --------------------
badcoadd = BYTARR(n_sky)
IF (badcoadd_rhs(0) ge 0) THEN badcoadd(badcoadd_rhs) = 1B
;

; Build BADCAL Index
; --------------------
badcal = BYTARR(n_cal)
IF (badcal_rhs(0) ge 0) THEN badcal(badcal_rhs) = 1B
;

; Initialize Destriper Coadd Weights
; ----------------------------------
sky_wgts_ds = sky_wgts
cal_wgts_ds = cal_wgts
;

; De-weight Bad Sky Coadds
; ------------------------
bad = WHERE(badcoadd eq 1,cbad)
IF (cbad GT 0) THEN sky_wgts_ds(bad) = 0.
;

; Deweight Sky Spectra with Temps out-of-Range
; --------------------------------------------
bad = WHERE( sky_dihd LT good_sky_dihd(0) OR sky_dihd GT good_sky_dihd(1) )
IF (bad(0) ge 0) THEN sky_wgts_ds(bad) = 0.
;

; De-weight Bad Cal Coadds
; ------------------------
bad = WHERE(badcal eq 1,cbad)
IF (cbad GT 0) THEN cal_wgts_ds(bad) = 0.
;

; De-weight Calibration Spectra Early in Mission
; ----------------------------------------------
bad = WHERE(cal_tm lt 28.35,cbad)
IF (cbad GT 0) THEN cal_wgts_ds(bad) = 0.
;

; Deweight Calibration Spectra with Temps out-of-Range
; ----------------------------------------------------
bad = WHERE(xcal LE tmin(0) OR xcal GE tmax(0) OR $
           ical LE tmin(1)  OR ical GE tmax(1) OR $
           skyh LE tmin(2)  OR skyh GE tmax(2) OR $
           refh LE tmin(3)  OR refh GE tmax(3) OR $
           dihd LT good_cal_dihd(0) OR dihd GT good_cal_dihd(1) )
IF (bad(0) ge 0) THEN cal_wgts_ds(bad) = 0.
;


model_var = STRUPCASE(STRMID(model_var,0,1))
;
IF (model_var EQ 'Y') THEN BEGIN  ; Coadds Will Be Re-Weighted

 ; List of Model Qualifiers
 ; ------------------------
 PRINT,' '
 PRINT,'Variance Model Qualifiers are :'
 PRINT,' '
 PRINT,'LP_ORDER =',lp_order   ;  Order of Legendre Polynomials of Time
 PRINT,' '
 PRINT,'TM_UP =',tm_up         ;  Time STEP_UP
 PRINT,' '
 PRINT,'TM_DN =',tm_dn         ;  Time STEP_DN
 PRINT,' '
 PRINT,'NSCAN =',nscan         ;  Number of COS,SIN Scan Angle Functions
 PRINT,' '
 PRINT,'GPOW =',gpow           ;  Order of Glitch_Rate Polynomial
 PRINT,' '
 PRINT,'TMCUT =',tmcut         ;  Special Glitch Time Range
 PRINT,' '
 PRINT,'GLMAX =',glmax         ; Maximum Glitch-Rate
 PRINT,' '
 PRINT,'GLCUT =',glcut         ;  Minimum Special Glitch-Rate
 PRINT,' '
 PRINT,'BPOW =',bpow           ;  Order of Bolometer Polynomial
 PRINT,' '
 PRINT,'DELT =',delt           ;  Size (days) of Time Boxcar Function
 PRINT,' '
 PRINT,'SKUT =',skut           ;  Special Scan_Angle Correction Parameters
 PRINT,' '
;

 ; Variance Mask
 ; -------------
 ;
 vmask = BYTARR(N_ELEMENTS(px)) + 1B
 ;
 ll = COORCONV(px,infmt='p',inco='f',outfmt='l',outco='g')
 pix_glat = ll(*,1)
 pix_glon = ll(*,0)
 nq = WHERE(pix_glon GT 180.,cq)
 IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) - 360.
 nq = WHERE(pix_glon LT -180.,cq)
 IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) + 360.
 bad = WHERE((ABS(pix_glon) LE loncut)and(ABS(pix_glat) LE latcut),cbad)
 IF (cbad GT 0) THEN vmask(bad) = 0
 ;
 bad = WHERE(sky_glitch GT glmax,cbad)
 IF (cbad GT 0) THEN vmask(bad) = 0
 ;
 bad = WHERE((nifgs LT nr(0))or(nifgs GT nr(1)),cbad)
 IF (cbad GT 0) THEN vmask(bad) = 0
 ;

 ; Initialize FSL Weights
 ; ----------------------
 fsl_sky_wgts = sky_wgts_ds
 fsl_cal_wgts = cal_wgts_ds
 ;

 ; Restore Coadd Variances and D-Vector
 ; ------------------------------------
 PRINT,' '
 PRINT,'Restoring ' + intrans(0) + 'RHS_VAR.ISS'
 RESTORE,'csdr$firas_in:rhs_var.iss'
 ;

 ; Index of "Good" Frequencies
 ; ---------------------------
 fidx = WHERE((f_rhs GT freq_band(0))and(f_rhs LT freq_band(1)),nf)
 IF (nf LT 1) THEN BEGIN
  PRINT,'FMD_WGTS_RHS : No Data in Frequency Range !'
  RETURN
 ENDIF
 ;

; Frequency Averaged Variances
 ; ----------------------------
 ;
 vw = 1. / d_rhs(fidx) / d_rhs(fidx)  ;  Frequency Weight = 1. / DVEC^2.
 vt = ( TRANSPOSE(var_rhs(fidx,*)) # vw ) / TOTAL(vw)           ; SKY
 cal_vt = ( TRANSPOSE(cal_var_rhs(fidx,*)) # vw ) / TOTAL(vw)   ; CAL
 ;

 ; "Good" Coadds
 ; -------------
 good = WHERE((vt GT 0.)and(fsl_sky_wgts GT 0.)and(vmask EQ 1),cgood)
 IF (cgood LE 0) THEN BEGIN
  PRINT,'FMD_WGTS_RHS : No Good Coadds for Variance Model !'
  error = 1
  RETURN
 ENDIF
 ;

 ; Mean Variance per IFG for Good Unmasked Coadds
 ; ----------------------------------------------
 tmg = tm(good)
 var2 = vt(good) * nifgs(good)
 good2 = WHERE((tmg LT tmcut(0))or(tmg GT tmcut(1)),cgood2)
 avar2 = TOTAL(var2(good2)) / cgood2
 ;

 PRINT,' '
 PRINT,'Computing Variance Model Weights'
 PRINT,' '
 slat = STRING(FLOAT(latcut))
 PRINT,'LATCUT =' + STRMID(slat,STRPOS(slat,'.')-3,3)
 slon = STRING(FLOAT(loncut))
 PRINT,'LONCUT =' + STRMID(slon,STRPOS(slon,'.')-3,3)
 PRINT,' '
 ;

 ; Initialize Model Weights
 ; ------------------------
 sky_wgts_mod = fsl_sky_wgts
 cal_wgts_mod = fsl_cal_wgts
 ;

 ; Good Coadd Index
 ; ----------------
 goodv = good
 ;

 ; Good Normalized Variances
 ; -------------------------
 var2n = var2 / avar2
 ;

 sgood = STRCOMPRESS(STRING(cgood))
 ;
 PRINT,'Model Will Be Computed from' + sgood + ' Coadds'
 PRINT,' '
 ;

 ; Compute Variance Model Weights
 ; ------------------------------
 var_coeff_rhs = $
  FMD_MODEL_WGT(var2n=var2n,good=goodv,nifgs=nifgs,tm=tm,sky_glitch=sky_glitch,$
                sky_s0=sky_s0,scan=scan,cal_nifgs=cal_nifgs,cal_tm=cal_tm,$
                cal_glitch=cal_glitch,cal_s0=cal_s0,lp_order=lp_order,$
                tim_up=tm_up,tim_dn=tm_dn,nscan=nscan,gpow=gpow,bpow=bpow,$
                tcut=tmcut,func=sky_var_func,cal_func=cal_var_func,$
                condition=condition,inv_stat=inv_stat,delt=delt,dtmod=dtmod,$
                skut=skut,dscmod=dscmod,vfunc_descrip=vfunc_descrip,$
                sky_wgts=sky_wgts_mod,cal_wgts=cal_wgts_mod)
 ;

 chanscan = 'RHS'
 lp_orderv = lp_order
 sname = 'csdr$firas_out:rhs_var_model.iss'
 SAVE,filename=sname,chanscan,freq_band,lp_orderv,tm_up,tm_dn,gpow,bpow,$
                     nscan,glmax,glcut,nr,tmcut,latcut,loncut,goodv,$
                     var2n,sky_var_func,cal_var_func,vfunc_descrip,$
                     condition,inv_stat,var_coeff_rhs,delt,dtmod,skut,dscmod,$
                     sky_wgts_mod,cal_wgts_mod
                      
 PRINT,' '
 PRINT,'IDL Save Set "' + outtrans(0) + 'RHS_VAR_MODEL.ISS" Created.'
 PRINT,' '
 ;

 ; Apply Model Weights
 ; -------------------
 good = WHERE(fsl_sky_wgts GT 0.,cgood)
 IF (cgood GT 0) THEN fsl_sky_wgts(good) = sky_wgts_mod(good)
 ;
 good = WHERE(fsl_cal_wgts GT 0.,cgood)
 IF (cgood GT 0) THEN fsl_cal_wgts(good) = cal_wgts_mod(good)
 ;

 ; Check for Negative Weights
 ; --------------------------
 nq = WHERE(fsl_sky_wgts LT 0.,cq)
 IF (cq GT 0) THEN fsl_sky_wgts(nq) = 0.
 ;
 nq = WHERE(fsl_cal_wgts LT 0.,cq)
 IF (cq GT 0) THEN fsl_cal_wgts(nq) = 0.
 ;

 ; Re-Weighted Coadd Destriper Weights
 ; -----------------------------------
 ;
 fac = TOTAL(sky_wgts_ds) / TOTAL(fsl_sky_wgts)
 ;
 sky_wgts_ds = fsl_sky_wgts * fac
 cal_wgts_ds = fsl_cal_wgts * fac
 ;

ENDIF
;

;  Pixel Weights and Coadd Fractional Weights
;  ------------------------------------------
error = FMD_PIXEL_WGT(pixel=px,weights=sky_wgts_ds, $
                      frac_wgt=frac_wgt,px_wgt=pixel_wgt)
;

IF (error NE 0) THEN BEGIN
 PRINT,''
 PRINT,'FMD_PIXEL_WGT Returned With Error !'
 PRINT,''
 RETURN
ENDIF
;

; Make RHS_WEIGHTS.ISS Save Set
; --------------------------------
chanscan = 'RHS'
sname = 'csdr$firas_ref:rhs_weights.iss'
SAVE,filename=sname,chanscan,freq_band,px,sky_wgts,sky_wgts_ds,pixel_wgt,$
                    frac_wgt,cal_wgts,cal_wgts_ds,model_var,tm_up,tm_dn,$
                    gpow,bpow,nscan,delt,latcut,loncut,tmcut,lp_orderv
PRINT,' '
PRINT,'IDL Save Set "' + reftrans(0) + 'RHS_WEIGHTS.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
hot_cal=0 & gain_convg=0 & gain_iter=0 & sky_glitch=0 & tm=0 & glat=0 & glon=0
sky_s0=0 & sky_bol_volt=0 & cal_glitch=0 & cal_bol_volt=0 & chanscan_id=0
cal_s0=0 & del_temp=0 & freq_corr=0 & dirbe_array=0 & dirbe_cut=0 & maxfrac=0
xtm=0 & cmn_tm=0 & step_up=0 & step_dn=0 & chan_label=0 & f_hi3=0
ref_s0=0 & rms_s0=0 & cmn_bol=0 & dihed_pow=0 & ref_dihd=0 & n_stripes=0
dihed_cut=0 & dihed_cut_min=0 & cmn_dihd=0 & cvec_cut=0 & latcut=0 & loncut=0
;
zodi_array=0 & max_frac=0 & solution=0 & f=0 & galcut=0 & cal_lbl=0
sky_lbl=0 & scan=0 & time=0 & st_sub=0 & fsl_idx=0
;
badcoadd_lhs=0 & badcoadd_rhf=0 & badcoadd_lhf=0
badcal_lhs=0 & badcal_rhf=0 & badcal_lhf=0
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_lsf=0 & badcoadd_rsf=0
badcal_lls=0 & badcal_rls=0 & badcal_lsf=0 & badcal_rsf=0
badcoadd_llf=0 & badcoadd_rlf=0 & badcal_llf=0 & badcal_rlf=0
d_rhs=0 & n_rhs=0 & b_rhs=0 & l_rhs=0 & f_rhs=0 & galcut_rhs=0
var_tm=0 & cal_var_tm=0 & var_lbl=0 & cal_var_lbl=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
