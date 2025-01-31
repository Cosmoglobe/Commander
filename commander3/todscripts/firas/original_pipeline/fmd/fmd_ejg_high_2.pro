Pro FMD_EJG_HIGH_2,error
;

;
;  FMD_EJG_HIGH_2 drives the FMD_DESTRIPER procedure to create
;  IDL save sets of combined high channel Band_2 offset stripes.
;
;
;  ARGUMENTS (I/O)    :
;
;   ERROR (O)         :  Return Error Status
;
;
;  PROGRAMS Called    :  FMD_PIX_SUM
;                        FMD_FUNC
;                        FMD_DESTRIPER
;                        FMD_ERRMAT
;                        COORCONV
;                        PIX2XY
;
;
;  Required Logicals  :
;
;    CSDR$FIRAS_REF   =  Directory containing FMD_QUALS_HI_2.ISS,
;                        FMD_DIRBE_FUNC_HIGH.ISS, and FMD_BAD_COADD_HI.ISS .
;
;    CSDR$FIRAS_IN    =  Directory containing HIGH.ISS, HIGH_WEIGHTS.ISS,
;                        HIGH_2_CALSPEC_ZSUB.ISS, and HIGH_2_ZODI_SKYMAP.ISS .
;
;    CSDR$FIRAS_OUT   =  Directory where output IDL save sets will be sent.
;
;
;  Output  :
;
;    HIGH_2_EJG.ISS, contains destriper qualifiers, weights, and stripe spectra.
;
;    HIGH_2_ERRORS.ISS, contains destriper error matrices.
;
;    HIGH_2_SKYMAP.ISS, contains skymap spectra and associated errors.
;
;    HIGH_2_FUNC.ISS, contains destriper kernals.
;
;
;    EXAMPLE : $ uidl
;              UIDL> FMD_EJG_HIGH_2,error
;
;
;    HISTORY : Written by Ken Jensen, Hughes STX, 18-Jun-1997.
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
 PRINT,'FMD_EJG_HIGH_2 : Incorrect invocation !'
 PRINT,' '
 PRINT,'Correct invocation is :'
 PRINT,' '
 PRINT,'IDL> FMD_EJG_HIGH_2,error'
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
;

nx = WHERE(dirbe_array EQ 1,n_dirbe)
;

n_tophat = N_ELEMENTS(step_up)
IF (n_tophat NE N_ELEMENTS(step_dn)) THEN BEGIN
 PRINT,'FMD_EJG_HIGH_2 : STEP_UP and STEP_DN are NOT Compatible !'
 RETURN
ENDIF
;

; Time Kernels
; ------------
;
lp_order = FIX(lp_order)   ; Must be Integers
;
dummy = WHERE(lp_order GT 0,n_tm)
IF (n_TM GT 1) THEN BEGIN
;
 hlp = HISTOGRAM(lp_order,min=1)
 IF (MAX(hlp) GT 1) THEN BEGIN
  PRINT,'FMD_EJG_HIGH_2 : Elements of LP_ORDER Are Degenerate !'
  RETURN
 ENDIF
;
 lp_order = lp_order(SORT(lp_order))  ;  FMD_FUNC Requires Sorted LP_ORDER
;
ENDIF
;
n_time = n_tm
;

; Dihedrals
; ---------
n_dihd = 0
IF (TOTAL(dihed_pow) GT 0.) THEN BEGIN
 n_dihd = N_ELEMENTS(dihed_pow)
 IF (dihed_cut ge dihed_cut_min) THEN n_dihd = 2 * n_dihd
ENDIF
;

; Bolometers
; ----------
n_bol = 0
IF (rms_s0 GT 0.) THEN BEGIN
 n_bol = 1
 IF (cmn_bol eq 'Y') THEN n_bol = 2
ENDIF
;

; Bad Coadd Indices
; -----------------
RESTORE,'csdr$firas_ref:fmd_bad_coadd_hi.iss'
;

; List of FMD Qualifiers
; ----------------------
PRINT,' '
PRINT,'FMD Qualifiers are :'
PRINT,' '
PRINT,'FREQ_BAND =',freq_band          ;  Frequency Range (icm)
PRINT,'DEL_TEMP =',del_temp            ;  XCAL adjustment (K)
PRINT,'FREQ_CORR =',freq_corr          ;  Frequency scale correction
PRINT,' '
PRINT,'DIRBE_ARRAY =',dirbe_array      ;  DIRBE kernel selection
PRINT,' '
PRINT,'LP_ORDER =',lp_order            ;  Orders of Legendre kernels
PRINT,' '
PRINT,'XTM =',xtm                      ;  Defines X(t) : Pn(t) = Pn(X(t))
PRINT,' '
PRINT,'CMN_TM = ' + cmn_tm             ;  Common Time Flag
PRINT,' '
IF (TOTAL(dihed_pow) LE 0.) THEN PRINT,' No Dihedral Stripes'
IF (TOTAL(dihed_pow) GT 0.) THEN BEGIN
 npow = N_ELEMENTS(dihed_pow)
 spow1 = STRTRIM(STRCOMPRESS(STRING(npow)),2)
 spow2 = STRTRIM(STRCOMPRESS(STRING(2*npow)),2)
 IF (dihed_cut lt dihed_cut_min) THEN PRINT,spow1 + ' Dihedral Stripe'
 IF (dihed_cut ge dihed_cut_min) THEN PRINT,spow2 + ' Dihedral Stripes'
;
 PRINT,'DIHED_POW =',dihed_pow           ;  Dihed kernels = DIHED ** dihed_pow
 IF (dihed_cut ge dihed_cut_min) THEN $
  PRINT,'DIHED_CUT =',dihed_cut          ;  Dihedral cut-off
 PRINT,'REF_DIHD =',ref_dihd             ;  Reference Dihedral_Temp (K)
 PRINT,'CMN_DIHD = ' + cmn_dihd          ;  Common Dihedral Flag
ENDIF
PRINT,' '
IF (n_bol NE 1) THEN PRINT,'No Bolometer Stripes'
IF (n_bol EQ 1) THEN PRINT,'LHS, RHS, LHF, and RHF Bolometer Stripes'
IF (n_bol EQ 2) THEN PRINT,'LH and RH Bolometer Stripes'
IF (n_bol GT 0) THEN BEGIN
 PRINT,'REF_S0 =',ref_s0                ; Reference Bolometer Responsivity
 PRINT,'RMS_S0 =',rms_s0                ; Approx. RMS of Bolometer Responsivity
 PRINT,'CMN_BOL = ' + cmn_bol           ; Common Bolometer Flag (Left and Right)
ENDIF
PRINT,' '
PRINT,'STEP_UP =',step_up               ; Offset period jstart
PRINT,'STEP_DN =',step_dn               ; Offset period jstop
PRINT,' '
PRINT,'TMIN =',tmin                     ; Minimum controllable temps
PRINT,'TMAX =',tmax                     ; Maximum controllable temps
PRINT,'GOOD_CAL_DIHD =',good_cal_dihd   ; Good Cal DIHED temps
PRINT,'GOOD_SKY_DIHD =',good_sky_dihd   ; Good Sky DIHED temps
PRINT,' '
PRINT,'BADCOADD_LHS =',badcoadd_lhs     ; Bad LHS SKY Coadds
PRINT,' '
PRINT,'BADCOADD_RHS =',badcoadd_rhs     ; Bad RHS SKY Coadds
PRINT,' '
PRINT,'BADCOADD_LHF =',badcoadd_lhf     ; Bad LHF SKY Coadds
PRINT,' '
PRINT,'BADCOADD_RHF =',badcoadd_rhf     ; Bad RHF SKY Coadds
PRINT,' '
PRINT,'BADCAL_LHS =',badcal_lhs         ; Bad LHS CAL Coadds
PRINT,' '
PRINT,'BADCAL_RHS =',badcal_rhs         ; Bad RHS CAL Coadds
PRINT,' '
PRINT,'BADCAL_LHF =',badcal_lhf         ; Bad LHF CAL Coadds
PRINT,' '
PRINT,'BADCAL_RHF =',badcal_rhf         ; Bad RHF CAL Coadds
PRINT,' '
PRINT,'DIRBE_CUT =',dirbe_cut           ; Maximum Allowed DIRBE Gradients
PRINT,' '
PRINT,'LATCUT,LONCUT =',[latcut,loncut] ; Latitude/Longitude Mask for Destriping
PRINT,' '
PRINT,'CVEC_CUT =',cvec_cut             ; Latitude/Longitude Mask for C-Vector
PRINT,' '
PRINT,'N_STRIPES =',n_stripes           ; Number of Stripes
PRINT,' '
PRINT,'MAX_FRAC =',max_frac             ; Maximum Allowed Frac_Wgt for C-Vector
                                        ; and Chi-Squared
PRINT,' '
;

; STRIPE_ORDER (Common Stripes First)
; -----------------------------------
;
stripe_id = ['DIRBE','TIME','DIHEDRAL','BOLOMETER','TOPHAT']
;
stripe_order = INTARR(5)
;
;
IF (n_dirbe GT 0) THEN stripe_order(0) = 1
;
IF ((n_time GT 0)and(cmn_tm EQ 'Y')) THEN $
     stripe_order(1) = MAX(stripe_order) + 1
;
IF ((n_dihd GT 0)and(cmn_dihd EQ 'Y')) THEN $
     stripe_order(2) = MAX(stripe_order) + 1
;
IF ((n_bol GT 0)and(cmn_bol EQ 'Y')) THEN $
     stripe_order(3) = MAX(stripe_order) + 1
;
IF ((n_time GT 0)and(cmn_tm NE 'Y')) THEN $
     stripe_order(1) = MAX(stripe_order) + 1
;
IF ((n_dihd GT 0)and(cmn_dihd NE 'Y')) THEN $
     stripe_order(2) = MAX(stripe_order) + 1
;
IF ((n_bol GT 0)and(cmn_bol NE 'Y')) THEN $
     stripe_order(3) = MAX(stripe_order) + 1
;
IF (n_tophat GT 0) THEN stripe_order(4) = MAX(stripe_order) + 1
;

PRINT,' '
PRINT,'STRIPE_ORDER =',stripe_order
PRINT,' '
;

; Check STRIPE_ORDER
; ------------------
IF (N_ELEMENTS(stripe_order) NE 5) THEN BEGIN
 PRINT,'FMD_EJG_HIGH_2 : STRIPE_ORDER Must Have 5 Elements !'
; STOP
 RETURN
ENDIF
;
stripe_err = 0
;
n1 = WHERE(stripe_order EQ 1,c1)
IF (c1 NE 1) THEN stripe_err = 1
n2 = WHERE(stripe_order EQ 2,c2)
IF (c2 GT 1) THEN stripe_err = 1
n3 = WHERE(stripe_order EQ 3,c3)
IF (c3 GT c2) THEN stripe_err = 1
n4 = WHERE(stripe_order EQ 4,c4)
IF (c4 GT c3) THEN stripe_err = 1
n5 = WHERE(stripe_order EQ 5,c5)
IF (c5 GT c4) THEN stripe_err = 1
;
IF (stripe_err EQ 1) THEN BEGIN
 PRINT,'FMD_EJG_HIGH_2 : Error in STRIPE_ORDER !'
; STOP
 RETURN
ENDIF
;

IF (dihed_cut lt dihed_cut_min) THEN dihed_cut = 0.
;

; Stripes Description String
; --------------------------
stripes_descrip = ''
;
IF (dirbe_array(0) EQ 1) THEN stripes_descrip = stripes_descrip + 'DIRBE_10,'
IF (dirbe_array(1) EQ 1) THEN stripes_descrip = stripes_descrip + 'DIRBE_9,'
IF (dirbe_array(2) EQ 1) THEN stripes_descrip = stripes_descrip + 'DIRBE_8,'
;
IF (n_tm GT 0) THEN BEGIN
 stripes_descrip = stripes_descrip + 'LEGENDRE_'
 FOR i=0,n_tm-1 DO stripes_descrip = $
  stripes_descrip + STRTRIM(STRCOMPRESS(STRING(lp_order(i))),2) + ','
ENDIF
;
IF (MAX(dihed_pow) GT 0.) THEN BEGIN
 IF (cmn_dihd eq 'Y') THEN stripes_descrip = stripes_descrip + 'COMMON_'
 stripes_descrip = stripes_descrip + 'DIHEDRAL'
 IF (dihed_cut LT dihed_cut_min) THEN stripes_descrip = stripes_descrip + ',' 
 IF (dihed_cut GE dihed_cut_min) THEN stripes_descrip = stripes_descrip + '2,' 
ENDIF
;
IF (rms_s0 GT 0.) THEN BEGIN
 IF (cmn_bol eq 'Y') THEN stripes_descrip = stripes_descrip + 'COMMON_'
 stripes_descrip = stripes_descrip + 'BOLOMETER,'
ENDIF
;
IF (N_ELEMENTS(step_up) gt 1) THEN BEGIN
 stripes_descrip = stripes_descrip + '6K,'
 IF (N_ELEMENTS(step_up) gt 2) THEN stripes_descrip = stripes_descrip + '4K,'
ENDIF
stripes_descrip = stripes_descrip + 'MISSION'
;

; Print more Info RE this run
; ---------------------------
PRINT,' '
PRINT, 'HIGH_2 Offset Destriping -- Pass 4'
PRINT,' '
PRINT, 'XCAL temperature adjustment (mK): ' + strtrim(string(del_temp*1000.),2)
PRINT,' '
PRINT, 'Frequency Correction: ' + strtrim(string(freq_corr),2)
PRINT,' '
PRINT, 'Real Part of Spectra will be destriped'
PRINT,' '
;

; Restore Coadd Data
; ------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH.ISS'
RESTORE,'csdr$firas_in:high.iss'
;

; Restore Undestriped ZODI-Subtracted Spectra
; -------------------------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH_2_CALSPEC_ZSUB.ISS'
RESTORE,'csdr$firas_in:high_2_calspec_zsub.iss'
;

; Frequency Band
; --------------
nf = WHERE((f_hi2 GE freq_band(0))and(f_hi2 LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(f_hi2)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_EJG_HIGH_2 : Error in Frequency Array !'
 PRINT,''
 RETURN
ENDIF
;

; Restore ZODI Pixel Spectra
; --------------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH_2_ZODI_SKYMAP.ISS'
RESTORE,'csdr$firas_in:high_2_zodi_skymap.iss'
nf = WHERE((fz_hi2 GE freq_band(0))and(fz_hi2 LE freq_band(1)),cf)
IF (cf NE N_ELEMENTS(fz_hi2)) THEN BEGIN
 PRINT,''
 PRINT,'FMD_EJG_HIGH_2 : Error in ZODI_SKYMAP Frequency Array !'
 RETURN
ENDIF
;

; Restore Coadd Weights
; ---------------------
PRINT,' '
PRINT,'Restoring ' + intrans(0) + 'HIGH_WEIGHTS.ISS'
RESTORE,'csdr$firas_in:high_weights.iss'
;

; Restore DIRBE Gradient Functions
; --------------------------------
PRINT,' '
PRINT,'Restoring ' + reftrans(0) + 'FMD_DIRBE_FUNC_HIGH.ISS'
RESTORE,'csdr$firas_ref:fmd_dirbe_func_high.iss'
;

n_sky = N_ELEMENTS(px)     ; Number of SKY Coadds
n_cal = N_ELEMENTS(xcal)   ; Number of CAL Coadds
;

; Initialize DIRBE Mask
; ---------------------
dirbe_mask = BYTARR(n_sky) + 1B
;

; Apply DIRBE_CUT to Coadds
; -------------------------
IF (dirbe_array(0) eq 1) THEN BEGIN   ; Band 10
 nw = WHERE(ABS(g10) gt dirbe_cut(0),cw)
 IF (cw gt 0) THEN dirbe_mask(nw) = 0       ;  Large Band 10 gradients masked
ENDIF
;
IF (dirbe_array(1) eq 1) THEN BEGIN   ; Band 9
 nw = WHERE(ABS(g9) gt dirbe_cut(1),cw)
 IF (cw gt 0) THEN dirbe_mask(nw) = 0        ;  Large Band 9 gradients masked
ENDIF
;
IF (dirbe_array(2) eq 1) THEN BEGIN   ; Band 8
 nw = WHERE(ABS(g8) gt dirbe_cut(1),cw)
 IF (cw gt 0) THEN dirbe_mask(nw) = 0        ;  Large Band 8 gradients masked
ENDIF
;

nmask = WHERE(dirbe_mask eq 0,cmask)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cmask)) + ' Coadds Masked by DIRBE_CUT'
PRINT,' '
;

; Sky Mask (DIRBE+LATCUT/LONCUT)
; ------------------------------
;
sky_mask = dirbe_mask
;
ll = COORCONV(px,infmt='p',inco='f',outfmt='l',outco='g')
pix_glat = ll(*,1)
pix_glon = ll(*,0)
nq = WHERE(pix_glon gt 180.,cq)
IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) - 360.
nq = WHERE(pix_glon lt -180.,cq)
IF (cq gt 0) THEN pix_glon(nq) = pix_glon(nq) + 360.
bad_gal = WHERE((ABS(pix_glon) LE loncut)and(ABS(pix_glat) LE latcut),cbad)
IF (cbad GT 0) THEN sky_mask(bad_gal) = 0
PRINT,' '
PRINT,STRCOMPRESS(STRING(cbad)) + ' Coadds Masked by LATCUT,LONCUT'
PRINT,' '
;

nmask = WHERE(sky_mask eq 0,cmask)
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
IF (good_cut eq 'Y') THEN BEGIN
 bad_cvec = $
   WHERE((ABS(pix_glat) LE cvec_cut(0))and(ABS(pix_glon) LE cvec_cut(1)),cbad)
 IF (cbad GT 0) THEN cvec_mask(bad_cvec) = 0
;
ENDIF
;
bad_cvec = WHERE(cvec_mask ne 1,cbad)
PRINT,' '
PRINT,STRCOMPRESS(STRING(cbad)) + ' Coadds Masked from C-VECTOR'
PRINT,' '
;

; Coadd Weight Status
; -------------------
PRINT,' '
stot = STRCOMPRESS(STRING(LONG(TOTAL(sky_wgts_ds))))
PRINT,'Total SKY Weight           =' + stot
stot = STRCOMPRESS(STRING(LONG(TOTAL(sky_wgts_ds*sky_mask))))
PRINT,'Total SKY Destriper Weight =' + stot
PRINT,' '
stot = STRCOMPRESS(STRING(LONG(TOTAL(sky_wgts_ds*cvec_mask))))
PRINT,'Total C-Vector Weight =' + stot
PRINT,' '
stot = STRCOMPRESS(STRING(LONG(TOTAL(cal_wgts_ds))))
PRINT,'Total CAL Weight                 =' + stot
stot = STRCOMPRESS(STRING(LONG(TOTAL(cal_wgts_ds))))
PRINT,'Total CAL Destriper Weight       =' + stot
PRINT,' '
;

; Compute n_high, l_high, b_high
; ------------------------------
uv = COORCONV([[glon],[glat]],infmt='l',outfmt='u')

uvavg0 = FMD_PIX_SUM(pix=px,data=uv(*,0),wgt=sky_wgts_ds,cmp_px=cmp_px)
uvavg1 = FMD_PIX_SUM(pix=px,data=uv(*,1),wgt=sky_wgts_ds,cmp_px=cmp_px)
uvavg2 = FMD_PIX_SUM(pix=px,data=uv(*,2),wgt=sky_wgts_ds,cmp_px=cmp_px)

ll = COORCONV([[uvavg0],[uvavg1],[uvavg2]],infmt='u',outfmt='l')

PIX2XY,cmp_px,data=ll(*,0),res=6,/six,ras=l_high
PIX2XY,cmp_px,data=ll(*,1),res=6,/six,ras=b_high

n_tot = FMD_PIX_SUM(pix=px,data=sky_wgts_ds)
PIX2XY,cmp_px,data=n_tot,res=6,/six,ras=n_high


; Build Destriper Kernels
; -----------------------
n_chan = 4
func = FMD_FUNC(g8=g8,g9=g9,g10=g10,dirbe_array=dirbe_array, $
                tm=tm,cal_tm=cal_tm,lp_order=lp_order,xtm=xtm, $
                step_up=step_up,step_dn=step_dn, $
                sky_dihd=sky_dihd,cal_dihd=dihd,ref_dihd=ref_dihd, $
                dihed_pow=dihed_pow,dihed_cut=dihed_cut, $
                sky_s0=sky_s0,cal_s0=cal_s0,ref_s0=ref_s0,rms_s0=rms_s0, $
                sky_idx=sky_idx,cal_idx=cal_idx, $
                n_dirbe=n_dirbe,n_time=n_time,n_tophat=n_tophat, $
                n_dihd=n_dihd,n_bol=n_bol,n_chan=n_chan, $
                stripe_order=stripe_order)
;

; Number of Kernels
; -----------------
n_krnl = n_dirbe + n_time + n_dihd + n_bol + n_tophat
;

; Number of Common Kernels
; ------------------------
n_cmn = n_dirbe
IF (cmn_tm EQ 'Y') THEN n_cmn = n_cmn + n_time
IF (cmn_dihd EQ 'Y') THEN n_cmn = n_cmn + n_dihd
IF (cmn_bol EQ 'Y') THEN n_cmn = n_cmn + n_bol
;

; Number of Stripes
; -----------------
n_stripes = n_cmn + n_chan * (n_krnl - n_cmn)
;

;  Run Destriper
;  -------------
c_high_2 = $
  FMD_DESTRIPER(freq=f_hi2,pixel=px,tm=tm,cal_tm=cal_tm, $
                spec=spec,cal_spec=cal_spec, $
                sky_wgts=sky_wgts_ds,cal_wgts=cal_wgts_ds, $
                frac_wgt=frac_wgt,max_frac=max_frac, $
 	        sky_mask=sky_mask,cvec_mask=cvec_mask, $
                sky_idx=sky_idx,cal_idx=cal_idx, $
                xcal=xcal,del_temp=del_temp, $
	        ejg=ejg_high_2,afp=s_high_2,chi2_chan=chi2_hi2ch, $
                func=func,sky_func=sky_func_high_2,cal_func=cal_func_high_2, $
                pcvr=pcvr,rect=rect,diag=diag,chi2_map=chi2_high_2, $
                n_krnl=n_krnl,n_cmn=n_cmn,square=square,d_inv=d_inv, $
                num_coadd=nc_high,ndf=ndf_ejg_high_2,/pure,/mjy) 
;

IF (MAX(c_high_2) LT 0.) THEN BEGIN
 PRINT,''
 PRINT,'FMD_DESTRIPER Returned With Error !'
 PRINT,''
 error = 1
 RETURN
ENDIF
;

; Restore ZODI to skymap
; ----------------------
s_high_2 = s_high_2 + z_high_2
;

; List of Pixels with Good Data
; -----------------------------
good = WHERE(sky_wgts_ds GT 0.)
out = FMD_PIX_SUM(pixel=px(good),data=0*px(good),cmp_px=cmp_px)
;

; Compute the auxiliary weights and covariance matrices
; -----------------------------------------------------
PRINT,'Computing the auxiliary matrices.'
PRINT,' '
FMD_ERRMAT,pcvr=pcvr,rect=rect,diag=diag,beta=beta,omega=omega,$
           stripe_conv=stripe_conv,ejg=ejg_high_2,gamma=gamma_high_2,lmat=lmat
;

; Destriper Pixel Weight and Stripe_Contribution
; ----------------------------------------------
wmask = sky_wgts_ds(good) * sky_mask(good)
destriper_wgt = FMD_PIX_SUM(pixel=px(good),data=wmask)
;
stripe_contrib = 0 * destriper_wgt
stripe_contrib(WHERE(destriper_wgt GT 0.)) = 1
;

; Make HIGH_2_EJG Save Set
; ------------------------
chan_label = 'LHS_RHS_LHF_RHF'
sname='csdr$firas_out:high_2_ejg.iss'
SAVE,filename=sname,chan_label,stripes_descrip,stripe_order,stripe_id, $
 dirbe_array,lp_order,xtm,cmn_tm,dihed_pow,ref_dihd,dihed_cut,cmn_dihd, $
 ref_s0,rms_s0,cmn_bol,step_up,step_dn,n_krnl,n_cmn,dirbe_cut,latcut,loncut, $
 del_temp,freq_corr,good_sky_dihd,sky_wgts_ds,sky_mask,cvec_mask,cvec_cut, $
 dirbe_mask,tmin,tmax,good_cal_dihd,cal_wgts_ds,f_hi2,ejg_high_2,gamma_high_2, $
 c_high_2,ndf_ejg_high_2,chi2_hi2ch,dihed_cut_min,max_frac, $
 n_stripes,freq_band
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_2_EJG.ISS" Created.'
;

; Make HIGH_2_ERRORS.ISS Save Set
; -------------------------------
sname = 'csdr$firas_out:high_2_errors.iss'
SAVE,filename=sname,chanscan,stripes_descrip,cmp_px,pcvr,rect,diag,beta,omega, $
                    square,d_inv,stripe_conv,destriper_wgt,stripe_contrib, $
                    lmat,f_hi2,c_high_2
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_2_ERRORS.ISS" Created.'
;

; Make HIGH_2_SKYMAP.ISS Save Set
; -----------------------------
sname = 'csdr$firas_out:high_2_skymap.iss'
freq_high_2 = f_hi2
SAVE,filename=sname,chanscan,stripes_descrip,freq_high_2,l_high,b_high,n_high, $
                    nc_high,c_high_2,s_high_2,z_high_2,chi2_high_2, $
                    cmp_px,pcvr,rect,beta,stripe_contrib
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_2_SKYMAP.ISS" Created.'
PRINT,' '
;

; Make HIGH_2_FUNC.ISS Save Set
; -----------------------------
sname = 'csdr$firas_out:high_2_func.iss'
SAVE,filename=sname,chan_label,stripes_descrip,sky_func_high_2,$
                    cal_func_high_2,sky_idx,cal_idx,chan_label
PRINT,' '
PRINT,'IDL Save Set "'+outtrans(0)+'HIGH_2_FUNC.ISS" Created.'
PRINT,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
idd=0 & hot_cal=0 & gain_convg=0 & gain_iter=0 & sky_glitch=0 & cal_glitch=0
frac_wgt=0 & pixel_wgt=0 & zodi_array=0 & scan=0 & sky_wgts=0
nifgs=0 & ical=0 & refh=0 & skyh=0 & cal_nifgs=0 & cal_wgts=0
st_sub=0 & fsl_idx=0 & sky_lbl=0 & cal_lbl=0
freq_high=0 & scale_factor=0 & g8x=0 & g9x=0 & g10x=0 & f8=0 & f9=0 & f10=0
;

; Return with No Error
; --------------------
error = 0
;

RETURN
END
