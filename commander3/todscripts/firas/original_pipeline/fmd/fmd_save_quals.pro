FUNCTION FMD_Save_Quals,bselect,fb=fb,dirb=dirb,db_cut=db_cut,zod=zod,$
         lp=lp,x_tm=x_tm,horn=horn,pow=pow,cut=cut,refd=refd,refb=refb,rmsb=rmsb,$
         sky_dihd=sky_dihd,tminx=tminx,tmaxx=tmaxx,cal_dihd=cal_dihd,$
         hot=hot,cv_cut=cv_cut,lat=lat,lon=lon,convg=convg,$
         iter=iter,cmn=cmn,nstripes=nstripes,maxfrac=maxfrac
;

;
;  FMD_SAVE_QUALS restores default qualifiers from a reference data set,
;                 modifies selected qualifiers, and creates a new IDL
;                 save set of FMD qualifiers.
;
;  Written By  :  Ken Jensen,  Hughes STX,  13-Jun-97
;
;

error = 1
;

; Logical Translations
; --------------------
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
reftrans = STRUPCASE(reftrans)
;
PRINT,' '
PRINT,'Logical Translations are :'
PRINT,' '
PRINT,'CSDR$FIRAS_REF    == ' + reftrans
PRINT,' '
;

; Restore Reference File of Default Qualifiers
; --------------------------------------------
RESTORE,'csdr$firas_ref:fmd_quals_default.iss'
;

; Change Qualifiers
; -----------------
;
freq_band = fb
dirbe_array = dirb
zodi_array = zod
dirbe_cut = db_cut
lp_order = lp
xtm = x_tm
;
step_up = ['89325']  &  step_dn = ['90270']
IF (horn(1) eq 1) THEN BEGIN
 step_up = ['901931850',step_up]  &  step_dn = ['902081120',step_dn]
ENDIF
IF (horn(0) eq 1) THEN BEGIN
 step_up = ['901391535',step_up]  &  step_dn = ['901931850',step_dn]
ENDIF
;
dihed_pow = pow
dihed_cut = cut
ref_dihd = refd
ref_s0 = refb
rms_s0 = rmsb
good_sky_dihd = sky_dihd
tmin = tminx
tmax = tmaxx
good_cal_dihd = cal_dihd
hot_cal = hot
cvec_cut = cv_cut
latcut = lat
loncut = lon
gain_convg = convg
gain_iter = iter
cmn_tm = cmn(0)
cmn_dihd = cmn(1)
cmn_bol = cmn(2)
n_stripes = nstripes
max_frac = maxfrac
;


; Make New QUALS Save Set
; -----------------------
;

IF (bselect eq 0) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_lo.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
;
 PRINT,' '
 PRINT,'IDL Save Set "' + reftrans(0) + 'FMD_QUALS_LO.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (bselect eq 1) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_hi_1.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
 PRINT,' '
 PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_QUALS_HI_1.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (bselect eq 2) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_hi_2.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
 PRINT,' '
 PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_QUALS_HI_2.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (bselect eq 3) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_hi_3.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
 PRINT,' '
 PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_QUALS_HI_3.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (bselect eq 4) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_hi_4.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
 PRINT,' '
 PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_QUALS_HI_4.ISS" Created.'
 PRINT,' '
;
ENDIF
;

IF (bselect eq 5) THEN BEGIN
;
 sname = 'csdr$firas_ref:fmd_quals_hr.iss'
 SAVE,filename=sname,del_temp,freq_corr,freq_band,dirbe_array,dirbe_cut, $
               lp_order,xtm,cmn_tm,step_up,step_dn,ref_s0,rms_s0,cmn_bol, $
               dihed_pow,ref_dihd,dihed_cut,dihed_cut_min,cmn_dihd, $
               zodi_array,good_sky_dihd,tmin,tmax,good_cal_dihd,hot_cal, $
               cvec_cut,latcut,loncut,n_stripes,max_frac,gain_convg,gain_iter
;
 PRINT,' '
 PRINT,'IDL Save Set "'+ reftrans(0) + 'FMD_QUALS_HR.ISS" Created.'
 PRINT,' '
;
ENDIF
;

badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_lsf=0 & badcoadd_rsf=0
badcal_lls=0 & badcal_rls=0 & badcal_lsf=0 & badcal_rsf=0
badcoadd_lhs=0 & badcoadd_rhs=0 & badcoadd_lhf=0 & badcoadd_rhf=0
badcal_lhs=0 & badcal_rhs=0 & badcal_lhf=0 & badcal_rhf=0
badcoadd_llf=0 & badcoadd_rlf=0 & badcal_llf=0 & badcal_rlf=0
;
error = 0
;

RETURN,error
END
