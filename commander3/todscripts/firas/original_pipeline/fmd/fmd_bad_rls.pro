Pro FMD_Bad_RLS,badcoadd_rls,badcal_rls,error
;

if N_Params() ne 3 then begin
 print,'FMD_Bad_RLS,badcoadd_rls,badcal_rls,error'
 error = 1
 return
endif
;

error = 1
;

RESTORE,'csdr$firas_ref:fmd_quals_default.iss'
badcoadd_old = badcoadd_rls
badcal_old = badcal_rls
;

RESTORE,'csdr$firas_out:rlsx_sky_chi2.iss'
wd0 = TEMPORARY(sky_wgts_ds)
dummy = MAX(ABS(badcoadd_old-WHERE(wd0 EQ 0.)))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD_RLS : Error in SKY_CHI2 Weights !'
 RETURN
ENDIF
;

RESTORE,'csdr$firas_out:rlsx_cal_chi2.iss'
cwd0 = TEMPORARY(cal_wgts_ds)
IF (MAX(badcal_old) GE 0) THEN BEGIN
 dummy = MAX(ABS(cwd0(badcal_old)))
 IF (dummy NE 0.) THEN BEGIN
  PRINT,'FMD_BAD_RLS : Error in CAL_CHI2 Weights !'
  RETURN
 ENDIF
ENDIF
;

ct = TOTAL(sky_chi2_rls,1)
cct = TOTAL(cal_chi2_rls,1)
;
good = WHERE((wd0 GT 0.)and(frac_wgt LT 1.)and(sky_mask EQ 1),cg)
goodc = WHERE(cwd0 GT 0.,ccg)
;

ctg = ct(good) / 43.
;
ctmax = 2.5
bad = WHERE (ctg GT ctmax,cbad)
;
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcoadd_old,good(bad)]
 badcoadd_rls = bad_new(SORT(bad_new))
 nq = WHERE(badcoadd_rls GE 0,cq)
 IF (cq LE 0) THEN badcoadd_rls = [-1]
 IF (cq GT 0) THEN badcoadd_rls = badcoadd_rls(nq)
ENDIF
;

ctg = cct(goodc) / 43.
;
bad = WHERE (ctg GT ctmax,cbad)
;
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcal_old,goodc(bad)]
 badcal_rls = bad_new(SORT(bad_new))
 nq = WHERE(badcal_rls GE 0,cq)
 IF (cq LE 0) THEN badcal_rls = [-1]
 IF (cq GT 0) THEN badcal_rls = badcal_rls(nq)
ENDIF
;

chanscan=0 & px=0 & sky_wgts=0 & pixel_wgt=0 & cal_wgts=0 & chan_label=0
stripes_descrip=0 & freq_dependence=0 & tm=0 & sky_s0=0 & sky_dihd=0
sky_glitch=0 & cvec_mask=0 & xcal=0 & cal_tm=0 & ical=0 & skyh=0 & refh=0
dihd=0 & cal_glitch=0 & cal_s0=0
del_temp=0 & freq_corr=0 & dihed_cut_min=0
freq_rls=0 & cvec_rls=0
badcoadd_lhf=0 & badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_rlf=0
badcoadd_lhs=0 & badcoadd_lls=0 & badcoadd_rsf=0 & badcoadd_lsf=0
badcoadd_llf=0 & badcal_lhf=0 & badcal_rhs=0 & badcal_rhf=0 & badcal_rlf=0 
badcal_llf=0 & badcal_lhs=0 & badcal_lls=0 & badcal_rsf=0 & badcal_lsf=0
;

error = 0
;

RETURN
END
