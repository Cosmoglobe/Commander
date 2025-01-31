Pro FMD_Bad_LHS,badcoadd_lhs,badcal_lhs,error
;

if N_Params() ne 3 then begin
 print,'FMD_Bad_LHS,badcoadd_lhs,badcal_lhs,error'
 error = 1
 return
endif
;

error = 1
;

RESTORE,'csdr$firas_ref:fmd_quals_default.iss'
badcoadd_old = badcoadd_lhs
badcal_old = badcal_lhs
;

RESTORE,'csdr$firas_out:lhs_2x_sky_chi2.iss'
wd0 = TEMPORARY(sky_wgts_ds)
IF (MAX(badcoadd_old) GE 0) THEN BEGIN
 dummy = MAX(ABS(wd0(badcoadd_old)))
 IF (dummy NE 0.) THEN BEGIN
  PRINT,'FMD_BAD_LHS : Error in SKY_CHI2 Weights !'
  RETURN
 ENDIF
ENDIF
;

RESTORE,'csdr$firas_out:lhs_2x_cal_chi2.iss'
cwd0 = TEMPORARY(cal_wgts_ds)
IF (MAX(badcal_old) GE 0) THEN BEGIN
 dummy = MAX(ABS(cwd0(badcal_old)))
 IF (dummy NE 0.) THEN BEGIN
  PRINT,'FMD_BAD_LHS : Error in CAL_CHI2 Weights !'
  RETURN
 ENDIF
ENDIF
;

ct2 = TOTAL(sky_chi2_lhs_2,1)
cct2 = TOTAL(cal_chi2_lhs_2,1)
;
good = WHERE((wd0 GT 0.)and(frac_wgt LT 1.)and(sky_mask EQ 1),cg)
goodc = WHERE(cwd0 GT 0.,ccg)
;

RESTORE,'csdr$firas_out:lhs_3x_sky_chi2.iss'
dummy = MAX(ABS(wd0-sky_wgts_ds))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD_LHS : Error in SKY_CHI2 Weights !'
 RETURN
ENDIF
;
RESTORE,'csdr$firas_out:lhs_3x_cal_chi2.iss'
dummy = MAX(ABS(cwd0-cal_wgts_ds))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD_LHS : Error in CAL_CHI2 Weights !'
 RETURN
ENDIF
;
ct3 = TOTAL(sky_chi2_lhs_3,1)
cct3 = TOTAL(cal_chi2_lhs_3,1)
;

RESTORE,'csdr$firas_out:lhs_4x_sky_chi2.iss'
dummy = MAX(ABS(wd0-sky_wgts_ds))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD_LHS : Error in SKY_CHI2 Weights !'
 RETURN
ENDIF
;
RESTORE,'csdr$firas_out:lhs_4x_cal_chi2.iss'
dummy = MAX(ABS(cwd0-cal_wgts_ds))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD_LHS : Error in CAL_CHI2 Weights !'
 RETURN
ENDIF
;
ct4 = TOTAL(sky_chi2_lhs_4,1)
cct4 = TOTAL(cal_chi2_lhs_4,1)
;

ctg = (ct2(good) + ct3(good) + ct4(good)) / 170.
ctg2 = ct2(good) / 55.
ctg3 = ct3(good) / 55.
ctg4 = ct4(good) / 60.
;
ctmax = 1. + 1.5 * SQRT(43./170)
ctmax2 = 1. + 1.5 * SQRT(43./55)
ctmax3 = 1. + 1.5 * SQRT(43./55)
ctmax4 = 1. + 1.5 * SQRT(43./60)
;
bad = $
 WHERE((ctg GT ctmax)or(ctg2 GT ctmax2)or(ctg3 GT ctmax3)or(ctg4 GT ctmax4),cbad)
;

IF (cbad GT 0) THEN BEGIN
 bad_new = [badcoadd_old,good(bad)]
 badcoadd_lhs = bad_new(SORT(bad_new))
 nq = WHERE(badcoadd_lhs GE 0,cq)
 IF (cq LE 0) THEN badcoadd_lhs = [-1]
 IF (cq GT 0) THEN badcoadd_lhs = badcoadd_lhs(nq)
ENDIF
;

ctg = (cct2(goodc) + cct3(goodc) + cct4(goodc)) / 170.
ctg2 = cct2(goodc) / 55.
ctg3 = cct3(goodc) / 55.
ctg4 = cct4(goodc) / 60.
;

bad = $
 WHERE((ctg GT ctmax)or(ctg2 GT ctmax2)or(ctg3 GT ctmax3)or(ctg4 GT ctmax4),cbad)
;

IF (cbad GT 0) THEN BEGIN
 bad_new = [badcal_old,goodc(bad)]
 badcal_lhs = bad_new(SORT(bad_new))
 nq = WHERE(badcal_lhs GE 0,cq)
 IF (cq LE 0) THEN badcal_lhs = [-1]
 IF (cq GT 0) THEN badcal_lhs = badcal_lhs(nq)
ENDIF
;

chanscan=0 & px=0 & sky_wgts=0 & pixel_wgt=0 & cal_wgts=0 & chan_label=0
stripes_descrip=0 & freq_dependence=0 & tm=0 & sky_s0=0 & sky_dihd=0
sky_glitch=0 & cvec_mask=0 & xcal=0 & cal_tm=0 & ical=0 & skyh=0 & refh=0
dihd=0 & cal_glitch=0 & cal_s0=0
del_temp=0 & freq_corr=0 & dihed_cut_min=0
freq_lhs_2=0 & cvec_lhs_2=0
freq_lhs_3=0 & cvec_lhs_3=0
freq_lhs_4=0 & cvec_lhs_4=0
badcoadd_lhf=0 & badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_rlf=0
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lsf=0
badcoadd_llf=0 & badcal_lhf=0 & badcal_rhs=0 & badcal_rhf=0 & badcal_rlf=0 
badcal_llf=0 & badcal_lls=0 & badcal_rls=0 & badcal_rsf=0 & badcal_lsf=0
;

error = 0
;

RETURN
END
