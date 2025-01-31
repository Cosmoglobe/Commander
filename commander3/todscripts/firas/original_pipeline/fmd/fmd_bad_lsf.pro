Pro FMD_Bad_LSF,badcoadd_lsf,badcal_lsf,error
;

if N_Params() ne 3 then begin
 print,'FMD_Bad_LSF,badcoadd_lsf,badcal_lsf,error'
 error = 1
 return
endif
;

error = 1
;

RESTORE,'csdr$firas_ref:fmd_quals_default.iss'
badcoadd_old = badcoadd_lsf
badcal_old = badcal_lsf
;

RESTORE,'csdr$firas_out:lsfx_sky_chi2.iss'
wd0 = TEMPORARY(sky_wgts_ds)
bad = WHERE(wd0 LE 0.,cbad)
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcoadd_old,bad]
 badcoadd_old = bad_new(SORT(bad_new))
 nq = WHERE(badcoadd_old GE 0,cq)
 IF (cq LE 0) THEN badcoadd_old = [-1]
 IF (cq GT 0) THEN badcoadd_old = badcoadd_old(nq)
ENDIF
;

RESTORE,'csdr$firas_out:lsfx_cal_chi2.iss'
cwd0 = TEMPORARY(cal_wgts_ds)
IF (MAX(badcal_old) GE 0) THEN BEGIN
 dummy = MAX(ABS(cwd0(badcal_old)))
 IF (dummy NE 0.) THEN BEGIN
  PRINT,'FMD_BAD_LSF : Error in CAL_CHI2 Weights !'
  RETURN
 ENDIF
ENDIF
;

ct = TOTAL(sky_chi2_lsf,1)
cct = TOTAL(cal_chi2_lsf,1)
;
good = WHERE((wd0 GT 0.)and(frac_wgt LT 1.)and(sky_mask EQ 1),cg)
goodc = WHERE(cwd0 GT 0.,ccg)
;

ctg = ct(good) / 43.
;
ctmax = 2.5
;
bad = WHERE (ctg GT ctmax,cbad)
;
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcoadd_old,good(bad)]
 badcoadd_lsf = bad_new(SORT(bad_new))
 nq = WHERE(badcoadd_lsf GE 0,cq)
 IF (cq LE 0) THEN badcoadd_lsf = [-1]
 IF (cq GT 0) THEN badcoadd_lsf = badcoadd_lsf(nq)
ENDIF
;

ctg = cct(goodc) / 43.
;
bad = WHERE (ctg GT ctmax,cbad)
;
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcal_old,goodc(bad)]
 badcal_lsf = bad_new(SORT(bad_new))
 nq = WHERE(badcal_lsf GE 0,cq)
 IF (cq LE 0) THEN badcal_lsf = [-1]
 IF (cq GT 0) THEN badcal_lsf = badcal_lsf(nq)
ENDIF
;

chanscan=0 & px=0 & sky_wgts=0 & pixel_wgt=0 & cal_wgts=0 & chan_label=0
stripes_descrip=0 & freq_dependence=0 & tm=0 & sky_s0=0 & sky_dihd=0
sky_glitch=0 & cvec_mask=0 & xcal=0 & cal_tm=0 & ical=0 & skyh=0 & refh=0
dihd=0 & cal_glitch=0 & cal_s0=0
del_temp=0 & freq_corr=0 & dihed_cut_min=0
freq_lsf=0 & cvec_lsf=0
badcoadd_rhf=0 & badcoadd_lhs=0 & badcoadd_rhs=0 & badcoadd_rlf=0
badcoadd_lls=0 & badcoadd_rls=0 & badcoadd_rsf=0 & badcoadd_lhf=0
badcoadd_llf=0 & badcal_rhf=0 & badcal_lhs=0 & badcal_rhs=0 & badcal_rlf=0 
badcal_llf=0 & badcal_lls=0 & badcal_rls=0 & badcal_rsf=0 & badcal_lhf=0
;

error = 0
;

RETURN
END
