Pro FMD_Bad2_RSF,badcoadd_rsf,badcal_rsf,error
;

if N_Params() ne 3 then begin
 print,'FMD_Bad2_RSF,badcoadd_rsf,badcal_rsf,error'
 error = 1
 return
endif
;

error = 1
;

RESTORE,'csdr$firas_ref:fmd_bad_coadd_lo.iss'
badcoadd_old = badcoadd_rsf
badcal_old = badcal_rsf
;

RESTORE,'csdr$firas_out:lowf_sky_chi2.iss'
nx = WHERE(sky_idx EQ 3,cx)
wd0 = TEMPORARY(sky_wgts_ds(nx))
dummy = MAX(ABS(badcoadd_old-WHERE(wd0 EQ 0.)))
IF (dummy NE 0.) THEN BEGIN
 PRINT,'FMD_BAD2_RSF : Error in SKY_CHI2 Weights !'
 RETURN
ENDIF
;

RESTORE,'csdr$firas_out:lowf_cal_chi2.iss'
nxx = WHERE(cal_idx EQ 3,cxx)
cwd0 = TEMPORARY(cal_wgts_ds(nxx))
IF (MAX(badcal_old) GE 0) THEN BEGIN
 dummy = MAX(ABS(cwd0(badcal_old)))
 IF (dummy NE 0.) THEN BEGIN
  PRINT,'FMD_BAD2_RSF : Error in CAL_CHI2 Weights !'
  RETURN
 ENDIF
ENDIF
;

ct = TOTAL(sky_chi2_lowf(*,nx),1)
cct = TOTAL(cal_chi2_lowf(*,nxx),1)
;
good = WHERE((wd0 GT 0.)and(frac_wgt(nx) LT 1.)and(sky_mask(nx) EQ 1),cg)
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
 badcoadd_rsf = bad_new(SORT(bad_new))
 nq = WHERE(badcoadd_rsf GE 0,cq)
 IF (cq LE 0) THEN badcoadd_rsf = [-1]
 IF (cq GT 0) THEN badcoadd_rsf = badcoadd_rsf(nq)
ENDIF
;

ctg = cct(goodc) / 43.
;
bad = WHERE (ctg GT ctmax,cbad)
;
IF (cbad GT 0) THEN BEGIN
 bad_new = [badcal_old,goodc(bad)]
 badcal_rsf = bad_new(SORT(bad_new))
 nq = WHERE(badcal_rsf GE 0,cq)
 IF (cq LE 0) THEN badcal_rsf = [-1]
 IF (cq GT 0) THEN badcal_rsf = badcal_rsf(nq)
ENDIF
;

chanscan=0 & px=0 & sky_wgts=0 & pixel_wgt=0 & cal_wgts=0 & chan_label=0
stripes_descrip=0 & freq_dependence=0 & tm=0 & sky_s0=0 & sky_dihd=0
sky_glitch=0 & cvec_mask=0 & xcal=0 & cal_tm=0 & ical=0 & skyh=0 & refh=0
dihd=0 & cal_glitch=0 & cal_s0=0
del_temp=0 & freq_corr=0 & dihed_cut_min=0 & cal_wgts_fsl=0
freq_lowf=0 & cvec_lowf=0
badcoadd_lhf=0 & badcoadd_rhs=0 & badcoadd_rhf=0 & badcoadd_rlf=0
badcoadd_lhs=0 & badcoadd_rls=0 & badcoadd_lsf=0 & badcoadd_lls=0
badcoadd_llf=0 & badcal_lhf=0 & badcal_rhs=0 & badcal_rhf=0 & badcal_rlf=0 
badcal_llf=0 & badcal_lhs=0 & badcal_rls=0 & badcal_lsf=0 & badcal_lls=0
;
error = 0
;

RETURN
END
