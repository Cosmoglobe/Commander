Pro FMD_DIRBE_FUNCR
;

;
;  FMD_DIRBE_FUNCR Makes the FMD_DIRBE_FUNC_HRES Reference Data Set
;
;
;  Required Logicals :
;
;     CSDR$FIRAS_IN  : Directory containing HRES.ISS and HRES_WEIGHTS.ISS.
;
;     CSDR$FIRAS_REF : Directory containing FMD_DIRBE_FUNC0_LHRES.ISS, and
;                      where FMD_DIRBE_FUNC_HRES.ISS will be sent.
;
;
;  Written by : Ken Jensen, Hughes STX, 18-Jun-1997
;
;

RESTORE,'csdr$firas_in:hres.iss'
;
RESTORE,'csdr$firas_in:hres_weights.iss'
;
RESTORE,'csdr$firas_ref:fmd_dirbe_func0_hres.iss'
;

g8x = g8
g9x = g9
g10x = g10
;
gx = 0. * g8
ng = WHERE(sky_wgts_ds GT 0.,cg)
pxg = px(ng)
wg = sky_wgts_ds(ng)
apx = FLTARR(cg)
h = 0. * apx
;

h(*) = g8(ng(*)) + f8(pxg)
;
dummy = FMD_PIX_SUM(pixel=pxg,data=h*wg,cmp_px=cmp_px)
dummy2 = FMD_PIX_SUM(pixel=pxg,data=wg,cmp_px=cmp_px)
ap = dummy / dummy2
;
FOR i = LONG(0),LONG(cg)-1 DO apx(i) = ap(WHERE(cmp_px eq pxg(i)))
;
gx(ng) = h - apx
g8 = gx
;

h(*) = g9(ng(*)) + f9(pxg)
;
dummy = FMD_PIX_SUM(pixel=pxg,data=h*wg,cmp_px=cmp_px)
dummy2 = FMD_PIX_SUM(pixel=pxg,data=wg,cmp_px=cmp_px)
ap = dummy / dummy2
;
FOR i = LONG(0),LONG(cg)-1 DO apx(i) = ap(WHERE(cmp_px eq pxg(i)))
;
gx(ng) = h - apx
g9 = gx
;

h(*) = g10(ng(*)) + f10(pxg)
;
dummy = FMD_PIX_SUM(pixel=pxg,data=h*wg,cmp_px=cmp_px)
dummy2 = FMD_PIX_SUM(pixel=pxg,data=wg,cmp_px=cmp_px)
ap = dummy / dummy2
;
FOR i = LONG(0),LONG(cg)-1 DO apx(i) = ap(WHERE(cmp_px eq pxg(i)))
;
gx(ng) = h - apx
g10 = gx
;

idd = 'HRES'
sname = 'csdr$firas_ref:fmd_dirbe_func_hres.iss'
SAVE,filename=sname,idd,g8x,g9x,g10x,g8,g9,g10,f8,f9,f10,sky_idx
PRINT,' '
PRINT,'IDL Save Set "' + STRUPCASE(sname) + '" Created.'
PRINT,' '
;

chanscan=0 & sky_wgts=0 & pixel_wgt=0 & frac_wgt=0 & cal_wgts=0 & cal_wgts_ds=0
chan_label=0 & sky_idx=0 & tm=0 & glat=0 & glon=0 & sky_dihd=0 & sky_glitch=0
sky_s0=0 & scan=0 & sky_wgts=0 & nifgs=0 & cal_idx=0 & xcal=0 & cal_tm=0
ical=0 & refh=0 & skyh=0 & dihd=0 & cal_glitch=0 & cal_s0=0 & cal_wgts=0
cal_nifgs=0 & freq_hres=0 & scale_factor=0 & solution=0 & st_sub=0
fsl_idx=0 & sky_lbl=0 & cal_lbl=0
;

RETURN
END
