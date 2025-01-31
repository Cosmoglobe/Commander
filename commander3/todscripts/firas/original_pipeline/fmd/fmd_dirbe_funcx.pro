Pro FMD_DIRBE_FUNCX,idd
;

if N_Params() ne 1 then begin
 print,'FMD_DIRBE_FUNCX,idd'
 return
endif
;

RESTORE,'csdr$firas_ref:' + idd + 'x_weights.iss'
;

iddx = STRMID(idd,0,3)
RESTORE,'csdr$firas_ref:fmd_dirbe_func0_' + iddx + '.iss'
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

sname = 'csdr$firas_ref:fmd_dirbe_funcx_' + idd + '.iss'
SAVE,filename=sname,idd,g8x,g9x,g10x,g8,g9,g10,f8,f9,f10
PRINT,' '
PRINT,'IDL Save Set "' + STRUPCASE(sname) + '" Created.'
PRINT,' '
;

chanscan=0 & sky_wgts=0 & pixel_wgt=0 & frac_wgt=0 & cal_wgts=0 & cal_wgts_ds=0
solution=0 & fsl_idx=0 & sky_lbl=0
;

RETURN
END
