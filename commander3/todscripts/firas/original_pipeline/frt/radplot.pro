Pro Radplot, DEV=pltdev

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Author:  Fred Shuman,  STX,  1991 Sep 23
; Purpose:  Plot the Van Allen Belts and South Atlantic Anomaly on a terrestrial
;           longitude & latitude map.
; Parameters:
;         pltdev:  plot device (= 'PS' or 'REGIS' (Default).)
; Input:
;       CSDR$FIRAS_REF:FEX_VABSAA.DAT
; Output:
;       Plot; if dev='ps' (postscript) is chosen, default filename is RAD.PS
;------------------------------------------------------------------------------

If N_Elements(pltdev) eq 0 Then pltdev = 'REGIS'
pltdev = strupcase(pltdev)  &  set_plot, pltdev
If (pltdev eq 'PS') Then device, filename='RAD.PS', /landscape

dset = "CSDR$FIRAS_REF:FEX_VABSAA.DAT"  &   OpenR, lunr, dset, /get_lun, /share
finfo = FStat(lunr)                     &   trec=finfo.size/finfo.rec_len
r0 = Assoc(lunr, fex_vabsaa_st(trec))   &   rad = r0(0)

; Note: number of points is one more than number of intervals
nvab     = Fix(360./rad.vab(0).lonstep + .5)
vablon   = FIndgen(nvab+1)*rad.vab(0).lonstep - 180.
vabnlatn = rad.vab(0).latn  &  vabnlats = rad.vab(0).lats
vabslatn = rad.vab(1).latn  &  vabslats = rad.vab(1).lats
nsaa     = Fix((rad.saa.lonmax - rad.saa.lonmin)/rad.saa.lonstep + .5)
saalon   = FIndgen(nsaa+1)*rad.saa.lonstep + rad.saa.lonmin
; Extend the SAA arrays by the last SAA-North point so as to close the SAA plot
saalon   = [saalon, saalon(nsaa)]
saalatn  = [rad.saa.latn, rad.saa.latn(nsaa)]
saalats  = [rad.saa.lats, rad.saa.latn(nsaa)]

Plot, vablon, vabnlatn, title='EARTH RADIATION REGIONS from REF FILES', $
    xrange=[-180.,180.], yrange=[-90.,90.], ystyle=1, position=[0,0,1,.82], $
    xstyle=1
OPlot, vablon, vabnlats  &  OPlot, vablon, vabslatn  &  OPlot, vablon, vabslats
OPlot, saalon, saalatn   &  OPlot, saalon, saalats
;
if(pltdev eq 'PS')then device,/close
;
Return
End
