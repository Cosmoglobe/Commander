Pro ave_cal_res, pltdev, field, range
;+
; PRO AVE_CAL_RES, PLTDEV, FIELD, RANGE
;
; Plots the average of any calibration resistors for the three time periods.
; The data file names are hard coded.
;
;    PLTDEV: plot device: 'REGIS', 'PS', or 'TEK'
;
;    FIELD: must be one of these:
;       'CALRES_AVE_A_LO(0)',
;       'CALRES_AVE_A_LO(1)',
;       'CALRES_AVE_A_LO(2)',
;       'CALRES_AVE_A_LO(3)',
;       'CALRES_AVE_A_HI(0)',
;       'CALRES_AVE_A_HI(1)',
;       'CALRES_AVE_A_HI(2)',
;       'CALRES_AVE_A_HI(3)',
;       'CALRES_AVE_B_LO(0)',
;       'CALRES_AVE_B_LO(1)',
;       'CALRES_AVE_B_LO(2)',
;       'CALRES_AVE_B_LO(3)',
;       'CALRES_AVE_B_HI(0)',
;       'CALRES_AVE_B_HI(1)',
;       'CALRES_AVE_B_HI(2)',
;       'CALRES_AVE_B_HI(3)',
;       'CALRES_DEV_A_LO(0)',
;       'CALRES_DEV_A_LO(1)',
;       'CALRES_DEV_A_LO(2)',
;       'CALRES_DEV_A_LO(3)',
;       'CALRES_DEV_A_HI(0)',
;       'CALRES_DEV_A_HI(1)',
;       'CALRES_DEV_A_HI(2)',
;       'CALRES_DEV_A_HI(3)',
;       'CALRES_DEV_B_LO(0)',
;       'CALRES_DEV_B_LO(1)',
;       'CALRES_DEV_B_LO(2)',
;       'CALRES_DEV_B_LO(3)',
;       'CALRES_DEV_B_HI(0)',
;       'CALRES_DEV_B_HI(1)',
;       'CALRES_DEV_B_HI(2)',
;       'CALRES_DEV_B_HI(3)',
;       Optional. If no field name specified, user will be prompted to
;       select from a menu. Menu will allow user to select ALL fields.
;
;    RANGE: time range:  1 - first, 2 - second, 3 - third, 4 - forth
;           5 - all 4 periods, 6 - concatenation of all periods (entire mission).
;       Optional.  If no range specified, entire mission is plotted.
;
; Author:  Harte Wang.  September 1991
; Modified: Larry Rosen, September 1991 - add comments, RANGE option, 3 periods.
; Modified: K.Jensen (STX), October 4, 1991 - IDL Structure Names and Field
;           Names changed. 'TEK' plot option added. Field Name menu added.
; Modified: Harte Wang, Oct 7, 1991 - RANGE option expanded to 4 periods.
;
;-
print,' '
print,'AVE_CAL_RES, PLTDEV ( , FIELD, RANGE )'
;
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='AVE_CAL_RES.PS',/landscape
If (pltdev eq 'TEK') then begin
 openw,lunp,'AVE_CAL_RES.TEK',/get_lun
 device,plot_to=lunp
Endif
;
FS = SIZE (FIELD)
IF (FS(1) EQ 0)THEN BEGIN
 field=''
 print,' '
 print,'You must Enter a Field name from the Following List :'
 print,"       CALRES_AVE_A_LO(0)"
 print,"       CALRES_AVE_A_LO(1)"
 print,"       CALRES_AVE_A_LO(2)"
 print,"       CALRES_AVE_A_LO(3)"
 print,"       CALRES_AVE_A_HI(0)"
 print,"       CALRES_AVE_A_HI(1)"
 print,"       CALRES_AVE_A_HI(2)"
 print,"       CALRES_AVE_A_HI(3)"
 print,"       CALRES_AVE_B_LO(0)"
 print,"       CALRES_AVE_B_LO(1)"
 print,"       CALRES_AVE_B_LO(2)"
 print,"       CALRES_AVE_B_LO(3)"
 print,"       CALRES_AVE_B_HI(0)"
 print,"       CALRES_AVE_B_HI(1)"
 print,"       CALRES_AVE_B_HI(2)"
 print,"       CALRES_AVE_B_HI(3)"
 print,"       CALRES_DEV_A_LO(0)"
 print,"       CALRES_DEV_A_LO(1)"
 print,"       CALRES_DEV_A_LO(2)"
 print,"       CALRES_DEV_A_LO(3)"
 print,"       CALRES_DEV_A_HI(0)"
 print,"       CALRES_DEV_A_HI(1)"
 print,"       CALRES_DEV_A_HI(2)"
 print,"       CALRES_DEV_A_HI(3)"
 print,"       CALRES_DEV_B_LO(0)"
 print,"       CALRES_DEV_B_LO(1)"
 print,"       CALRES_DEV_B_LO(2)"
 print,"       CALRES_DEV_B_LO(3)"
 print,"       CALRES_DEV_B_HI(0)"
 print,"       CALRES_DEV_B_HI(1)"
 print,"       CALRES_DEV_B_HI(2)"
 print,"       CALRES_DEV_B_HI(3)"
 print,"       ALL"
 print,' '
 print,'Enter Desired Field Name Now'
 read,field
 print,' '
 print,'Selected Field is :'
 print,field
 print,' '
Endif
R = SIZE (RANGE)
IF (R(1) EQ 0) THEN  RANGE = 6		; Sets default to concatenation
;                                         of all 4 periods. (Entire mission).
;
; Hard coded data sets for each time range.
;
dataset_fex = strarr(4)
dataset_fex(0) = "CSDR$FIRAS_REF:FEX_AV_CALRS.DAT;-3"
dataset_fex(1) = "CSDR$FIRAS_REF:FEX_AV_CALRS.DAT;-2"
dataset_fex(2) = "CSDR$FIRAS_REF:FEX_AV_CALRS.DAT;-1"
dataset_fex(3) = "CSDR$FIRAS_REF:FEX_AV_CALRS.DAT"
;
; Set period to plot.
;
CASE RANGE of
   1: Begin & i1 = 0 & i2 = 0 & end
   2: Begin & i1 = 1 & i2 = 1 & end
   3: Begin & i1 = 2 & i2 = 2 & end
   4: Begin & i1 = 3 & i2 = 3 & end
   5: Begin & i1 = 0 & i2 = 3 & end
   6: Begin & i1 = 0 & i2 = 3 & end
ENDCASE
If ( RANGE NE 6 ) Then Begin
 For i = i1,i2 Do Begin
   openr, lunr, dataset_fex(i), /get_lun, /share
   finfo = fstat(lunr)
   trec = finfo.size/finfo.rec_len
   r0 = assoc(lunr,fex_av_calrs_st(trec))
   f0_fex = r0(0)
   nf=1
   IF(strupcase(field) eq 'ALL')then nf=32
   time_fex=conv_adt_sec(f0_fex.ct_head.time)
   time_fex=(time_fex-time_fex(0))/86400.
   emax=max(time_fex) & imax=max(time_fex)
   xmax=max(emax,imax)
   xlab = 'Time Range : ' + STRING (f0_fex(0).ct_head.gmt) + ' to ' + $
      STRING (f0_fex(trec-1).data_stop) +' (day)'
   If(nf eq 1)then begin
     code ='data_cal=f0_fex'+'.'+field
     status=execute(code)
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = field
   Endif
   If(nf eq 32)then begin
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(0)
     status=execute(code)
     yfield='calres_ave_a_lo(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(1)
     status=execute(code)
     yfield='calres_ave_a_lo(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(2)
     status=execute(code)
     yfield='calres_ave_a_lo(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(3)
     status=execute(code)
     yfield='calres_ave_a_lo(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(0)
     status=execute(code)
     yfield='calres_ave_a_hi(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(1)
     status=execute(code)
     yfield='calres_ave_a_hi(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(2)
     status=execute(code)
     yfield='calres_ave_a_hi(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(3)
     status=execute(code)
     yfield='calres_ave_a_hi(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(0)
     status=execute(code)
     yfield='calres_ave_b_lo(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(1)
     status=execute(code)
     yfield='calres_ave_b_lo(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(2)
     status=execute(code)
     yfield='calres_ave_b_lo(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(3)
     status=execute(code)
     yfield='calres_ave_b_lo(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(0)
     status=execute(code)
     yfield='calres_ave_b_hi(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(1)
     status=execute(code)
     yfield='calres_ave_b_hi(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(2)
     status=execute(code)
     yfield='calres_ave_b_hi(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(3)
     status=execute(code)
     yfield='calres_ave_b_hi(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(0)
     status=execute(code)
     yfield='calres_dev_a_lo(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(1)
     status=execute(code)
     yfield='calres_dev_a_lo(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(2)
     status=execute(code)
     yfield='calres_dev_a_lo(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(3)
     status=execute(code)
     yfield='calres_dev_a_lo(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(0)
     status=execute(code)
     yfield='calres_dev_a_hi(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(1)
     status=execute(code)
     yfield='calres_dev_a_hi(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(2)
     status=execute(code)
     yfield='calres_dev_a_hi(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(3)
     status=execute(code)
     yfield='calres_dev_a_hi(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(0)
     status=execute(code)
     yfield='calres_dev_b_lo(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(1)
     status=execute(code)
     yfield='calres_dev_b_lo(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(2)
     status=execute(code)
     yfield='calres_dev_b_lo(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(3)
     status=execute(code)
     yfield='calres_dev_b_lo(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(0)
     status=execute(code)
     yfield='calres_dev_b_hi(0)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(1)
     status=execute(code)
     yfield='calres_dev_b_hi(1)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(2)
     status=execute(code)
     yfield='calres_dev_b_hi(2)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(3)
     status=execute(code)
     yfield='calres_dev_b_hi(3)'
     plot, time_fex,data_cal, title= dataset_fex(i), xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
   Endif
 EndFor
Endif
If ( RANGE EQ 6 ) Then Begin
 For i = i1,i2 Do Begin
   openr, lunr, dataset_fex(i), /get_lun, /share
   finfo = fstat(lunr)
   trec = finfo.size/finfo.rec_len
   r0 = assoc(lunr,fex_av_calrs_st(trec))
   if(i gt i1)then begin
     f0_tmp=r0(0)
     f0_fex = [f0_fex,f0_tmp]
     time_fex=[time_fex,conv_adt_sec(f0_tmp.ct_head.time)]
   endif
   if(i eq i1)then begin
     f0_fex=r0(0)
     time_fex=conv_adt_sec(f0_fex.ct_head.time)
   endif
 Endfor
 trec=n_elements(time_fex)
 nf=1
 if(strupcase(field) eq 'ALL')then nf=32
 time_fex=(time_fex-time_fex(0))/86400.
 emax=max(time_fex) & imax=max(time_fex)
 xmax=max(emax,imax)
 xlab = 'Time Range : ' + STRING (f0_fex(0).ct_head.gmt) + ' to ' + $
    STRING (f0_fex(trec-1).data_stop) +' (day)'
 dataset_all = 'CSDR$FIRAS_REF:FEX_AV_CALRS.DAT;*'
 If(nf eq 1)then begin
     code='data_cal=f0_fex'+'.'+field
     status = execute(code)
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = field
 Endif
 If(nf eq 32)then begin
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(0)
     status=execute(code)
     yfield='calres_ave_a_lo(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(1)
     status=execute(code)
     yfield='calres_ave_a_lo(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(2)
     status=execute(code)
     yfield='calres_ave_a_lo(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_lo(3)
     status=execute(code)
     yfield='calres_ave_a_lo(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(0)
     status=execute(code)
     yfield='calres_ave_a_hi(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(1)
     status=execute(code)
     yfield='calres_ave_a_hi(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(2)
     status=execute(code)
     yfield='calres_ave_a_hi(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_a_hi(3)
     status=execute(code)
     yfield='calres_ave_a_hi(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(0)
     status=execute(code)
     yfield='calres_ave_b_lo(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(1)
     status=execute(code)
     yfield='calres_ave_b_lo(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(2)
     status=execute(code)
     yfield='calres_ave_b_lo(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_lo(3)
     status=execute(code)
     yfield='calres_ave_b_lo(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(0)
     status=execute(code)
     yfield='calres_ave_b_hi(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(1)
     status=execute(code)
     yfield='calres_ave_b_hi(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(2)
     status=execute(code)
     yfield='calres_ave_b_hi(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_ave_b_hi(3)
     status=execute(code)
     yfield='calres_ave_b_hi(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(0)
     status=execute(code)
     yfield='calres_dev_a_lo(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(1)
     status=execute(code)
     yfield='calres_dev_a_lo(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(2)
     status=execute(code)
     yfield='calres_dev_a_lo(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_lo(3)
     status=execute(code)
     yfield='calres_dev_a_lo(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(0)
     status=execute(code)
     yfield='calres_dev_a_hi(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(1)
     status=execute(code)
     yfield='calres_dev_a_hi(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(2)
     status=execute(code)
     yfield='calres_dev_a_hi(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_a_hi(3)
     status=execute(code)
     yfield='calres_dev_a_hi(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(0)
     status=execute(code)
     yfield='calres_dev_b_lo(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(1)
     status=execute(code)
     yfield='calres_dev_b_lo(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(2)
     status=execute(code)
     yfield='calres_dev_b_lo(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_lo(3)
     status=execute(code)
     yfield='calres_dev_b_lo(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(0)
     status=execute(code)
     yfield='calres_dev_b_hi(0)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(1)
     status=execute(code)
     yfield='calres_dev_b_hi(1)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(2)
     status=execute(code)
     yfield='calres_dev_b_hi(2)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
     code ='data_cal=f0_fex'+'.calres_dev_b_hi(3)
     status=execute(code)
     yfield='calres_dev_b_hi(3)'
     plot, time_fex,data_cal, title= dataset_all, xrange=[0.,xmax*1.1],$
        ystyle=3, xstyle=1, xtitle= xlab, ytitle = yfield
 Endif
Endif
;
if(pltdev eq 'PS')then device,/close
if(pltdev eq 'TEK')then begin
   device,plot_to=0
   close,lunp
endif
;
return
end
