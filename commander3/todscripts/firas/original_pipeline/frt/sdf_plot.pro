Pro SDF_PLOT, DATE, CHAN, FIELD_SDF, PLTDEV, F0_SDF, READ
;+
; Pro SDF_PLOT, DATE, CHAN, FIELD_SDF, PLTDEV, F0_SDF, READ
; Author:  Harte Wang, STX, August 1991
; Purpose:  Plot field values from FDQ_SDF_<xx> for the date.
; Parameters:
;         DATE: day to check data, ie. '901100000'
;         CHAN: 'RH','RL','LH', or 'LL'
;         FIELD_SDF: field in FDQ_SDF, ie. 'ATTITUDE.PIXEL_NO'
;         PLTDEV: plot device, 'PS' or 'REGIS'
;         F0_SDF: IDL structure of FDQ_SDF records
;         READ: Data read control ('READ' to read records). Not required
;               if F0_SDF structure has already been read during the
;               current IDL session.
;
; Input:
;       CSDR$FIRAS_OUT:FDQ_SDF_<xx>.ED_...
; Output:
;       F0_SDF: IDL structure of FDQ_SDF records.
;       Plot to screen or post script file SDF_PLOT.PS in local directory.
;
; Modifications: K.Jensen (STX) - Oct 7, 1991 - Added F0_SDF to command line.
;      Program now returns this as output.
;      Added 'READ' option to command line. Program can now skip the data
;      read if the F0_SDF structure already exists in memory.
;
; Note: Excludes data with bad midpoint of collect time.
;
;      SPR 9170 - K.Jensen (STX) - Oct 28, 1991 - Correct read statment.
;-
print,' '
print,"SDF_PLOT, 'DATE', 'CHAN', 'FIELD_SDF', 'PLTDEV', F0_SDF (,'READ')"
print,' '
!p.multi=0
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='SDF_PLOT.PS',/landscape
DATE = STRMID (DATE, 0, 5) + '0000'
dataset_sdf="CSDR$FIRAS_OUT:FDQ_SDF_"+CHAN+".ED_"+DATE
nread='NOREAD'
RS = SIZE (READ)
IF (RS(1) NE 0) THEN NREAD=READ
if(strupcase(nread) eq 'READ')then begin
 openr,lunr,dataset_sdf,/get_lun,/share
 finfo=fstat(lunr)
 trec=finfo.size/finfo.rec_len
 r1=assoc(lunr,fdq_sdf_st(trec))
 srec1=strcompress(string(trec))
 print,' '
 print,strupcase(dataset_sdf)
 print,'Reading'+srec1+' FDQ_SDF Records.'
 print,' '
 f0_sdf=r1(0)
 f0_sdf=f0_sdf(where(f0_sdf.collect_time.badtime_flag eq 0))
endif
code ='data_sdf=f0_sdf'+'.'+field_sdf
status=execute(code)
time_sdf=conv_adt_sec(f0_sdf.ct_head.time)
time_sdf=(time_sdf-time_sdf(0))/3600.
emax=max(time_sdf) & imax=max(time_sdf)
xmax=max(emax,imax)
plot, time_sdf,data_sdf,xrange=[0.,xmax*1.1],xstyle=1,$
    title="FDQ_SDF_"+STRUPCASE(CHAN)+".ED_"+DATE,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_sdf + ' -- ' + strupcase(chan)
;
if ( pltdev eq 'PS') then device,/close
set_plot,'regis'
;
return
end
