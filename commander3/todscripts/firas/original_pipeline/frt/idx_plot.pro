Pro IDX_PLOT, DATE, FIELD_IDX, PLTDEV, F0_IDX, READ
;+
; Pro IDX_PLOT, DATE, FIELD_IDX, PLTDEV, F0_IDX, READ
; Author:  Harte Wang, STX, August 1991
; Purpose:  Plot engineering field values from FDQ_IDX.
; Parameters:
;         DATE: file extension; day to check data, ie. '901100000'
;         FIELD_IDX: index to field, ie.  'TEMP_CTRL.GROUP2(7)'
;         PLTDEV: plot device, 'PS' or 'REGIS'
;         F0_IDX: IDL structure of FDQ_IDX records
;         READ: Data read control ('READ' to read records). Not required
;               if F0_IDX structure has already been read during the
;               current IDL session.
;
; Input:
;       CSDR$FIRAS_OUT:FDQ_IDX.ED_...
; Output:
;       F0_IDX: IDL structure of FDQ_IDX records.
;       Plot to screen or post script file IDX_PLOT.PS in local directory.
;
; Modifications: K.Jensen (STX) - Oct 7, 1991 - Added F0_IDX to command line.
;      Program now returns this as output.
;      Added 'READ' option to command line. Program can now skip the data
;      read if the F0_IDX structure already exists in memory.
;
;      SPR 9170 - K.Jensen (STX) - Oct 28, 1991 - Correct read statment.
;-
print,' '
print,"IDX_PLOT, 'DATE', 'FIELD_IDX', 'PLTDEV', F0_IDX (,'READ')"
print,' '
!p.multi=0
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='IDX_PLOT.PS',/landscape
DATE = STRMID (DATE, 0, 5) + '0000'
dataset_idx="CSDR$FIRAS_OUT:FDQ_IDX.ED_"+DATE
nread='NOREAD'
RS = SIZE (READ)
IF (RS(1) NE 0) THEN NREAD=READ
if(strupcase(nread) eq 'READ')then begin
 openr,lunr_idx,dataset_idx,/get_lun,/share
 finfo_idx=fstat(lunr_idx)
 trec_idx=finfo_idx.size/finfo_idx.rec_len
 r2=assoc(lunr_idx,fdq_idx_st(trec_idx))
 srec1=strcompress(string(trec_idx))
 print,' '
 print,strupcase(dataset_idx)
 print,'Reading'+srec1+' FDQ_IDX Records.'
 print,' '
 f0_idx = r2(0)
endif
code='data_IDX=f0_idx'+'.'+field_idx
status=execute(code)
time_idx=conv_adt_sec(f0_idx.header.bnry_start_time)
time_idx=(time_idx-time_idx(0))/3600.
emax=max(time_idx) & imax=max(time_idx)
xmax=max(emax,imax)
plot, time_idx,data_idx,xrange=[0.,xmax*1.1],xstyle=1,$
title=dataset_idx,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_idx
;
if ( pltdev eq 'PS') then device,/close
set_plot,'regis'
;
return
end
