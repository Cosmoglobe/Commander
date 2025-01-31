Pro ENG_IDX_PLOT, DATE, FIELD_ENG, FIELD_IDX, PLTDEV, F0_ENG, F0_IDX, READ
;+
; Pro ENG_IDX_PLOT, DATE, FIELD_ENG, FIELD_IDX, PLTDEV, F0_ENG, F0_IDX, READ
; Author:  Harte Wang, STX, August 1991
; Purpose:  Plot engineering field values from FDQ_ENG and plot corresponding
;    field in FDQ_IDX file to show that FDQ correctly generates an IDX record
;    whenever a field value changes.
; Parameters:
;         DATE: day to check data, ie. '901100000'
;         FIELD_ENG: engineering field, ie. 'EN_ANALOG.TEMP_CTRL(7)'
;         FIELD_IDX: index to field, ie.  'TEMP_CTRL.GROUP2(7)'
;         PLTDEV: plot device, 'PS' or 'REGIS'
;         F0_ENG: IDL structure of FDQ_ENG records
;         F0_IDX: IDL structure of FDQ_IDX records
;         READ: Data read control ('READ' to read records). Not required
;               if F0_ENG and F0_IDX structures have already been read
;               during the current IDL session.
;
; Input:
;       CSDR$FIRAS_OUT:FDQ_ENG.ED_...
;       CSDR$FIRAS_OUT:FDQ_IDX.ED_...
; Output:
;       F0_ENG: IDL structure of FDQ_ENG records
;       F0_IDX: IDL structure of FDQ_IDX records
;       Plot to screen or post script file ENG_IDX_PLOT.PS in local directory.
;
; Modifications: K.Jensen (STX) - Oct 7, 1991 - Added F0_ENG and F0_IDX
;      as Command Line Parameters. Program now returns these as output.
;      Added 'READ' option to command line. Program can now skip the data
;      read if the F0_ENG and F0_IDX structures already exist in memory.
;
;      SPR 9170 - K.Jensen (STX) - Oct 28, 1991 - Correct read statment.
;-
print,' '
print,"ENG_IDX_PLOT,'DATE','FIELD_ENG','FIELD_IDX','PLTDEV',F0_ENG,F0_IDX,('READ')"
print,' '
!P.MULTI(2) = 3
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='ENG_IDX_PLOT.PS',/landscape
DATE = STRMID (DATE, 0, 5) + '0000'
dataset="CSDR$FIRAS_OUT:FDQ_ENG.ED_"+DATE
dataset_idx="CSDR$FIRAS_OUT:FDQ_IDX.ED_"+DATE
dataset_idx2="FDQ_IDX.ED_"+DATE
nread='NOREAD'
RS = SIZE (READ)
IF (RS(1) NE 0) THEN NREAD=READ
if(strupcase(nread) eq 'READ')then begin
 openr,lunr,dataset,/get_lun,/share
 finfo=fstat(lunr)
 trec=finfo.size/finfo.rec_len
 r1=assoc(lunr,fdq_eng_st(trec))
 print,' '
 print,strupcase(dataset)
 srec1=strcompress(string(trec))
 print,'Reading'+srec1+' FDQ_ENG Records.'
 print,' '
 f0_eng=r1(0)
 openr,lunr_idx,dataset_idx,/get_lun,/share
 finfo_idx=fstat(lunr_idx)
 trec_idx=finfo_idx.size/finfo_idx.rec_len
 r2=assoc(lunr_idx,fdq_idx_st(trec_idx))
 print,' '
 print,strupcase(dataset_idx)
 srec1=strcompress(string(trec_idx))
 print,'Reading'+srec1+' FDQ_IDX Records.'
 print,' '
 f0_idx = r2(0)
endif
code ='data_eng=f0_eng'+'.'+field_eng
status=execute(code)
time_eng=conv_adt_sec(f0_eng.ct_head.time)
code='data_IDX=f0_idx'+'.'+field_idx
status=execute(code)
time_idx=conv_adt_sec(f0_idx.header.bnry_start_time)
time_eng=(time_eng-time_eng(0))/3600.
time_idx=(time_idx-time_idx(0))/3600.
emax=max(time_eng) & imax=max(time_idx)
xmax=max(emax,imax)
plot, time_eng,data_eng,xrange=[0.,xmax*1.1],xstyle=1,$
title=dataset,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_eng
plot, time_idx,data_idx,xrange=[0.,xmax*1.1],xstyle=1,$
title=dataset_idx,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_idx
plot, time_eng,data_eng,xrange=[0.,xmax*1.1],xstyle=1,$
title=dataset+ ' & ' +dataset_idx2,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_eng + ' & ' + field_idx
oplot, time_idx,data_idx,linestyle=1
;oplot, data2,yrange=[2.5,3.0]
;
If ( pltdev eq 'PS' ) then device,/close
;
return
end
