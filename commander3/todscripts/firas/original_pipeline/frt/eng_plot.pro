Pro ENG_PLOT, DATE, FIELD_ENG, PLTDEV, F0_ENG, READ
;+
; Pro ENG_PLOT, DATE, FIELD_ENG, PLTDEV, F0_ENG, READ
; Author:  Harte Wang, STX, August 1991
; Purpose:  Plot engineering field values from FDQ_ENG.
; Parameters:
;         DATE: day to check data, ie. '901100000'
;         FIELD_ENG: engineering field, ie. 'EN_ANALOG.TEMP_CTRL(7)'
;         PLTDEV: plot device, 'PS' or 'REGIS'
;         F0_ENG: IDL structure of FDQ_ENG records
;         READ: Data read control ('READ' to read records). Not required
;               if F0_ENG structure has already been read during the
;               current IDL session.
;
; Input:
;       CSDR$FIRAS_OUT:FDQ_ENG.ED_...
; Output:
;       F0_ENG: IDL structure of FDQ_ENG records.
;       Plot to screen or post script file ENG_PLOT.PS in local directory.
;
; Modifications: K.Jensen (STX) - Oct 7, 1991 - Added F0_ENG to command line.
;      Program now returns this as output.
;      Added 'READ' option to command line. Program can now skip the data
;      read if the F0_ENG structure already exists in memory.
;
;      SPR 9170 - K.Jensen (STX) - Oct 28, 1991 - Correct read statment.
;-
print,' '
print,"ENG_PLOT, 'DATE', 'FIELD_ENG', 'PLTDEV', F0_ENG (,'READ')"
print,' '
!p.multi=0
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='ENG_PLOT.PS',/landscape
DATE = STRMID (DATE, 0, 5) + '0000'
dataset="CSDR$FIRAS_OUT:FDQ_ENG.ED_"+DATE
nread='NOREAD'
RS = SIZE (READ)
IF (RS(1) NE 0) THEN NREAD=READ
if(strupcase(nread) eq 'READ')then begin
 openr,lunr,dataset,/get_lun,/share
 finfo=fstat(lunr)
 trec=finfo.size/finfo.rec_len
 r1=assoc(lunr,fdq_eng_st(trec))
 srec1=strcompress(string(trec))
 print,' '
 print,strupcase(dataset)
 print,'Reading'+srec1+' FDQ_ENG Records.'
 print,' '
 f0_eng=r1(0)
 f0_eng = f0_eng(where(f0_eng.en_head.dq_sum_flag(0) ne 127 and $
          f0_eng.en_head.dq_sum_flag(1) ne 127 and $
          f0_eng.en_head.dq_sum_flag(2) ne 127 and $
          f0_eng.en_head.dq_sum_flag(3) ne 127))
endif
code ='data_eng=f0_eng'+'.'+field_eng
status=execute(code)
data_eng = data_eng(where(data_eng gt -9999))
time_eng=conv_adt_sec(f0_eng.ct_head.time)
time_eng=(time_eng-time_eng(0))/3600.
emax=max(time_eng) & imax=max(time_eng)
xmax=max(emax,imax)
plot, time_eng,data_eng,xrange=[0.,xmax*1.1],xstyle=1,$
title=dataset,$
    xtitle=date + ' -- ' + strmid(date,0,5) +'235959990',$
    ytitle = field_eng
;
if ( pltdev eq 'PS') then device,/close
set_plot,'regis'
;
return
end
