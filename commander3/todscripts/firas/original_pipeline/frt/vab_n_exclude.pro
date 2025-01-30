Pro VAB_N_EXCLUDE, DATE, CHAN, PLTDEV, F0_SDF, READ
;+
; Pro VAB_N_EXCLUDE, DATE, CHAN, PLTDEV, F0_SDF, READ
; Author:  Harte Wang, STX, August 1991
; Purpose:  Plot FDQ_SDF data that is outside the North Van Allen Belt, on
;     a terrestrial longitude & lattitude map.
; Parameters:
;         DATE: filename extension; date to check data, ie. '901100000'
;         CHAN: 'RH', 'RL', 'LH', or 'LL'
;         PLTDEV: plot device, 'PS' or 'REGIS'
;         F0_SDF: IDL structure of FDQ_SDF records
;         READ: Data read control ('READ' to read records). Not required
;               if F0_SDF structure has already been read during the
;               current IDL session.
;
; Input:
;       CSDR$FIRAS_REF:FEX_VABSAA.DAT
;	CSDR$FIRAS_OUT:FDQ_SDF_xx.ED_...
; Output:
;       F0_SDF: IDL structure of FDQ_SDF records.
;       Plot to screen or post script file VAB_N_EXCLUDE.PS in local directory.
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
print,"VAB_N_EXCLUDE, 'DATE', 'CHAN', 'PLTDEV', F0_SDF (, 'READ')"
print,' '
!P.MULTI(2) =1
chan=strupcase(chan)
pltdev = strupcase(pltdev)
set_plot, pltdev
If (pltdev eq 'PS') then device, filename='VAB_N_EXCLUDE.PS',/landscape
dataset_fex="CSDR$FIRAS_REF:FEX_VABSAA.DAT"
openr,lunr,dataset_fex,/get_lun,/share
finfo=fstat(lunr)
trec=finfo.size/finfo.rec_len
r0=assoc(lunr,fex_vabsaa_st(trec))
f0_fex=r0(0)
data_vab1y = f0_fex.vab(0).latn
data_vab1x = FIndgen(41)*f0_fex.vab(0).lonstep - 180.
data_vab2y = f0_fex.vab(0).lats
data_vab2x = FIndgen(41)*f0_fex.vab(0).lonstep - 180.
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
f0_tmp=f0_sdf(where(f0_sdf.attitude.terr_rad_byte ne 2))
;code ='data_fdq=f0_tmp'+'.'+field_fdq
;status=execute(code)
data_fdq_y = f0_tmp.attitude.terr_lat*.018/!PI		;  *180/pi/10000
data_fdq_x = f0_tmp.attitude.terr_long*.018/!PI
;time_fdq=conv_adt_sec(f0_tmp.ct_head.time)
;time_fdq=(time_fdq-time_fdq(0))/3600.
;emax=max(time_fdq) & imax=max(time_fdq)
;xmax=max(emax,imax)
plot, data_vab1x,data_vab1y,$
title=dataset_sdf,$
    yrange=[30,90],$
    xtitle=date + '--  TERRESTRIAL LONGITUDE ',$
    ytitle = ' TERRESTRIAL LATITUDE '
oplot, data_vab2x,data_vab2y
oplot, data_fdq_x,data_fdq_y,psym=3
;xrange=[0.,xmax*1.1],xstyle=1
;
if ( pltdev eq 'PS') then device,/close
set_plot,'regis'
;
return
end
