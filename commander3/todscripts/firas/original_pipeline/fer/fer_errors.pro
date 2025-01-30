;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FER_ERRORS extracts the FIRAS calibration errors from their IDL savesets.
;
;DESCRIPTION:
;     An IDL procedure to extract the FIRAS calibration errors from their
;     various IDL savesets and to write the these errors to the binary file 
;     FER_ERR_CCSS.VVV_XXXXXXX.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS
;     CHANSCAN  (I)  :  String defining channel and scan mode (e.g. 'RHSS')
;     ERROR     (O)  :  Error status for procedure
;
;WARNINGS:
;     For the merged short fast and long fast scan modes, the program reads the
;     FMS_VAR_CCSS.VVV_XXXXXXXX file rather than the FEX_VAR_CCSS.VVV_XXXXXXXX
;     file.
;
;     The IDL path must be set before invoking this procedure.
;
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_VAR   directory containing FEX_VAR binary files
;	CSDR$FIRAS_CVS   directory containing FEX_CVS binary files
;	CSDR$FIRAS_PEP   directory containing PEP error IDL savesets
;	CSDR$FIRAS_JCJ   directory containing JCJ error IDL savesets
;	CSDR$FIRAS_PUP   directory containing PUP IDL savesets
;	CSDR$FIRAS_PTP   directory containing PTP IDL savesets
;	CSDR$FIRAS_OUT	 directory where FER_ERR error file will be written
;
;EXAMPLE:
;     $ UIDL
;     !path = 'source directory,' + !path
;     setlog,'csdr$firas_var','a'
;     setlog,'csdr$firas_cvs','b'
;     setlog,'csdr$firas_pep','c'
;     setlog,'csdr$firas_jcj','d'
;     setlog,'csdr$firas_pup','e'
;     setlog,'csdr$firas_ptp','f'
;     setlog,'csdr$firas_out','g'
;     chanscan='RHSS'
;     fer_errors,chanscan,error
;     if(error ne 0)then print,'Error Returned from FER_ERRORS !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL function FER_ERR_STR to define the FER_ERR record structure.
;     Calls IDL functions FER_GET_VAR and FER_FEX_VAR to get the D-Vector.
;     Calls IDL functions FER_GET_CVS and FER_FEX_CVS to get the C-Vector.
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 26 August 1994, SER 11894
;-
;______________________________________________________________________________
;
Pro FER_Errors,chanscan,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FER_ERRORS : Called Incorrectly : FER_Errors,chanscan,error'
 print,' '
 return
endif
;
; Define CHAN and SMODE strings, and make them uppercase, as well as CHANSCAN
;
chan=strupcase(strmid(chanscan,0,2))
smode=strupcase(strmid(chanscan,2,2))
chanscan=chan+smode
;
print,' '
print,'FER_ERRORS  :  CHANSCAN='+chanscan
print,' '
;
;  Is CHANSCAN an acceptable string ? If not, Return with Error.
;
if (chan eq 'RH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FER_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'RL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FER_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FER_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FER_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else begin
  print,'FER_ERRORS : Invalid Channel'
  print,' '
  return
endelse
;
; Print Logical Translations.
;
print,'Logical Translations :'
ret=trnlog('csdr$firas_var',intrans,/full,/issue_error)
print,'CSDR$FIRAS_VAR  == '+strupcase(intrans)
ret=trnlog('csdr$firas_cvs',intrans,/full,/issue_error)
print,'CSDR$FIRAS_CVS  == '+strupcase(intrans)
ret=trnlog('csdr$firas_pep',intrans,/full,/issue_error)
print,'CSDR$FIRAS_PEP  == '+strupcase(intrans)
ret=trnlog('csdr$firas_jcj',intrans,/full,/issue_error)
print,'CSDR$FIRAS_JCJ  == '+strupcase(intrans)
ret=trnlog('csdr$firas_pup',intrans,/full,/issue_error)
print,'CSDR$FIRAS_PUP  == '+strupcase(intrans)
ret=trnlog('csdr$firas_ptp',intrans,/full,/issue_error)
print,'CSDR$FIRAS_PTP  == '+strupcase(intrans)
ret=trnlog('csdr$firas_out',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
print,' '
;
;  Initialize the FER_ERR record.
;
print,'Initializing the FER_ERR record.'
print,' '
err_rec = fer_err_str(1)
if (chanscan eq 'RHSS') then begin
   err_rec.channel     = 1
   err_rec.scan_length = 0
   err_rec.scan_speed  = 0
   start = 4
   stop  = 170
endif
if (chanscan eq 'RHLF') then begin
   err_rec.channel     = 1
   err_rec.scan_length = 1
   err_rec.scan_speed  = 1
   start  = 4
   stop   = 170
endif
if (chanscan eq 'RLSS') then begin
   err_rec.channel     = 2
   err_rec.scan_length = 0
   err_rec.scan_speed  = 0
   start = 4
   stop  = 37
endif
if (chanscan eq 'RLSF') then begin
   err_rec.channel     = 2
   err_rec.scan_length = 0
   err_rec.scan_speed  = 1
   start  = 4
   stop   = 37
endif
if (chanscan eq 'RLLF') then begin
   err_rec.channel     = 2
   err_rec.scan_length = 1
   err_rec.scan_speed  = 1
   start = 8
   stop  = 155
endif
if (chanscan eq 'LHSS') then begin
   err_rec.channel     = 3
   err_rec.scan_length = 0
   err_rec.scan_speed  = 0
   start = 4
   stop  = 170
endif
if (chanscan eq 'LHLF') then begin
   err_rec.channel     = 3
   err_rec.scan_length = 1
   err_rec.scan_speed  = 1
   start  = 4
   stop   = 170
endif
if (chanscan eq 'LLSS') then begin
   err_rec.channel     = 4
   err_rec.scan_length = 0
   err_rec.scan_speed  = 0
   start = 4
   stop  = 37
endif
if (chanscan eq 'LLSF') then begin
   err_rec.channel     = 4
   err_rec.scan_length = 0
   err_rec.scan_speed  = 1
   start  = 4
   stop   = 37
endif
if (chanscan eq 'LLLF') then begin
   err_rec.channel     = 4
   err_rec.scan_length = 1
   err_rec.scan_speed  = 1
   start = 8
   stop  = 155
endif
;
;  Get the run time information
;
ctime=!stime
ctime=strmid(ctime,0,11)+':'+strmid(ctime,12,11)
ztime=timeconv(ctime,infmt='v',outfmt='z')
atime=timeconv(ctime,infmt='v',outfmt='a')
err_rec.gmt = ztime
err_rec.time(0)=atime(0)
err_rec.time(1)=atime(1)

;
;  Get the PEP error.
;
print,'Getting the PEP Errors.'
print,' '
restore,file='csdr$firas_pep:pep_errors.'+chanscan
err_rec.pep_gain(start:stop)   = gain
err_rec.pep_offset(start:stop) = pep
;
;  Get the JCJ errors.
;
print,'Getting the JCJ Errors.'
print,' '
restore,file='csdr$firas_jcj:jcj_spec.'+chanscan
for j=0,13 do err_rec.jcj_gain(j,start:stop) = float(ag(j,*))
for j=0,13 do err_rec.jcj_offset(j,start:stop) = float(ao(j,*))
;
;  Get the PUP error.
;
print,'Getting the PUP Error.'
print,' '
restore,file='csdr$firas_pup:pup_errors.'+chanscan
err_rec.pup_temp             = float(dt)
err_rec.pup_spec(start:stop) = float(pup)
;
;  Get the PTP error.
;
print,'Getting the PTP Error.'
print,' '
restore,file='csdr$firas_ptp:ptp_errors.'+chanscan
err_rec.ptp_temp             = float(dt)
err_rec.ptp_spec(start:stop) = float(ptp)

;
;  Deassign the unneeded variables.
;
dbdt=0 & dp=0 & f=0 & gl=0 & idh=0 & idl=0 & jcj=0 & pixel=0
;
;  Get the model solution
;
i = strpos(model,' ')
s1 = strtrim(strmid(model,i,35),2)
s1 = strcompress(s1)
i = strpos(s1,' ')
len = strlen(s1)
solution = strtrim(strmid(s1,i+1,len-1)+'_'+strmid(s1,0,i-1),2)

;
;  Get the D-Vector from the FEX_VAR file
;
print,'Getting the D-Vector.'
print,' '
fer_get_var,chanscan,solution,var_rec
err_rec.model_label = string(var_rec.model_label)
err_rec.galatexc    = var_rec.galat_exc
err_rec.d_vector(start:stop) = sqrt(var_rec.variances(start:stop)) 
;
;  Get the C-Vector from the FEX_CVS file
;
print,'Getting the C-Vector.'
print,' '
fer_get_cvs,chanscan,solution,cvs_rec
err_rec.c_vector(start:stop) = sqrt(cvs_rec.cvector(start:stop)) 

;
;  Write the FER_ERR file.
;
print,'Writing the FER_ERR file.'
print,' '
outfile = 'csdr$firas_out:fer_err_' + chanscan + '.' + solution
rec_len = 29184                 ; fixed length record size in bytes
openw,1,outfile, rec_len, /fixed
writeu,1,err_rec
close,1
;
print,'File "'+strupcase(outfile)+'" Written.'
print,' '
;
; Set Error Status to NO Error.
;
error=0
;
return
end
