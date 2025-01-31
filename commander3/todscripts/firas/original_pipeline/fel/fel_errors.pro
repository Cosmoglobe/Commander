;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FEL_ERRORS extracts the FIRAS calibration errors from their IDL savesets.
;
;DESCRIPTION:
;     An IDL procedure to extract the FIRAS calibration errors from their
;     various IDL savesets and to write the these errors to the binary file 
;     FEL_ERR_CCSS.XXXXX.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS
;     CHANSCAN  (I)  :  String defining channel and scan mode (e.g. 'RHSS')
;     ERROR     (O)  :  Error status for procedure
;
;WARNINGS:
;     The IDL path must be set before invoking this procedure.
;
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_PEP   directory containing PEP error IDL savesets
;	CSDR$FIRAS_JCJ   directory containing JCJ error IDL savesets
;	CSDR$FIRAS_PUP   directory containing PUP IDL savesets
;	CSDR$FIRAS_PTP   directory containing PTP IDL savesets
;	CSDR$FIRAS_OUT	 directory where FEL_ERR error file will be written
;
;EXAMPLE:
;     $ UIDL
;     setlog,'csdr$firas_pep','b'
;     setlog,'csdr$firas_jcj','c'
;     setlog,'csdr$firas_pup','d'
;     setlog,'csdr$firas_ptp','e'
;     setlog,'csdr$firas_out','f'
;     chanscan='RHSS'
;     fel_errors,chanscan,error
;     if(error ne 0)then print,'Error Returned from FEL_ERRORS !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL function FEL_ERR_STR to define the FEL_ERR record structure.
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 26 August 1994, SER 11894
;     Modified by S. Brodd, HSTX, 6/5/96, for Pass 4, SER 12332
;     Modified by S. Brodd, HSTX, 6/20/97, remove C and D vectors, SPR 12350.
;-
;______________________________________________________________________________
;
Pro FEL_Errors,chanscan,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FEL_ERRORS : Called Incorrectly : FEL_Errors,chanscan,error'
 print,' '
 return
endif
;
; Define CHAN and SCAN strings, and make them uppercase, as well as CHANSCAN
;
chan=strupcase(strmid(chanscan,0,2))
scan=strupcase(strmid(chanscan,2,2))
chanscan=chan+scan
;
print,' '
print,'FEL_ERRORS  :  CHANSCAN='+chanscan
print,' '
;
;  Is CHANSCAN an acceptable string ? If not, Return with Error.
;
if (chan eq 'RH') then begin
  if ((scan ne 'SS')and(scan ne 'LF')) then begin
    print,'FEL_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'RL') then begin
  if ((scan ne 'SS')and(scan ne 'FS')and(scan ne 'LF')) then begin
    print,'FEL_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LH') then begin
  if ((scan ne 'SS')and(scan ne 'LF')) then begin
    print,'FEL_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LL') then begin
  if ((scan ne 'SS')and(scan ne 'FS')and(scan ne 'LF')) then begin
    print,'FEL_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else begin
  print,'FEL_ERRORS : Invalid Channel'
  print,' '
  return
endelse
;
; Print Logical Translations.
;
print,'Logical Translations :'
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
;  Initialize the FEL_ERR record.
;
print,'Initializing the FEL_ERR record.'
print,' '
err_rec = fel_err_str(1)
;
;  Define channel/scan-mode specific parameters
;
if (strmid(chanscan,1,1) eq 'H') then begin
   stop = 209
   max_pem = 12
endif else if (strmid(chanscan,1,2) eq 'LL') then begin
   stop = 181
   max_pem = 10
endif else begin
   stop = 42
   max_pem = 10
endelse
;
;  Get the PEP error.
;
print,'Getting the PEP Errors.'
print,' '
restore,file='csdr$firas_pep:pep_errors.'+chanscan
err_rec.pep_gain(0:stop)   = gain
err_rec.pep_offset(0:stop) = pep
;
;  Get the JCJ errors.
;
print,'Getting the JCJ Errors.'
print,' '
restore,file='csdr$firas_jcj:jcj_spec.'+chanscan
for j=0,max_pem do err_rec.jcj_gain(j,0:stop) = float(ag(j,*))
for j=0,max_pem do err_rec.jcj_offset(j,0:stop) = float(ao(j,*))
;
;  Get the PUP error.
;
print,'Getting the PUP Error.'
print,' '
restore,file='csdr$firas_pup:pup_errors.'+chanscan
err_rec.pup_temp         = float(dt)
err_rec.pup_spec(0:stop) = float(pup)
;
;  Get the PTP error.
;
print,'Getting the PTP Error.'
print,' '
restore,file='csdr$firas_ptp:ptp_errors.'+chanscan
err_rec.ptp_temp         = float(dt)
err_rec.ptp_spec(0:stop) = float(ptp)
;
;  Deassign the unneeded variables.
;
dbdt=0 & dp=0 & f=0 & gl=0 & idh=0 & idl=0 & jcj=0 & pixel=0
;
;  Write the FEL_ERR file.
;
print,'Writing the FEL_ERR file.'
print,' '
outfile = 'csdr$firas_out:fel_err_' + chanscan + '.pass4'
rec_len = 25208                 ; fixed length record size in bytes
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
