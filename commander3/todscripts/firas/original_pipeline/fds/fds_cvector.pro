;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FDS_CVECTOR extracts the C-Vector from the FDS sky IDL savesets.
;
;DESCRIPTION:
;     An IDL command file to extract the C-Vector from the FDS IDL savesets
;     and to write the C-Vector^2 to the binary file FEX_CVS_CCSS.VVV_XXXXXXX.
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
;	CSDR$FIRAS_IN    directory containing FDS IDL savesets
;	CSDR$FIRAS_OUT	 directory where FEX_CVS_CCSS file will be written
;
;     As written, FDS_CVECTOR.PRO will run on all FIRAS channels, but only if
;     Short-Slow, Short-Fast, or Long-Fast MTM scan modes are specified.
;
;EXAMPLE:
;     $ UIDL
;     !path = 'source directory,' + !path
;     setlog,'csdr$firas_in','x'
;     setlog,'csdr$firas_out','y'
;     chanscan='LHLF'
;     fds_cvector,chanscan,error
;     if(error ne 0)then print,'Error Returned from FDS_CVECTOR !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL function FDS_FEX_CVS to define the FEX_CVS record structure.
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 7 June 1993, SER 11039
;     Modified by Ken Jensen, Hughes STX, 8 June 1993  to change it
;       from a COM file to an IDL procedure for batch processing. The
;       global symbol CHANSCAN is now an input argument
;     Modified by Ken Jensen, Hughes STX, 8 June 1993  to return an error
;       status and to check for valid CHANSCAN.
;     Modified by Ken Jensen, Hughes STX, 8 June 1993  to re-define
;       unneeded restored parameters.
;       status and to check for valid CHANSCAN.
;     Modified by Ken Jensen, Hughes STX, 8 June 1993  to add PRINT
;       statements for the Batch Log File.
;     Modified by Gene Eplee, GSC, 29 September 1993 to add low frequency
;       channels.
;     Modified by Gene Eplee, GSC, 25 October 1993 to process the low frequency
;       short fast scan mode.  SER 11394
;     Modified by Gene Eplee, GSC, 27 October 1993 to write the C-Vector
;       into the FEX_CVS record structure and to write this record structure
;       to the FEX_CVS_CCSS binary file.  SER 11402
;     Modified by Gene Eplee, GSC, 4 February 1994 to process the gain/offset
;       destriper savesets.
;     Modified by Ken Jensen, Hughes STX, 31 March 1994  to re-define
;       additional unneeded restored parameters.
;     Modified by Ken Jensen, Hughes STX, 18 April 1994. Start and Stop
;       indices for RLLF and LLLF changed from (9:156) to (8:155), since
;       they are IDL indices, not Fortran.
;     Modified by Ken Jensen, Hughes STX, 10 May 1994. Modified the
;       re-definition of unused parameters to account for modified
;     Modified by Ken Jensen, Hughes STX, 20 June 1994. Modified the
;       re-definition of unused parameters to account for additional
;       saved arrays.
;-
;______________________________________________________________________________
;
Pro FDS_CVector,chanscan,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FDS_CVECTOR : Called Incorrectly : FDS_CVector,chanscan,error'
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
print,'FDS_CVECTOR  :  CHANSCAN='+chanscan
print,' '
;
;  Is CHANSCAN an acceptable string ? If not, Return with Error.
;
if (chan eq 'RH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FDS_CVECTOR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'RL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FDS_CVECTOR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FDS_CVECTOR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FDS_CVECTOR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else begin
  print,'FDS_CVECTOR : Invalid Channel'
  print,' '
  return
endelse
;
; Print Logical Translations.
;
print,'Logical Translations :'
ret=trnlog('csdr$firas_in',intrans,/full,/issue_error)
print,'CSDR$FIRAS_IN  == '+strupcase(intrans)
ret=trnlog('csdr$firas_out',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
print,' '
;
;  Determine which savesets to restore.
;
if (chanscan eq 'RHSS') then infile1 = 'csdr$firas_in:rhs.iss'
if (chanscan eq 'RHSS') then infile2 = 'csdr$firas_in:rhs_destriped.iss'
;
if (chanscan eq 'RHLF') then infile1 = 'csdr$firas_in:rhf.iss'
if (chanscan eq 'RHLF') then infile2 = 'csdr$firas_in:rhf_destriped.iss'
;
if (chanscan eq 'RLSS') then infile1 = 'csdr$firas_in:rls.iss'
if (chanscan eq 'RLSS') then infile2 = 'csdr$firas_in:rls_destriped.iss'
;
if (chanscan eq 'RLSF') then infile1 = 'csdr$firas_in:rsf.iss'
if (chanscan eq 'RLSF') then infile2 = 'csdr$firas_in:rsf_destriped.iss'
;
if (chanscan eq 'RLLF') then infile1 = 'csdr$firas_in:rlf.iss'
if (chanscan eq 'RLLF') then infile2 = 'csdr$firas_in:rlf_destriped.iss'
;
if (chanscan eq 'LHSS') then infile1 = 'csdr$firas_in:lhs.iss'
if (chanscan eq 'LHSS') then infile2 = 'csdr$firas_in:lhs_destriped.iss'
;
if (chanscan eq 'LHLF') then infile1 = 'csdr$firas_in:lhf.iss'
if (chanscan eq 'LHLF') then infile2 = 'csdr$firas_in:lhf_destriped.iss'
;
if (chanscan eq 'LLSS') then infile1 = 'csdr$firas_in:lls.iss'
if (chanscan eq 'LLSS') then infile2 = 'csdr$firas_in:lls_destriped.iss'
;
if (chanscan eq 'LLSF') then infile1 = 'csdr$firas_in:lsf.iss'
if (chanscan eq 'LLSF') then infile2 = 'csdr$firas_in:lsf_destriped.iss'
;
if (chanscan eq 'LLLF') then infile1 = 'csdr$firas_in:llf.iss'
if (chanscan eq 'LLLF') then infile2 = 'csdr$firas_in:llf_destriped.iss'
;
;  Restore the savesets.
;
print,'Restoring IDL Save Set '+strupcase(infile1)
print,' '
restore,infile1
print,'Restoring IDL Save Set '+strupcase(infile2)
print,' '
restore,infile2
;
; Re-Define UnNeeded Restored Parameters.
;
cal_glitch=1 & cal_sp=1 & cal_tm=1 & cal_wgts=1 & cal_lbl=1 & del_temp=1 & f=1
gn_step_dn=1 & gn_step_up=1 & glat=1 & glon=1 & ical=1 & nifgs=1
order=1 & px=1 & refh=1 & skyh=1 & sky_glitch=1 & sky_wgts=1 & sky_lbl=1 & sp=1
step_dn=1 & step_up=1 & tm=1 & xcal=1 & ndf_off=1 & ndf_none=1 & csp=1 & index=1
dihd=0 & sky_wgts_ds=0 & ndf_ejv=0 & sqr_mat=0 & rect=0 & diag=0 & c_calsp=0
cal_wgts_ds=0 & pcvr=0 & vib=0 & corr_index=0
;
; Re-Define UnNeeded Mode-Specific Restored Parameters.
;
if(chanscan eq 'RHSS')then begin
 b_rhs=1 & d_rhs=1 & ejv_rhs=1 & l_rhs=1 & n_rhs=1 & s_rhs=1 & tau_rhs=1 & urhs=1
endif
if(chanscan eq 'RHLF')then begin
 b_rhf=1 & d_rhf=1 & ejv_rhf=1 & l_rhf=1 & n_rhf=1 & s_rhf=1 & tau_rhf=1 & urhf=1
 d_rhsf=1 & d_rhlf=1 & sl_weight_rat=1
endif
if(chanscan eq 'RLSS')then begin
 b_rls=1 & d_rls=1 & ejv_rls=1 & l_rls=1 & n_rls=1 & s_rls=1 & tau_rls=1 & urls=1
endif
if(chanscan eq 'RLSF')then begin
 b_rsf=1 & d_rsf=1 & ejv_rsf=1 & l_rsf=1 & n_rsf=1 & s_rsf=1 & tau_rsf=1 & ursf=1
 d_rlsf=1 & d_rlfl=1 & sl_weight_rat=1
endif
if(chanscan eq 'RLLF')then begin
 b_rlf=1 & d_rlf=1 & ejv_rlf=1 & l_rlf=1 & n_rlf=1 & s_rlf=1 & tau_rlf=1 & urlf=1
endif
if(chanscan eq 'LHSS')then begin
 b_lhs=1 & d_lhs=1 & ejv_lhs=1 & l_lhs=1 & n_lhs=1 & s_lhs=1 & tau_lhs=1 & ulhs=1
endif
if(chanscan eq 'LHLF')then begin
 b_lhf=1 & d_lhf=1 & ejv_lhf=1 & l_lhf=1 & n_lhf=1 & s_lhf=1 & tau_lhf=1 & ulhf=1
 d_lhsf=1 & d_lhlf=1 & sl_weight_rat=1
endif
if(chanscan eq 'LLSS')then begin
 b_lls=1 & d_lls=1 & ejv_lls=1 & l_lls=1 & n_lls=1 & s_lls=1 & tau_lls=1 & ulls=1
endif
if(chanscan eq 'LLSF')then begin
 b_lsf=1 & d_lsf=1 & ejv_lsf=1 & l_lsf=1 & n_lsf=1 & s_lsf=1 & tau_lsf=1 & ulsf=1
 d_llsf=1 & d_llfl=1 & sl_weight_rat=1
endif
if(chanscan eq 'LLLF')then begin
 b_llf=1 & d_llf=1 & ejv_llf=1 & l_llf=1 & n_llf=1 & s_llf=1 & tau_llf=1 & ullf=1
endif
;
;  Extract the real part of the C-Vector and square the result.
;
print,'Converting the C-Vector to Real.'
print,' '
if (chanscan eq 'RHSS') then c = (crhs(*,0))^2
if (chanscan eq 'RHLF') then c = (crhf(*,0))^2
if (chanscan eq 'RLSS') then c = (crls(*,0))^2
if (chanscan eq 'RLSF') then c = (crsf(*,0))^2
if (chanscan eq 'RLLF') then c = (crlf(*,0))^2
if (chanscan eq 'LHSS') then c = (clhs(*,0))^2
if (chanscan eq 'LHLF') then c = (clhf(*,0))^2
if (chanscan eq 'LLSS') then c = (clls(*,0))^2
if (chanscan eq 'LLSF') then c = (clsf(*,0))^2
if (chanscan eq 'LLLF') then c = (cllf(*,0))^2
;
;  Define and fill in the FEX_CVS record.
;
print,'Filling in the FEX_CVS record.'
print,' '
fex_rec = fds_fex_cvs(1)
;
;  Get the run time information
;
ctime=!stime
ctime=strmid(ctime,0,11)+':'+strmid(ctime,12,11)
ztime=timeconv(ctime,infmt='v',outfmt='z')
atime=timeconv(ctime,infmt='v',outfmt='a')
fex_rec.gmt = ztime
fex_rec.time(0)=atime(0)
fex_rec.time(1)=atime(1)
;
;  Set the scanmode information
;
if (chanscan eq 'RHSS') then begin
   fex_rec.channel     = 1
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 0
endif
if (chanscan eq 'RHLF') then begin
   fex_rec.channel     = 1
   fex_rec.scan_length = 1
   fex_rec.scan_speed  = 1
endif
if (chanscan eq 'RLSS') then begin
   fex_rec.channel     = 2
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 0
endif
if (chanscan eq 'RLSF') then begin
   fex_rec.channel     = 2
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 1
endif
if (chanscan eq 'RLLF') then begin
   fex_rec.channel     = 2
   fex_rec.scan_length = 1
   fex_rec.scan_speed  = 1
endif
if (chanscan eq 'LHSS') then begin
   fex_rec.channel     = 3
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 0
endif
if (chanscan eq 'LHLF') then begin
   fex_rec.channel     = 3
   fex_rec.scan_length = 1
   fex_rec.scan_speed  = 1
endif
if (chanscan eq 'LLSS') then begin
   fex_rec.channel     = 4
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 0
endif
if (chanscan eq 'LLSF') then begin
   fex_rec.channel     = 4
   fex_rec.scan_length = 0
   fex_rec.scan_speed  = 1
endif
if (chanscan eq 'LLLF') then begin
   fex_rec.channel     = 4
   fex_rec.scan_length = 1
   fex_rec.scan_speed  = 1
endif
;
;  Set the remaining header information
;
fex_rec.model_label = string(solution,format='(A40)')
fex_rec.galat_exc   = float(galcut)
fex_rec.ncal_ifgs   = total(cal_nifgs)
;
;  Get the frequency range of the C-Vector and
;  write the C-Vector to the FEX record
;
if (strmid(chanscan,1,1) eq 'H') then begin
   jstart = 4
   jstop  = 170
endif else if (strmid(chanscan,1,2) eq 'LS') then begin
   jstart = 4
   jstop  = 37
endif else if (strmid(chanscan,1,2) eq 'LL') then begin
   jstart = 8
   jstop  = 155
endif
for  j=jstart,jstop do fex_rec.cvector(j) = c(j-jstart)
;
;  Write the FEX_CVS file.
;
print,'Writing the FEX_CVS file.'
print,' '
fext = solution
i = strpos(fext,' ')
s1 = strtrim(strmid(fext,i,35),2)
fext = strtrim(s1 + '_' + strmid(fext,0,i),2)
outfile = 'csdr$firas_out:fex_cvs_' + chanscan + '.' + fext
rec_len = 2176                 ; fixed length record size in bytes
openw,1,outfile, rec_len, /fixed
writeu,1,fex_rec
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
