;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FDS_ERRORS extracts the Destriper Errors from the destriped IDL savesets.
;
;DESCRIPTION:
;     An IDL command file to extract the Destriper Error Matrices from the FDS
;     IDL savesets and to write the these matrices to the binary files 
;     FIP_DSE_CCSS.VVV_XXXXXXX and FIP_DSQ_CCSS.VVV_XXXXXXX.
;     The command file also calls FDS_ERRMAT to compute the auxiliary FDS
;     weights and covariance matrices.
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
;	CSDR$FIRAS_OUT	 directory where FIP error matrix files will be written
;
;EXAMPLE:
;     $ UIDL
;     !path = 'source directory,' + !path
;     setlog,'csdr$firas_in','x'
;     setlog,'csdr$firas_out','y'
;     chanscan='LHLF'
;     fds_errors,chanscan,error
;     if(error ne 0)then print,'Error Returned from FDS_ERRORS !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL procedure FDS_ERRMAT to compute the auxiliary matrices.
;     Calls IDL function FDS_FIP_DSE to define the FIP_DSE record structure.
;     Calls IDL function FDS_FIP_DSQ to define the FIP_DSQ record structure.
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 4 August 1994, SER 11858
;
;     Modified by Ken Jensen, Hughes STX, 9 August 1994, STRIPE_CONTRIB field
;     set to zero for all records where the total SKY_WGTS_DS weight = 0.
;
;     Modified by Gene Eplee, GSC, 9 September 1994 to compute the auxiliary
;     matrices and write them out to the FIP_DSQ and FIP_DSE matrix files,
;     SER 11893
;-
;______________________________________________________________________________
;
Pro FDS_Errors,chanscan,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FDS_ERRORS : Called Incorrectly : FDS_Errors,chanscan,error'
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
print,'FDS_ERRORS  :  CHANSCAN='+chanscan
print,' '
;
;  Is CHANSCAN an acceptable string ? If not, Return with Error.
;
if (chan eq 'RH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FDS_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'RL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FDS_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LH') then begin
  if ((smode ne 'SS')and(smode ne 'LF')) then begin
    print,'FDS_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FDS_ERRORS : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else begin
  print,'FDS_ERRORS : Invalid Channel'
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
cal_glitch=0 & cal_nifgs=0 & cal_sp=0 & cal_tm=0 & cal_wgts=0 & dihd=0
glat=0 & glon=0 & ical=0 & nifgs=0 & refh=0 & skyh=0 & sky_glitch=0
sky_wgts=0 & sp=0 & tm=0
;
cal_wtds_ds=0 & corr_index=0 & csp=0 & c_calsp=0 & del_temp=0 & f=0
index=0 & ndf_ejv=0 & ndf_none=0 & cal_wgts_ds = 0 & sqr_mat=0
step_dn=0 & step_up=0 & vib=0 & xcal=0
;
; Re-Define UnNeeded Mode-Specific Restored Parameters.
;
if(chanscan eq 'RHSS')then begin
 b_rhs = 1 & d_rhs=1 & l_rhs=1 & n_rhs=1 & urhs=1
 crhs=1 & ejv_rhs=1 & s_rhs=1 & tau_rhs=1
endif
if(chanscan eq 'RHLF')then begin
 b_rhf = 1 & d_rhf=1 & d_rhsf=1 & d_rhlf=1 & l_rhf=1 & n_rhf=1 & urhf=1
 cal_lbl = 1 & sky_lbl=1 & sl_weight_ratio=1
 crhf=1 & ejv_rhf=1 & s_rhf=1 & tau_rhf=1
endif
if(chanscan eq 'RLSS')then begin
 b_rls = 1 & d_rls=1 & l_rls=1 & n_rls=1 & urls=1
 crls=1 & ejv_rls=1 & s_rls=1 & tau_rls=1
endif
if(chanscan eq 'RLSF')then begin
 b_rsf = 1 & d_rsf=1 & d_rlsf=1 & d_rlfl=1 & l_rsf=1 & n_rsf=1 & ursf=1
 cal_lbl = 1 & sky_lbl=1 & sl_weight_ratio=1
 crsf=1 & ejv_rsf=1 & s_rsf=1 & tau_rsf=1
endif
if(chanscan eq 'RLLF')then begin
 b_rlf = 1 & d_rlf=1 & l_rlf=1 & n_rlf=1 & urlf=1
 crlf=1 & ejv_rlf=1 & & s_rlf=1 & tau_rlf=1
endif
if(chanscan eq 'LHSS')then begin
 b_lhs = 1 & d_lhs=1 & l_lhs=1 & n_lhs=1 & ulhs=1
 clhs=1 & ejv_lhs=1 & s_lhs=1 & tau_lhs=1
endif
if(chanscan eq 'LHLF')then begin
 b_lhf = 1 & d_lhf=1 & d_rhsf=1 & d_rhlf=1 & l_lhf=1 & n_lhf=1 & ulhf=1
 cal_lbl = 1 & sky_lbl=1 & sl_weight_ratio=1
 clhf=1 & ejv_lhf=1 & s_lhf=1 & tau_lhf=1
endif
if(chanscan eq 'LLSS')then begin
 b_lls = 1 & d_lls=1 & l_lls=1 & n_lls=1 & ulls=1
 clls=1 & ejv_lls=1 & s_lls=1 & tau_lls=1
endif
if(chanscan eq 'LLSF')then begin
 b_lsf = 1 & d_lsf=1 & d_llsf=1 & d_llfl=1 & l_lsf=1 & n_lsf=1 & ulsf=1
 cal_lbl = 1 & sky_lbl=1 & sl_weight_ratio=1
 clsf=1 & ejv_lsf=1 & s_lsf=1 & tau_lsf=1
endif
if(chanscan eq 'LLLF')then begin
 b_llf = 1 & d_llf=1 & l_llf=1 & n_llf=1 & ullf=1
 cllf=1 & ejv_llf=1 & s_llf=1 & tau_llf=1
endif
;
;  Compute the auxiliary weights and covariance matrices
;
print,'Computing the auxiliary matrices.'
print,' '
fds_errmat,pcvr=pcvr,rect=rect,diag=diag,beta=beta,omega=omega,$
    strp_conv=strp_conv
;
;
;  Define and fill in the FIP_DSQ record.
;
print,' '
print,'Filling in the FIP_DSQ record.'
print,' '
dsq_rec = fds_fip_dsq(1)
;
;  Set the header information
;
dsq_rec.model_label = string(solution,format='(A40)')
dsq_rec.chanscan = string(chanscan,format='(A4)')
;
;  Fill in the Stripe Covariance Matrix
;
qsize = size(pcvr)
nj = qsize(1)-1 & nk = qsize(2)-1
for j=0,nj do begin
   for k=0,nk do begin
      dsq_rec.q_inverse(j,k) = pcvr(j,k)
      dsq_rec.stripe_convert(j,k) = strp_conv(j,k)
   endfor
endfor
;
;  Write the FIP_DSQ file.
;
print,' '
print,'Writing the FIP_DSQ file.'
print,' '
fext = solution
i = strpos(fext,' ')
s1 = strtrim(strmid(fext,i,35),2)
fext = strtrim(s1 + '_' + strmid(fext,0,i),2)
outfile1 = 'csdr$firas_out:fip_dsq_' + chanscan + '.' + fext
rec_len = 444                 ; fixed length record size in bytes
openw,1,outfile1, rec_len, /fixed
writeu,1,dsq_rec
close,1
;
print,'File "'+strupcase(outfile1)+'" Written.'
;
;  Define and fill in the FIP_DSE records.
;
print,' '
print,'Filling in the FIP_DSE records.'
print,' '
;
;  Initialze the data records.
;
rsize = size(rect)
npix = rsize(1) & nstr = rsize(2)
dse_rec = fds_fip_dse(npix)
;
;  Get the pixel numbers and the pixel-ordered attitude information.
;
hist = histogram(min=0,px)
pixel = where(hist gt 0)
attitude = coorconv(pixel,infmt='P',outfmt='L',inco='F',outco='G')
;
;  Loop over the pixels.
;
for j=0,npix-1 do begin
   dse_rec(j).pixel = pixel(j)
   dse_rec(j).model_label = string(solution,format='(A40)')
   dse_rec(j).chanscan = string(chanscan,format='(A4)')
   dse_rec(j).diag_element = diag(j)
   for k=0,nstr-1 do begin
      dse_rec(j).rect_element(k)  = rect(j,k)
      dse_rec(j).beta_element(k)  = beta(j,k)
   endfor
   dse_rec(j).stripe_contrib = 0
   if (abs(attitude(j,1)) ge galcut) then begin
      nx=where(px eq pixel(j))
      wtot=total(sky_wgts_ds(nx))
      if(wtot gt 0.)then dse_rec(j).stripe_contrib = 1
   endif
endfor
;
;  Write the FIP_DSE file.
;
print,' '
print,'Writing the FIP_DSE file.'
print,' '
outfile2 = 'csdr$firas_out:fip_dse_' + chanscan + '.' + fext
rec_len = 140                 ; fixed length record size in bytes
openw,2,outfile2, rec_len, /fixed
for j=0,npix-1 do writeu,2,dse_rec(j)
close,2
;
print,'File "'+strupcase(outfile2)+'" Written.'
print,' '
;
; Set Error Status to NO Error.
;
error=0
;
return
end
