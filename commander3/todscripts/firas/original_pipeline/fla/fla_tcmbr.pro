;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FLA_TCMBR fits a dipole to the FIRAS CMBR temperature map.
;
;DESCRIPTION:
;     An IDL program to extract the FIRAS CMBR temperature map from the line
;     emission IDL savesets, to fit a dipole to the temperature map, and to 
;     write the temperature list and the residual temperature list to the IDL
;     saveset FLA_TCMBR_CCSS.ISS and to the FIRAS project dataset binary file
;     FIP_TCB_CCSS.ISS.  The results of the fit are written to the text file
;     FLA_TCMBR_CCSS.TXT for use in constructing FITS headers.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS
;     CHANSCAN   (I)  :  String defining channel and scan mode (e.g. 'RHSS')
;     ERROR      (O)  :  Error status for procedure
;
;WARNINGS:
;     The IDL path must be set before invoking this procedure.
;     This procedure will only run on low frequency channels and scan modes.
;
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_SAVE  directory containing line emission IDL savesets
;                        and where the temperature saveset will be written
;	CSDR$FIRAS_OUT	 directory where the temperature project dataset
;                        will be written
;
;EXAMPLE:
;     $ UIDL
;     !path = 'source directory,' + !path
;     setlog,'csdr$firas_save','x'
;     setlog,'csdr$firas_out','y'
;     chanscan='LLSS'
;     FLA_TCMBR_G,chanscan,error
;     if(error ne 0)then print,'Error Returned from FLA_TCMBR !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL procedure FLA_DIPOLE to perform the dipole fit and compute the
;     residual temperature map.
;     Calls IDL function FLA_FIP_TCB to define the FIP_TCB record structure.
;     Calls IDL function COORCONV to perform coordinate conversions.
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 13 October 1994
;     Modified by Ken Jensen, Hughes STX, 24 October 1994, CSDR$FIRAS_IN 
;      replaced by CSDR$FIRAS_SAVE.
;     Modified by Joel Gales, Applied Research Corp., 10 January 1995
;       Change GALATEXC from 20 degrees to 15 degrees.
;
;-
;______________________________________________________________________________
;
Pro FLA_TCMBR,chanscan,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FLA_TCMBR : Called Incorrectly : FLA_TCMBR,chanscan,error'
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
print,'FLA_TCMBR  :  CHANSCAN='+chanscan
print,' '
;
;  Is CHANSCAN an acceptable string ? If not, Return with Error.
;
if (chan eq 'RH') then begin
    print,'FLA_TCMBR : Invalid Channel/ScanMode'
    print,' '
    return
endif else if (chan eq 'RL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FLA_TCMBR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else if (chan eq 'LH') then begin
    print,'FLA_TCMBR : Invalid Channel/ScanMode'
    print,' '
    return
endif else if (chan eq 'LL') then begin
  if ((smode ne 'SS')and(smode ne 'SF')and(smode ne 'LF')) then begin
    print,'FLA_TCMBR : Invalid Channel/ScanMode'
    print,' '
    return
  endif
endif else begin
  print,'FLA_TCMBR : Invalid Channel'
  print,' '
  return
endelse
;
; Print Logical Translations.
;
print,'Logical Translations :'
ret=trnlog('csdr$firas_out',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
ret=trnlog('csdr$firas_save',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_SAVE == '+strupcase(outtrans)
print,' '
;

if (chanscan eq 'RLSS') then infile = 'csdr$firas_save:lines_rls.iss'
;
if (chanscan eq 'RLSF') then infile = 'csdr$firas_save:lines_rsf.iss'
;
if (chanscan eq 'RLLF') then infile = 'csdr$firas_save:lines_rlf.iss'
;
if (chanscan eq 'LLSS') then infile = 'csdr$firas_save:lines_lls.iss'
;
if (chanscan eq 'LLSF') then infile = 'csdr$firas_save:lines_lsf.iss'
;
if (chanscan eq 'LLLF') then infile = 'csdr$firas_save:lines_llf.iss'
;
;  Restore the savesets.
;
print,' '
print,'Restoring IDL Save Set '+strupcase(infile)
print,' '
restore,infile
;
; Re-Define UnNeeded Restored Parameters.
;
csqr=0 & f=0 & fit_func=0 & freq=0 & f_char=0 & gal_cut=0 & lfx=0
n_base=0 & pcvr=0 & mode=0 & dst_covar=0 & alpha=0
;
;  Extract the CMBR temperatures and temperature errors.
;    Compute the temperature weights from the errors.
temps = reform(dst(0,*))
temps_sig = reform(dst_sig(0,*))
temps_wgt = 1.0/temps_sig^2
tsize = size(temps)
npix = tsize(1)
nifgs = nifg_in_pix
;
;  Convert the galactic coordinates to pixel lists.
;     Convert the coordinates to ecliptic and equatorial.
;
glon = pix2dat(pixel=pix,raster=l_xxx)
glat = pix2dat(pixel=pix,raster=b_xxx)
galco = fltarr(npix,2)
galco(*,0) = glon
galco(*,1) = glat
eclco = coorconv(galco,infmt='l',inco='g',outfmt='l',outco='e') 
equco = coorconv(galco,infmt='l',inco='g',outfmt='l',outco='q') 
;
;  Perform the dipole fit.
;
print,' '
print,'Performing the dipole fit.'
print,' '
galatexc=15.0					; galactic cut for dipole fit
FLA_DIPOLE,temps,glon,glat,tcmbr,sig_tcmbr,$
           dip_amp,sig_dip_amp,dip_glon,dip_glat,sig_dip_dir,$
           weights=temps_wgt,galexc=galatexc,resid=resid_temp
rms_resid = sqrt(total(temps_wgt*resid_temp^2)/total(temps_wgt))
;
;  Save the results into an IDL saveset.
;
print,' '
print,'Saving the results to the TCB IDL saveset.'
print,' '
outfile1='csdr$firas_save:fla_tcmbr_' + chanscan + '.iss'
save,file=outfile1,$
     temps,temps_sig,glon,glat,nifgs,pix,galatexc,resid_temp,rms_resid,$
     tcmbr,sig_tcmbr,dip_amp,sig_dip_amp,dip_glon,dip_glat,sig_dip_dir
;
;  Write out the fits results for use in FITS headers.
;
print,' '
print,'Writing the results to the TCB text file.'
print,' '
outfile2 = 'csdr$firas_save:fla_tcmbr_' + chanscan + '.txt'
openw,2,outfile2
printf,2,chanscan + ' CMBR Temperature Dipole Fit Results'
printf,2,' '
printf,2,format='("T_cmbr         =   ",f9.7,"     K")',tcmbr
printf,2,format='("T_cmbr Sig     =  ",e14.7," K")',sig_tcmbr
printf,2,' '
printf,2,format='("Dipole Amp     =  ",e14.7," K")',dip_amp
printf,2,format='("Dipole Amp Sig =  ",e14.7," K")',sig_dip_amp
printf,2,' '
printf,2,format='("Dipole Lon     = ",f9.5,"       deg")',dip_glon
printf,2,format='("Dipole Lat     =  ",f9.6,"      deg")',dip_glat
printf,2,format='("Dipole Dir Sig =   ",f10.8,"    deg")',sig_dip_dir
printf,2,' '
printf,2,format='("RMS Temp Resid =  ",e14.7," K")',rms_resid
close,2
;

;  Define and fill in the FIP_TCB records.
;
print,' '
print,'Filling in the FIP_TCB records.'
print,' '
tcb_recs = fla_fip_tcb(npix)
;
;  Loop over the pixels.
;
FOR j=0,npix-1 DO BEGIN
   tcb_recs(j).pixel = pix(j)
   tcb_recs(j).eclon = eclco(j,0)
   tcb_recs(j).eclat = eclco(j,1)
   tcb_recs(j).temp = temps(j)
   tcb_recs(j).temp_sig = temps_sig(j)
   tcb_recs(j).resid_temp = resid_temp(j)
   tcb_recs(j).chanscan = string(chanscan,format='(A4)')
   tcb_recs(j).num_ifgs = nifgs(j)
   tcb_recs(j).galon = glon(j)
   tcb_recs(j).galat = glat(j)
   tcb_recs(j).ra = equco(j,0)
   tcb_recs(j).dec = equco(j,1)
   tcb_recs(j).galatexc = galatexc
ENDFOR
;
;  Write the FIP_TCB file.
;
print,' '
print,'Writing the FIP_TCB file.'
print,' '
outfile3 = 'csdr$firas_out:fip_tcb_' + chanscan + '.f16_93hybrid'
rec_len = 52                 ; fixed length record size in bytes
OPENW,3,outfile3, rec_len, /fixed
for j=0,npix-1 do writeu,3,tcb_recs(j)
CLOSE,3
;
print,' '
print,'File "'+strupcase(outfile3)+'" Written.'
print,' '
;
; Set Error Status to NO Error.
;
error=0
;
return
end
