;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FRD_MERGE_CSP computes and writes out the FIRAS merge correction spectra.
;
;DESCRIPTION:
;    An IDL program to compute the FIRAS merge correction spectra for use by
;    FMS in performing first, second, and third order skymap merges.  The program
;    reads the FEF difference and ratio spectra and the relative weights of
;    calibration model solutions (the C-Vectors), from which it computes the
;    correction spectra that are required by FMS.  The correction spectra are
;    written out to FEX_MCS_CCSS.F16_93HYBRID binary files for use by FMS and
;    are saved to the FEX_MCS_CCSS.ISS savesets for use by the tennis tree
;    programs.
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS
;     SMODE (I) : String defining channel and scan mode
;                 Allowed values are : 'HIGH' , 'LOWF' , 'HRES' ,
;                                      'HIG2' , 'LOW2' , 'LRES' .
;     ERROR (O) : Error status for procedure
;
;WARNINGS:
;     The IDL path must be set before invoking this procedure.
;
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_REF   directory containing the FEX_CVS_CCSS datasets
;	CSDR$FIRAS_IN    directory containing the FMS_CVS_CCSS datasets
;	CSDR$FIRAS_SAVE	 directory containing the FEF savesets
;	CSDR$FIRAS_OUT	 directory where the correction spectra
;                        will be written to the binary files and saved to the
;                        savesets
;
;EXAMPLE:
;     $ UIDL
;     !path = 'source directory,' + !path
;     setlog,'csdr$firas_ref','w'
;     setlog,'csdr$firas_in','x'
;     setlog,'csdr$firas_save','y'
;     setlog,'csdr$firas_out','z'
;     smode='LOWF'
;     frd_merge_csp,smode,error
;     if(error ne 0)then print,'Error Returned from FRD_Merge_CSP !'
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Calls IDL function FRD_CAL_WEIGHTS to get the calibration model solution
;    weights.
;    Calls IDL procedure FRD_COMPUTE_CSP to compute the correction spectra for
;    the RHSS, RHFA, RLSS, RLFA, LHSS, LHFA, LLSS, and LLFA skymaps.
;    Calls IDL procedure FRD_COMPUTE_HR to compute the correction spectra for
;    the RLLF, LLLF, LOSL, LOFA, HISL, and HIFA skymaps.
;    Calls IDL procedure FRD_COMPUTE_LR to compute the correction spectra for
;    the HIGH and LOWF skymaps.
;    Calls IDL function FRD_FEX_MCS to define the FEX_MCS record structure
;
;MODIFICATION HISTORY:
;     Written by Gene Eplee, General Sciences Corp., 28 October 1994
;     Modified by Ken Jensen, Hughes STX, 2 November 1994,  Arrays S31_D
;         and S31_R  renamed to S30_D and S30_R  for clarity.
;     Modified by Gene Eplee, General Science Corp., 4 November 1994
;         to separate HRES and LRES computations.  FRD_COMPUTE_RES.PRO is
;         renamed FRD_COMPUTE_HRES.PRO to do HRES computation and a new
;         FRD_COMPUTE_LRES.PRO will to LRES computation.
;     Modified by Ken Jensen, Hughes STX, 7 Nov 1994 to call FRD_COMPUTE_HR
;         instead of FRD_COMPUTE_HRES and FRD_COMPUTE_LR instead of
;         FRD_COMPUTE_LRES .
;     Modified by Ken Jensen, Hughes STX, 15 Nov 1994 to correctly compute
;     the correction spectra for the third order merge.
;     Modified by Ken Jensen, Hughes STX, 1 Dec 1994 to add the 'HIG2' and
;     'LOW2' options for calculating LOSL, LOFA, HISL and HIFA corrections.
;     Modified by Ken Jensen, Hughes STX, 9 Dec 1994 to write the
;     Gain CORRECTION spectra rather than the GAIN spectra.
;
;-
;______________________________________________________________________________
;
Pro FRD_Merge_CSP,smode,error
;
; Set Error Status
;
error=1
;
if N_Params() ne 2 then begin
 print,' '
 print,'FRD_Merge_CSP : Called Incorrectly : FRD_Merge_CSP,smode,error'
 print,' '
 return
endif
;
smode=strupcase(smode)
print,' '
print,'FRD_Merge_CSP  :  SMODE='+smode
print,' '
;
;  Is SMODE an acceptable string ? If not, Return with Error.
;
if ((smode ne 'LOWF')and(smode ne 'HIGH')and(smode ne 'HIG2')$
    and(smode ne 'HRES')and(smode ne 'LRES')and(smode ne 'LOW2')) then begin
    print,'FRD_Merge_Csp : ScanMode'
    print,' '
    return
endif
;
; Print Logical Translations.
;
print,'Logical Translations :'
ret=trnlog('csdr$firas_ref',intrans,/full,/issue_error)
print,'CSDR$FIRAS_REF  == '+strupcase(intrans)
ret=trnlog('csdr$firas_in',intrans,/full,/issue_error)
print,'CSDR$FIRAS_IN   == '+strupcase(intrans)
ret=trnlog('csdr$firas_save',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_SAVE == '+strupcase(outtrans)
ret=trnlog('csdr$firas_out',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_OUT  == '+strupcase(outtrans)
print,' '
;
; Define the file names.
;
if (smode eq 'LOWF') then begin
    infile1 = 'csdr$firas_ref:fex_cvs_llss.f16_93hybrid'
    infile2 = 'csdr$firas_ref:fex_cvs_rlss.f16_93hybrid'
    infile3 = 'csdr$firas_in:fms_cvs_rlfa.f16_93hybrid'
    infile4 = 'csdr$firas_in:fms_cvs_llfa.f16_93hybrid'
    savefile1 = 'csdr$firas_save:rls_fef_5.iss'
    savefile2 = 'csdr$firas_save:rls_fef_6.iss'
    savefile3 = 'csdr$firas_save:rsf_fef_5.iss'
    savefile4 = 'csdr$firas_save:lls_fef_6.iss'
    nrecs = 4
    chsn = ['LLSS','RLSS','RLFA','LLFA']
endif else if (smode eq 'HIGH') then begin
    infile1 = 'csdr$firas_ref:fex_cvs_lhss.f16_93hybrid'
    infile2 = 'csdr$firas_ref:fex_cvs_rhss.f16_93hybrid'
    infile3 = 'csdr$firas_in:fms_cvs_rhfa.f16_93hybrid'
    infile4 = 'csdr$firas_in:fms_cvs_lhfa.f16_93hybrid'
    savefile1 = 'csdr$firas_save:rhs_fef_5.iss'
    savefile2 = 'csdr$firas_save:rhs_fef_6.iss'
    savefile3 = 'csdr$firas_save:rhf_fef_5.iss'
    savefile4 = 'csdr$firas_save:lhs_fef_6.iss'
    nrecs = 4
    chsn = ['LHSS','RHSS','RHFA','LHFA']
endif else if (smode eq 'HRES') then begin
    infile1 = 'csdr$firas_ref:fex_cvs_lllf.f16_93hybrid'
    infile2 = 'csdr$firas_ref:fex_cvs_rllf.f16_93hybrid'
    savefile1 = 'csdr$firas_save:rlf_fef_5.iss'
    nrecs = 2
    chsn = ['LLLF','RLLF']
endif else if (smode eq 'HIG2') then begin
    infile1 = 'csdr$firas_in:fms_cvs_hisl.f16_93hybrid'
    infile2 = 'csdr$firas_in:fms_cvs_hifa.f16_93hybrid'
    savefile1 = 'csdr$firas_save:fef_8.iss'
    nrecs = 2
    chsn = ['HISL','HIFA']
endif else if (smode eq 'LOW2') then begin
    infile1 = 'csdr$firas_in:fms_cvs_losl.f16_93hybrid'
    infile2 = 'csdr$firas_in:fms_cvs_lofa.f16_93hybrid'
    savefile1 = 'csdr$firas_save:fef_9.iss'
    nrecs = 2
    chsn = ['LOSL','LOFA']
endif else if (smode eq 'LRES') then begin
    savefile1 = 'csdr$firas_save:fef_7.iss'
    nrecs = 2
    chsn = ['HIGH','LOWF']
endif
;
;  Get the Cal Weights.
;
print,' '
print,'Getting the Calibration Model Solution Weights.'
print,' '
if ((smode eq 'LOWF')or(smode eq 'HIGH')) then begin
    cw0 = frd_cal_weights(infile1,smode)
    cw1 = frd_cal_weights(infile2,smode)
    cw2 = frd_cal_weights(infile3,smode)
    cw3 = frd_cal_weights(infile4,smode)
endif else if ((smode eq 'HRES')or(smode eq 'HIG2')or(smode eq 'LOW2')) then begin
    cw0 = frd_cal_weights(infile1,smode)
    cw1 = frd_cal_weights(infile2,smode)
endif
;
;  Get the difference and ratio spectra.
;  Define away the unneeded constants.
;
print,' '
print,'Getting the difference and ratio spectra.'
print,' '
if ((smode eq 'LOWF')or(smode eq 'HIGH')) then begin
    restore,savefile1
        s01_d = -diff
        s01_r = 1.0/ratio
    restore,savefile2
        s12_d = diff
        s12_r = ratio
    restore,savefile3
        s23_d = diff
        s23_r = ratio
    restore,savefile4
        s30_d = -diff
        s30_r = 1.0/ratio
endif else if (smode eq 'HRES') then begin
    restore,savefile1
        s01_d = -diff
        s01_r = 1.0/ratio
endif else if (smode eq 'HIG2') then begin
    restore,savefile1
        s01_d = diff
        s01_r = ratio
endif else if (smode eq 'LOW2') then begin
    restore,savefile1
        s01_d = diff
        s01_r = float(diff)*0.+1.
endif else if (smode eq 'LRES') then begin
    restore,savefile1
        s01_d = diff
        s01_r = float(diff)*0.
endif
d_unc=0 & id1=0 & id2=0 & r_unc=0 & step_dn=0 & step_up=0 & tau=0
;
;  Compute the correction spectra
;
print,' '
print,'Computing the correction spectra.'
print,' '
if ((smode eq 'LOWF')or(smode eq 'HIGH')) then begin
    frd_compute_csp,s01_d,s12_d,s23_d,s30_d, $
                    s01_r,s12_r,s23_r,s30_r, $
                    cw0,cw1,cw2,cw3, $
                    off_cor,gain_cor
endif else if (smode eq 'HRES') then begin
    frd_compute_hr,s01_d,s01_r, $
                   cw0,cw1, $
                   off_cor,gain_cor
endif else if (smode eq 'HIG2') then begin
    frd_compute_hr,s01_d,s01_r, $
                   cw0,cw1, $
                   off_cor,gain_cor
    f = (4 + findgen(167)) * 144.981 / 256.
endif else if (smode eq 'LOW2') then begin
    frd_compute_hr,s01_d,s01_r, $
                   cw0,cw1, $
                   off_cor,gain_cor
    f = (4 + findgen(34)) * 144.981 / 256.
endif else if (smode eq 'LRES') then begin
    frd_compute_lr,s01_d,s01_r, $
                   off_cor,gain_cor
    f = (4 + findgen(34)) * 144.981 / 256.
endif
;

;  Returned GAIN_COR is actually the GAIN, NOT the correction.
;  The correction is the INVERSE of the gain
;
gain_cor = 1. / gain_cor
;

;  Save the correction spectra to the IDL savesets.
;
print,' '
print,'Saving the correction spectra.'
print,' '
save,file='csdr$firas_out:fex_mcs_'+smode+'.iss',f,off_cor,gain_cor,chsn
;
;
;  Define and fill in the FEX_MCS records.
;
print,' '
print,'Filling in the FEX_MCS records.'
print,' '
fex_recs = frd_fex_mcs(nrecs)
;
if ((smode eq 'LOWF')or(smode eq 'LOW2')or(smode eq 'LRES')) then begin
   jstart = 4
   jstop = 37
endif else if ((smode eq 'HIGH')or(smode eq 'HIG2')) then begin
   jstart = 4
   jstop = 170
endif else if (smode eq 'HRES') then begin
   jstart = 8
   jstop = 155
endif
;
for i=0,nrecs-1 do begin
    fex_recs(i).chanscan = STRING(chsn(i),FORMAT='(A4)')
    for j=jstart,jstop do fex_recs(i).offset_spec(j) = off_cor(j-jstart,i)
    for j=jstart,jstop do fex_recs(i).gain_spec(j) = gain_cor(j-jstart,i)
endfor

;
;  Write the FEX_MCS files.
;
print,' '
print,'Writing the FEX_MCS files.'
print,' '
rec_len = 3088              ;  fixed length record size in bytes
for j=0,nrecs-1 do begin
    lun=j+10
    outfile = 'csdr$firas_out:fex_mcs_' + chsn(j) + '.f16_93hybrid'
    openw,lun,outfile,rec_len,/fixed
    writeu,lun,fex_recs(j)
    close,lun
endfor
;
; Set Error Status to NO Error.
;
error=0
;
return
end
