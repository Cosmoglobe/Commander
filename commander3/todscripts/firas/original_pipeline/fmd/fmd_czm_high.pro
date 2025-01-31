;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_CZM_HIGH creates and saves the binary records of HIGH coadd ZODI spectra.
;
;  Written by :  Ken Jensen,  Hughes STX,  18-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_CZM_HIGH,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_CZM_HIGH : Called Incorrectly : FMD_CZM_HIGH,error'
 print,' '
 return
endif
;

; Logical Translations
; --------------------
PRINT,'Logical Translations :'
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
PRINT,' '
PRINT,'CSDR$FIRAS_IN  == '+strupcase(intrans)
PRINT,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
PRINT,' '
;

; Restore the savesets
; --------------------
PRINT,'Restoring HIGH.ISS'
PRINT,' '
RESTORE,'csdr$firas_in:high.iss'
;
PRINT,'Restoring HIGH_2_EJG.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:high_2_ejg.iss'
;
PRINT,'Restoring HIGH_ZODISPEC.ISS'
PRINT,' '
RESTORE,'csdr$firas_in:high_zodispec.iss'
;

; Define the FMD_CZM records
; --------------------------
nrec = N_ELEMENTS(px)
czm_rec = FMD_CZM_ST(nrec)
;

; Fill in the FMD_CZM records
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_CZM records.'
PRINT,' '
FOR i=0L,nrec-1 DO BEGIN
 czm_rec(i).pixel = px(i)
 czm_rec(i).chanscan = sky_lbl(i)
 czm_rec(i).rec_num = fsl_idx(i)
 czm_rec(i).weight = sky_wgts_ds(i)
 czm_rec(i).zodi_mod = zodi_spec(40:209,i)
ENDFOR
;

; Write the FMD_CZM file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_CZM file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_czm_high.pass4'
rec_len = 696            ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0L,nrec-1 DO WRITEU,1,czm_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_CZM_HIGH.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
alpha=0 & cal_glitch=0 & cal_idx=0 & cal_nifgs=0 & cal_s0=0 & cal_tm=0
cal_wgts=0 & cal_wgts_ds=0 & cmn_bol=0 & cmn_dihd=0 & cmn_tm=0 & cvec_cut=0
cvec_mask=0 & chanscan=0 & chan_label=0 & stripes_descrip=0 & chi2hich=0
delbeta=0 & del_temp=0 & dihd=0 & dihed_cut=0 & dihed_cut_min=0 & dihed_pow=0
dirbe_array=0 & dirbe_cut=0 & dirbe_mask=0 & freq_band=0 & freq_corr=0
glat=0 & glon=0 & good_cal_dihd=0 & good_sky_dihd=0 & ical=0 & latcut=0
loncut=0 & lp_order=0 & max_frac=0 & nifgs=0 & n_cmn=0 & n_krnl=0 & refh=0
n_stripes=0 & ref_dihd=0 & ref_s0=0 & rms_s0=0 & scan=0 & skyh=0 & sky_dihd=0
sky_glitch=0 & sky_mask=0 & sky_s0=0 & sky_wgts=0 & step_up=0 & step_dn=0
stripe_id=0 & stripe_order=0 & tm=0 & tmin=0 & tmax=0 & wdelbeta=0 & xcal=0
xtm=0 & zeta=0 & zodimodpars=0 & zodi_freq=0 & zodi_model=0 & st_sub=0
sky_idx=0 & cal_lbl=0 & chi2_hi2ch=0
f_hi2=0 & c_high_2=0 & ejg_high_2=0 & ndf_ejg_high_2=0 & gamma_high_2=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
