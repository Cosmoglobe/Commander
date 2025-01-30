;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_OST_LOWF creates and saves the binary records of LOWF stripe arrays
;               in the orthogonal basis.
;
;  Written by :  Ken Jensen,  Hughes STX,  23-May-1997
;-
;______________________________________________________________________________
;
Pro FMD_OST_LOWF,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_OST_LOWF : Called Incorrectly : FMD_OST_LOWF,error'
 print,' '
 return
endif
;

; Logical Translations
; --------------------
PRINT,'Logical Translations :'
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
PRINT,' '
PRINT,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
PRINT,' '
;

; Restore the savesets
; --------------------
PRINT,'Restoring LOWF_EJG.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:lowf_ejg.iss'
;
PRINT,'Restoring LOWF_ERRORS.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:lowf_errors.iss'
;

; Define the FMD_OST records
; --------------------------
sz = SIZE(beta)
nrec = sz(2)
ost_rec = FMD_OST_ST(nrec)
;

; Fill in the FMD_OST records
; ---------------------------
PRINT,' '
PRINT,'Filling in the FMD_OST records.'
PRINT,' '
sz = SIZE(gamma_lowf)
nf = sz(1)
FOR i=0,nrec-1 DO BEGIN
 ost_rec(i).str_gama(0:nf-1) = gamma_lowf(*,i)
 ost_rec(i).str_beta(cmp_px) = FLOAT(beta(*,i))
ENDFOR
;

; Write the FMD_OST file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_OST file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_ost_lowf.pass4'
rec_len = 25304             ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,ost_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_OST_LOWF.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
sky_wgts_ds=0 & chi2_loch=0 & pcvr=0 & rect=0 & diag=0 & omega=0 & square=0
d_inv=0 & stripe_conv=0 & destriper_wgt=0 & stripe_contrib=0 & lmat=0
cal_wgts_ds=0 & cmn_bol=0 & cmn_dihd=0 & cmn_tm=0 & cvec_cut=0
cvec_mask=0 & chanscan=0 & chan_label=0 & stripes_descrip=0 & chi2hich=0
del_temp=0 & dihed_cut=0 & dihed_cut_min=0 & dihed_pow=0
dirbe_array=0 & dirbe_cut=0 & dirbe_mask=0 & freq_band=0 & freq_corr=0
good_cal_dihd=0 & good_sky_dihd=0 & latcut=0 & loncut=0
lp_order=0 & max_frac=0 & n_cmn=0 & n_krnl=0 & n_stripes=0
ref_dihd=0 & ref_s0=0 & rms_s0=0 & sky_mask=0 & step_up=0 & step_dn=0
stripe_id=0 & stripe_order=0 & tmin=0 & tmax=0 & xtm=0
f_lo=0 & c_lowf=0 & ejg_lowf=0 & ndf_ejg_lowf=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
