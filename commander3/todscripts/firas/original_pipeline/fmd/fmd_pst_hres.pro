;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_PST_HRES creates and saves the binary records of HRES stripe arrays
;               in the original (physical) basis.
;
;  Written by :  Ken Jensen,  Hughes STX,  09-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_PST_HRES,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_PST_HRES : Called Incorrectly : FMD_PST_HRES,error'
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
PRINT,'Restoring HRES_EJG.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:hres_ejg.iss'
;
PRINT,'Restoring HRES_ERRORS.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:hres_errors.iss'
;

; Stripe Description
; ------------------
descrip = ['DIRBE_BAND_10','DIRBE_BAND_9','TIME_1','TIME_2','TIME_3']
descrip = [descrip,'TIME_4','TIME_5','TIME_6','TIME_7']
descrip = $
 [descrip,'LLLF_DIHEDRAL','LLLF_BOLOMETER','LLLF_6K_HORN','LLLF_MISSION']
descrip = $
 [descrip,'RLLF_DIHEDRAL','RLLF_BOLOMETER','RLLF_6K_HORN','RLLF_MISSION']
;

; Define the FMD_PST records
; --------------------------
sz = SIZE(beta)
nrec = sz(2)
;
IF (nrec NE N_ELEMENTS(descrip)) THEN BEGIN
 PRINT,'FMD_PST_HRES : Error in Number of Records !'
 RETURN
ENDIF
;
pst_rec = FMD_PST_ST(nrec)
;

; Fill in the FMD_PST records
; ---------------------------
PRINT,' '
PRINT,'Filling in the FMD_PST records.'
PRINT,' '
sz = SIZE(gamma_hres)
nf = sz(1)
FOR i=0,nrec-1 DO BEGIN
 slen = STRLEN(descrip(i))
 IF(slen LT 40) THEN FOR j=slen,39 DO descrip(i) = descrip(i) + ' '
 pst_rec(i).str_id = descrip(i)
 pst_rec(i).str_spec(0:nf-1) = FLOAT(ejg_hres(*,i))
 pst_rec(i).str_rect(cmp_px) = FLOAT(rect(*,i))
 pst_rec(i).str_covr(0:nrec-1) = FLOAT(pcvr(i,*))
ENDFOR
;

; Write the FMD_PST file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_PST file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_pst_hres.pass4'
rec_len = 25460             ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0,nrec-1 DO WRITEU,1,pst_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_PST_HRES.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
sky_wgts_ds=0 & chi2_hrch=0 & diag=0 & omega=0 & square=0
d_inv=0 & stripe_conv=0 & destriper_wgt=0 & stripe_contrib=0 & lmat=0
cal_wgts_ds=0 & cmn_bol=0 & cmn_dihd=0 & cmn_tm=0 & cvec_cut=0
cvec_mask=0 & chanscan=0 & chan_label=0 & stripes_descrip=0
del_temp=0 & dihed_cut=0 & dihed_cut_min=0 & dihed_pow=0
dirbe_array=0 & dirbe_cut=0 & dirbe_mask=0 & freq_band=0 & freq_corr=0
good_cal_dihd=0 & good_sky_dihd=0 & latcut=0 & loncut=0
lp_order=0 & max_frac=0 & n_cmn=0 & n_krnl=0 & n_stripes=0
ref_dihd=0 & ref_s0=0 & rms_s0=0 & sky_mask=0 & step_up=0 & step_dn=0
stripe_id=0 & stripe_order=0 & tmin=0 & tmax=0 & xtm=0
f_hr=0 & c_hres=0 & ndf_ejg_hres=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
