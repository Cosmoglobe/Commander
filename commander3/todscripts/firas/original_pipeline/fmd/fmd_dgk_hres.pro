;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FMD_DGK_HRES creates and saves the binary records of HRES coadd DIRBE
;               gradient kernals.
;
;  Written by :  Ken Jensen,  Hughes STX,  03-Jun-1997
;-
;______________________________________________________________________________
;
Pro FMD_DGK_HRES,error
;

; Set Error Status
; ----------------
error=1
;

; Correct Invocation ?
; --------------------
if N_Params() ne 1 then begin
 print,' '
 print,'FMD_DGK_HRES : Called Incorrectly : FMD_DGK_HRES,error'
 print,' '
 return
endif
;

; Logical Translations
; --------------------
PRINT,'Logical Translations :'
ret = TRNLOG('csdr$firas_ref',reftrans,/full,/issue_error)
ret = TRNLOG('csdr$firas_in',intrans,/full,/issue_error)
ret = TRNLOG('csdr$firas_out',outtrans,/full,/issue_error)
PRINT,' '
PRINT,'CSDR$FIRAS_REF == '+strupcase(reftrans)
PRINT,'CSDR$FIRAS_IN  == '+strupcase(intrans)
PRINT,'CSDR$FIRAS_OUT == '+strupcase(outtrans)
PRINT,' '
;

; Restore the savesets
; --------------------
PRINT,'Restoring HRES.ISS'
PRINT,' '
RESTORE,'csdr$firas_in:hres.iss'
;
PRINT,'Restoring HRES_EJG.ISS'
PRINT,' '
RESTORE,'csdr$firas_out:hres_ejg.iss'
;
PRINT,'Restoring FMD_DIRBE_FUNC_HRES.ISS'
PRINT,' '
RESTORE,'csdr$firas_ref:fmd_dirbe_func_hres.iss'
;

; Define the FMD_DGK records
; --------------------------
nrec = N_ELEMENTS(px)
dgk_rec = FMD_DGK_ST(nrec)
;

; Fill in the FMD_DGK records
; --------------------------
PRINT,' '
PRINT,'Filling in the FMD_DGK records.'
PRINT,' '
FOR i=0L,nrec-1 DO BEGIN
 dgk_rec(i).pixel = px(i)
 dgk_rec(i).chanscan = sky_lbl(i)
 dgk_rec(i).rec_num = fsl_idx(i)
 dgk_rec(i).weight = sky_wgts_ds(i)
 dgk_rec(i).band_8 = g8(i)
 dgk_rec(i).band_9 = g9(i)
 dgk_rec(i).band_10 = g10(i)
ENDFOR
;

; Write the FMD_DGK file
; ----------------------
PRINT,' '
PRINT,'Writing the FMD_DGK file.'
PRINT,' '
outfile = 'csdr$firas_out:fmd_dgk_hres.pass4'
rec_len = 28             ; fixed length record size in bytes
OPENW,1,outfile, rec_len, /fixed
FOR j=0L,nrec-1 DO WRITEU,1,dgk_rec(j)
CLOSE,1
;
fname = STRUPCASE(outtrans(0)) + ':FMD_DGK_HRES.PASS4'
PRINT,'File "' + fname + '" Written.'
PRINT,' '
;

; Re-Define Restored Fields
; -------------------------
cal_glitch=0 & cal_idx=0 & cal_nifgs=0 & cal_s0=0 & cal_tm=0
cal_wgts=0 & cal_wgts_ds=0 & cmn_bol=0 & cmn_dihd=0 & cmn_tm=0 & cvec_cut=0
cvec_mask=0 & chanscan=0 & chan_label=0 & stripes_descrip=0
del_temp=0 & dihd=0 & dihed_cut=0 & dihed_cut_min=0 & dihed_pow=0
dirbe_array=0 & dirbe_cut=0 & dirbe_mask=0 & freq_band=0 & freq_corr=0
glat=0 & glon=0 & good_cal_dihd=0 & good_sky_dihd=0 & ical=0 & latcut=0
loncut=0 & lp_order=0 & max_frac=0 & nifgs=0 & n_cmn=0 & n_krnl=0 & refh=0
n_stripes=0 & ref_dihd=0 & ref_s0=0 & rms_s0=0 & scan=0 & skyh=0 & sky_dihd=0
sky_glitch=0 & sky_mask=0 & sky_s0=0 & sky_wgts=0 & step_up=0 & step_dn=0
stripe_id=0 & stripe_order=0 & tm=0 & tmin=0 & tmax=0 & xcal=0 & xtm=0
sky_idx=0 & st_sub=0 & cal_lbl=0 & idd=0 & g8x=0 & g9x=0 & g10x=0
f8=0 & f9=0 & f10=0
f_hr=0 & c_hres=0 & ejg_hres=0 & ndf_ejg_hres=0 & gamma_hres=0 & chi2_hrch=0
;

; Set Error Status to NO Error
; ----------------------------
error=0
;

RETURN
END
