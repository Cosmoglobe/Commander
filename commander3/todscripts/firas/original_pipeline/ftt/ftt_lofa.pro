Pro FTT_LOFA,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_LOFA,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:lsf.iss'
restore,'csdr$firas_in:rsf.iss'
restore,'csdr$firas_in:lsf_destriped.iss'
restore,'csdr$firas_in:rsf_destriped.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / clsf(*,0)^2
r_wgt = 1 / crsf(*,0)^2
;
norm = l_wgt + r_wgt
l_wgt = l_wgt / norm
r_wgt = r_wgt / norm
;

; Frequency-Independent Weights
; -----------------------------
nf = N_ELEMENTS(f)
tl = TOTAL(l_wgt) / nf
tr = TOTAL(r_wgt) / nf
;

; Combined C and D Variances
; --------------------------
c_lofa = (l_wgt * l_wgt * clsf(*,0)^2) + (r_wgt * r_wgt * crsf(*,0)^2)
d_lofa = (l_wgt * l_wgt * d_lsf(*,0)^2) + (r_wgt * r_wgt * d_rsf(*,0)^2)
;

; Combined Num_IFGs
; -----------------
n_lofa  = tl*n_lsf + tr*n_rsf
;

; Combined Galactic Latitude
; --------------------------
b_lofa = (tl*n_lsf*b_lsf + tr*n_rsf*b_rsf) / (n_lofa + (n_lofa eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read RLFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_rlfa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_r = fex_struct.off_spec(4:37)
;

; Read LLFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_llfa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_l = fex_struct.off_spec(4:37)
;

; Combined Spectra
; ----------------
s_lofa = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_lofa(*,*,i)=(l_wgt(i)*(s_lsf(*,*,i)-off_l(i))*n_lsf + $
                r_wgt(i)*(s_rsf(*,*,i)-off_r(i))*n_rsf) / $
               (l_wgt(i)*n_lsf + r_wgt(i)*n_rsf + (n_lofa eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:lofa.iss'
save,file=sname,s_lofa,n_lofa,b_lofa,c_lofa,d_lofa
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
px=0 & tm=0 & sp=0 & nifgs=0 & glon=0 & glat=0 & solution=0 & galcut=0
cal_sp=0 & cal_nifgs=0 & cal_tm=0 & xcal=0 & ical=0 & skyh=0 & refh=0 & dihd=0
sky_glitch=0 & cal_glitch=0 & sky_wgts=0 & cal_wgts=0 & sky_wgts_ds=0
cal_wgts_ds=0 & del_temp=0 & csp=0 & step_up=0 & step_dn=0 & ndf_ejv=0
sqr_mat=0 & rect=0 & diag=0 & c_calsp=0 & ndf_none=0 & pcvr=0 & vib=0
corr_index=0 & index=0 & sl_weight_rat=0 & cal_lbl=0 & sky_lbl=0
;
l_rsf=0 & ursf=0 & tau_rsf=0 & ejv_rsf=0 & d_rlsf=0 & d_rlfl=0
l_lsf=0 & ulsf=0 & tau_lsf=0 & ejv_lsf=0 & d_llsf=0 & d_llfl=0
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
