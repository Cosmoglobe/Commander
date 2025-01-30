Pro FTT_HIFA,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_HIFA,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:lhf.iss'
restore,'csdr$firas_in:rhf.iss'
restore,'csdr$firas_in:lhf_destriped.iss'
restore,'csdr$firas_in:rhf_destriped.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / clhf(*,0)^2
r_wgt = 1 / crhf(*,0)^2
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
c_hifa = (l_wgt * l_wgt * clhf(*,0)^2) + (r_wgt * r_wgt * crhf(*,0)^2)
d_hifa = (l_wgt * l_wgt * d_lhf(*,0)^2) + (r_wgt * r_wgt * d_rhf(*,0)^2)
;

; Combined Num_IFGs
; -----------------
n_hifa  = tl*n_lhf + tr*n_rhf
;

; Combined Galactic Latitude
; --------------------------
b_hifa = (tl*n_lhf*b_lhf + tr*n_rhf*b_rhf) / (n_hifa + (n_hifa eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read RHFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_rhfa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_r = fex_struct.off_spec(4:170)
gain_r = fex_struct.gain_spec(4:170)
;
off_r(34:166)=off_r(34:166)*0. & gain_r(0:33)=1.
;

; Read LHFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_lhfa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_l = fex_struct.off_spec(4:170)
gain_l = fex_struct.gain_spec(4:170)
;
off_l(34:166)=off_l(34:166)*0. & gain_l(0:33)=1.
;

; Combined Spectra
; ----------------
s_hifa = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_hifa(*,*,i)=(l_wgt(i)*(s_lhf(*,*,i)*gain_l(i)-off_l(i))*n_lhf + $
                r_wgt(i)*(s_rhf(*,*,i)*gain_r(i)-off_r(i))*n_rhf) / $
               (l_wgt(i)*n_lhf + r_wgt(i)*n_rhf + (n_hifa eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:hifa.iss'
save,file=sname,s_hifa,n_hifa,b_hifa,c_hifa,d_hifa
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
l_rhf=0 & urhf=0 & tau_rhf=0 & ejv_rhf=0 & d_rhsf=0 & d_rhlf=0
l_lhf=0 & ulhf=0 & tau_lhf=0 & ejv_lhf=0 & d_lhsf=0 & d_lhlf=0
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
