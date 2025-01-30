Pro FTT_HISL,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_HISL,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:lhs.iss'
restore,'csdr$firas_in:rhs.iss'
restore,'csdr$firas_in:lhs_destriped.iss'
restore,'csdr$firas_in:rhs_destriped.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / clhs(*,0)^2
r_wgt = 1 / crhs(*,0)^2
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
c_hisl = (l_wgt * l_wgt * clhs(*,0)^2) + (r_wgt * r_wgt * crhs(*,0)^2)
d_hisl = (l_wgt * l_wgt * d_lhs(*,0)^2) + (r_wgt * r_wgt * d_rhs(*,0)^2)
;

; Combined Num_IFGs
; -----------------
n_hisl  = tl*n_lhs + tr*n_rhs
;

; Combined Galactic Latitude
; --------------------------
b_hisl = (tl*n_lhs*b_lhs + tr*n_rhs*b_rhs) / (n_hisl + (n_hisl eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read RHSS Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_rhss.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_r = fex_struct.off_spec(4:170)
gain_r = fex_struct.gain_spec(4:170)
;
off_r(34:166)=off_r(34:166)*0. & gain_r(0:33)=1.
;

; Read LHSS Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_lhss.f16_93hybrid',rec_len,/FIXED
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
s_hisl = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_hisl(*,*,i)=(l_wgt(i)*(s_lhs(*,*,i)*gain_l(i)-off_l(i))*n_lhs + $
                r_wgt(i)*(s_rhs(*,*,i)*gain_r(i)-off_r(i))*n_rhs) / $
               (l_wgt(i)*n_lhs + r_wgt(i)*n_rhs + (n_hisl eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:hisl.iss'
save,file=sname,s_hisl,n_hisl,b_hisl,c_hisl,d_hisl
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
corr_index=0 & index=0
;
l_rhs=0 & urhs=0 & tau_rhs=0 & ejv_rhs=0
l_lhs=0 & ulhs=0 & tau_lhs=0 & ejv_lhs=0
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
