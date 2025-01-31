Pro FTT_HRES,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_HRES,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:llf.iss'
restore,'csdr$firas_in:rlf.iss'
restore,'csdr$firas_in:llf_destriped.iss'
restore,'csdr$firas_in:rlf_destriped.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / cllf(*,0)^2
r_wgt = 1 / crlf(*,0)^2
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
c_hres = (l_wgt * l_wgt * cllf(*,0)^2) + (r_wgt * r_wgt * crlf(*,0)^2)
d_hres = (l_wgt * l_wgt * d_llf(*,0)^2) + (r_wgt * r_wgt * d_rlf(*,0)^2)
;

; Combined Num_IFGs
; -----------------
n_hres  = tl*n_llf + tr*n_rlf
;

; Combined Galactic Latitude
; --------------------------
b_hres = (tl*n_llf*b_llf + tr*n_rlf*b_rlf) / (n_hres + (n_hres eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read RLLF Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_rllf.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_r = fex_struct.off_spec(8:155)
;

; Read LLLF Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_lllf.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_l = fex_struct.off_spec(8:155)
;

; Combined Spectra
; ----------------
s_hres = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_hres(*,*,i)=(l_wgt(i)*(s_llf(*,*,i)-off_l(i))*n_llf + $
                r_wgt(i)*(s_rlf(*,*,i)-off_r(i))*n_rlf) / $
               (l_wgt(i)*n_llf + r_wgt(i)*n_rlf + (n_hres eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:hres.iss'
save,file=sname,s_hres,n_hres,b_hres,c_hres,d_hres
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
l_rlf=0 & urlf=0 & tau_rlf=0 & ejv_rlf=0
l_llf=0 & ullf=0 & tau_llf=0 & ejv_llf=0
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
