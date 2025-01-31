Pro FTT_LOSL,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_LOSL,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:lls.iss'
restore,'csdr$firas_in:rls.iss'
restore,'csdr$firas_in:lls_destriped.iss'
restore,'csdr$firas_in:rls_destriped.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / clls(*,0)^2
r_wgt = 1 / crls(*,0)^2
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
c_losl = (l_wgt * l_wgt * clls(*,0)^2) + (r_wgt * r_wgt * crls(*,0)^2)
d_losl = (l_wgt * l_wgt * d_lls(*,0)^2) + (r_wgt * r_wgt * d_rls(*,0)^2)
;

; Combined Num_IFGs
; -----------------
n_losl  = tl*n_lls + tr*n_rls
;

; Combined Galactic Latitude
; --------------------------
b_losl = (tl*n_lls*b_lls + tr*n_rls*b_rls) / (n_losl + (n_losl eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read RLSS Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_rlss.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_r = fex_struct.off_spec(4:37)
;

; Read LLSS Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_llss.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_l = fex_struct.off_spec(4:37)
;

; Combined Spectra
; ----------------
s_losl = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_losl(*,*,i)=(l_wgt(i)*(s_lls(*,*,i)-off_l(i))*n_lls + $
                r_wgt(i)*(s_rls(*,*,i)-off_r(i))*n_rls) / $
               (l_wgt(i)*n_lls + r_wgt(i)*n_rls + (n_losl eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:losl.iss'
save,file=sname,s_losl,n_losl,b_losl,c_losl,d_losl
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
l_rls=0 & urls=0 & tau_rls=0 & ejv_rls=0
l_lls=0 & ulls=0 & tau_lls=0 & ejv_lls=0
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
