Pro FTT_LOWF,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_LOWF,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:losl.iss'
restore,'csdr$firas_in:lofa.iss'
;

; Frequency Dependent Weights
; ---------------------------
s_wgt = 1 / c_losl
f_wgt = 1 / c_lofa
;
norm = s_wgt + f_wgt
s_wgt = s_wgt / norm
f_wgt = f_wgt / norm
;

; Frequency-Independent Weights
; -----------------------------
nf = N_ELEMENTS(s_wgt)
ts = TOTAL(s_wgt) / nf
tf = TOTAL(f_wgt) / nf
;

; Combined C and D Variances
; --------------------------
c_lowf = (s_wgt * s_wgt * c_losl) + (f_wgt * f_wgt * c_lofa)
d_lowf = (s_wgt * s_wgt * d_losl) + (f_wgt * f_wgt * d_lofa)
;

; Combined Num_IFGs
; -----------------
n_lowf  = ts*n_losl + tf*n_lofa
;

; Combined Galactic Latitude
; --------------------------
b_lowf = (ts*n_losl*b_losl + tf*n_lofa*b_lofa) / (n_lowf + (n_lowf eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read LOSL Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_losl.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_s = fex_struct.off_spec(4:37)
;

; Read LOFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_lofa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_f = fex_struct.off_spec(4:37)
;
; Combined Spectra
; ----------------
s_lowf = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_lowf(*,*,i)=(s_wgt(i)*(s_losl(*,*,i)-off_s(i))*n_losl + $
                f_wgt(i)*(s_lofa(*,*,i)-off_f(i))*n_lofa) / $
               (s_wgt(i)*n_losl + f_wgt(i)*n_lofa + (n_lowf eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:lowf.iss'
save,file=sname,s_lowf,n_lowf,b_lowf,c_lowf,d_lowf
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Set Status = No Error
; --------------------
error = 0
;

RETURN
END
