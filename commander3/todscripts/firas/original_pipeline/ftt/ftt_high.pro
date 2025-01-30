Pro FTT_HIGH,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_HIGH,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:hisl.iss'
restore,'csdr$firas_in:hifa.iss'
;

; Frequency Dependent Weights
; ---------------------------
s_wgt = 1 / c_hisl
f_wgt = 1 / c_hifa
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
c_high = (s_wgt * s_wgt * c_hisl) + (f_wgt * f_wgt * c_hifa)
d_high = (s_wgt * s_wgt * d_hisl) + (f_wgt * f_wgt * d_hifa)
;

; Combined Num_IFGs
; -----------------
n_high  = ts*n_hisl + tf*n_hifa
;

; Combined Galactic Latitude
; --------------------------
b_high = (ts*n_hisl*b_hisl + tf*n_hifa*b_hifa) / (n_high + (n_high eq 0))
;


; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read HISL Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_hisl.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_s = fex_struct.off_spec(4:170)
gain_s = fex_struct.gain_spec(4:170)
;
off_s(34:166)=off_s(34:166)*0. & gain_s(0:33)=1.
;

; Read HIFA Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_hifa.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_f = fex_struct.off_spec(4:170)
gain_f = fex_struct.gain_spec(4:170)
;
off_f(34:166)=off_f(34:166)*0. & gain_f(0:33)=1.
;

; Combined Spectra
; ----------------
s_high = COMPLEXARR(96,64,nf)
FOR i=0,nf-1 DO $
 s_high(*,*,i)=(s_wgt(i)*(s_hisl(*,*,i)*gain_s(i)-off_s(i))*n_hisl + $
                f_wgt(i)*(s_hifa(*,*,i)*gain_f(i)-off_f(i))*n_hifa) / $
               (s_wgt(i)*n_hisl + f_wgt(i)*n_hifa + (n_high eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:high.iss'
save,file=sname,s_high,n_high,b_high,c_high,d_high
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
