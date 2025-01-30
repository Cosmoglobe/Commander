Pro FTT_LRES,error
;

; Set Status = Error
; ------------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FTT_LRES,error'
 return
endif
;

; Restore Skymap Save Sets
; ------------------------
restore,'csdr$firas_in:lowf.iss'
restore,'csdr$firas_in:high.iss'
;

; Frequency Dependent Weights
; ---------------------------
l_wgt = 1 / c_lowf
h_wgt = 1 / c_high
;
norm = l_wgt + h_wgt
l_wgt = l_wgt / norm
h_wgt = h_wgt / norm
;

; Frequency-Independent Weights
; -----------------------------
tl = TOTAL(l_wgt) / 34.
th = TOTAL(h_wgt(0:33)) / 34.
;

; Combined C and D Variances
; --------------------------
c_lres = (l_wgt * l_wgt * c_lowf) + (h_wgt * h_wgt * c_high)
d_lres = (l_wgt * l_wgt * d_lowf) + (h_wgt * h_wgt * d_high)
;

; Combined Num_IFGs
; -----------------
n_lres  = tl*n_lowf + th*n_high
;

; Combined Galactic Latitude
; --------------------------
b_lres = (tl*n_lowf*b_lowf + th*n_high*b_high) / (n_lres + (n_lres eq 0))
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; Read HIGH Correction Spectra
; ----------------------------
OPENR,10,'csdr$firas_ref:fex_mcs_high.f16_93hybrid',rec_len,/FIXED
READU, 10, fex_struct
CLOSE,10
;
off_h = fex_struct.off_spec(4:37)
;

; Combined Spectra
; ----------------
s_lres = s_high
FOR i=0,33 DO $
 s_lres(*,*,i)=(l_wgt(i)*s_lowf(*,*,i)*n_lowf + $
                h_wgt(i)*(s_high(*,*,i)-off_h(i))*n_high) / $
               (l_wgt(i)*n_lowf + h_wgt(i)*n_high + (n_lres eq 0))
;

; Create Merged Save Set
; ----------------------
sname='csdr$firas_out:lres.iss'
save,file=sname,s_lres,n_lres,b_lres,c_lres,d_lres
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
