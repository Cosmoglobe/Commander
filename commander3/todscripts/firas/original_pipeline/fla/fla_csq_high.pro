;______________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  FLA_CSQ_HIGH computes the combined chi-squared per pixel for the
;               HIGH merged skymap.
;
;DESCRIPTION:
;     An IDL procedure to compute the combined chi-squared per pixel for
;     the merge of HISL and HIFA skymaps, and to write the computed values
;     and associated data to the FIP_CSQ_HIGH.F16_93HYBRID binary file and
;     to IDL Save Set FIP_CSQ_HIGH.ISS
;
;CALLING SEQUENCE:
;     Invoked directly from IDL.
;
; ARGUMENTS
;     ERROR     (O)  :  Error status for procedure
;
;WARNINGS:
;     The following logical pointers must be defined before using this 
;     procedure:
;	CSDR$FIRAS_REF   directory containing FMS_CVS files
;	CSDR$FIRAS_IN    directory containing FMS_SKY files
;	CSDR$FIRAS_MCS   directory containing FEX_MCS files
;	CSDR$FIRAS_SAVE	 directory where FIP_CSQ IDL save set will be sent
;	CSDR$FIRAS_OUT	 directory where FIP_CSQ file will be written
;
;EXAMPLE:
;     $ UPRIDL
;     setlog,'csdr$firas_ref','cobearcv1:[firas_ref]'
;     setlog,'csdr$firas_in','fircoadd:[skyspec]'
;     setlog,'csdr$firas_mcs','cobearcv1:[firas_ref]'
;     setlog,'csdr$firas_save','fircoadd:[errors]'
;     setlog,'csdr$firas_out','fircoadd:[errors]'
;     FLA_CSQ_HIGH,error
;     if(error ne 0)then print,'Error Returned from FLA_CSQ_HIGH !'
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls IDL procedure FLA_MATCH to match the skymap indices.
;     Calls IDL function FLA_FEX_CVS_ST to define the FEX_CVS record structure.
;     Calls IDL function FLA_FMS_SKY_ST to define the FMS_SKY record structure.
;     Calls IDL function FLA_FIP_CSQ to define the FIP_CSQ record structure.
;
;MODIFICATION HISTORY:
;     Written by Ken Jensen, Hughes STX, 19 January 1995, SER ?????
;
;-
;______________________________________________________________________________
;
Pro FLA_CSQ_HIGH,error
;

; Return With Error if Incorrectly Invoked
; ----------------------------------------
if N_Params() ne 1 then begin
 print,'Procedure Invoked Incorrectly'
 print,'Correct Invocation is "FLA_CSQ_HIGH,error" '
 error = 1
 return
endif
;

; Error Status
; ------------
error = 1
;

; Print Logical Translations
; --------------------------
print,'Logical Translations :'
ret=trnlog('csdr$firas_ref',intrans,/full,/issue_error)
print,'CSDR$FIRAS_REF   == '+strupcase(intrans)
ret=trnlog('csdr$firas_in',intrans,/full,/issue_error)
print,'CSDR$FIRAS_IN    == '+strupcase(intrans)
ret=trnlog('csdr$firas_mcs',intrans,/full,/issue_error)
print,'CSDR$FIRAS_MCS   == '+strupcase(intrans)
ret=trnlog('csdr$firas_save',intrans,/full,/issue_error)
print,'CSDR$FIRAS_SAVE  == '+strupcase(intrans)
ret=trnlog('csdr$firas_out',outtrans,/full,/issue_error)
print,'CSDR$FIRAS_OUT   == '+strupcase(outtrans)
print,' '
;

; Skymap IDs
; ----------
id = 'HIGH'
map1 = 'HISL'  &  map2 = 'HIFA'
;

; Frequencies
; -----------
freq = (4.+findgen(167))*144.981/256.
fidx = indgen(160)
num_combined = 160
freq_combined = freq(fidx)
;

; Fetch C-Variances
; -----------------
openr,lunr,'csdr$firas_ref:fms_cvs_hisl.f16_93hybrid',/get_lun,/share
r=assoc(lunr,fla_fex_cvs_st(1))
cvs=r(0)
c1=cvs.cvector(4:163)
close,lunr & free_lun,lunr
;
openr,lunr,'csdr$firas_ref:fms_cvs_hifa.f16_93hybrid',/get_lun,/share
r=assoc(lunr,fla_fex_cvs_st(1))
cvs=r(0)
c2=cvs.cvector(4:163)
close,lunr & free_lun,lunr
;

; Frequency-Dependent Weights
; ---------------------------
w1=c2/(c1+c2)
w2=c1/(c1+c2)
;

; Define, Load and Store FEX_MCS dataset structure
; ------------------------------------------------
fex_struct = {CHANSCAN: STRING(' ',FORMAT='(A4)'), $
              OFF_SPEC: COMPLEXARR(257), $
              GAIN_SPEC: FLTARR(257)}
;
rec_len = 3088
;

; HISL Correction Spectra
; -----------------------
filename = 'fex_mcs_hisl.f16_93hybrid'
OPENR,10,'csdr$firas_mcs:' + filename,rec_len,/FIXED
READU, 10, fex_struct
fex1 = fex_struct
CLOSE,10
;

; HIFA Correction Spectra
; -----------------------
filename = 'fex_mcs_hifa.f16_93hybrid'
OPENR,10,'csdr$firas_mcs:' + filename,rec_len,/FIXED
READU, 10, fex_struct
fex2 = fex_struct
CLOSE,10
;

; Offset Corrections
; ------------------
o1=float(fex1.off_spec(4:163))   ; HISL Offset Correction
o2=float(fex2.off_spec(4:163))   ; HIFA Offset Correction
;

; Gain Corrections
; ----------------
g1=float(fex1.gain_spec(4:163))   ; HISL Gain Correction
g2=float(fex2.gain_spec(4:163))   ; HIFA Gain Correction
;

; Gain = 1 for Low Frequencies
; ----------------------------
g1(0:33)=1. & g2(0:33)=1.
;

; Offset = 0 for High Frequencies
; -------------------------------
o1(34:159)=0. &  o2(34:159)=0.
;

; Fetch HISL Skymap Records
; -------------------------
fname='csdr$firas_in:fms_sky_hisl.fad_8932800_9026409'
openr,lunr,fname,/get_lun,/share
finfo=fstat(lunr)
nrec=finfo.size/finfo.rec_len
r=assoc(lunr,fla_fms_sky_st(nrec))
srec=strcompress(string(nrec))
print,'Reading'+srec+' HISL Sky Records'
f1=r(0)
close,lunr & free_lun,lunr
;

; Fetch HIFA Skymap Records
; -------------------------
fname='csdr$firas_in:fms_sky_hifa.fad_8932800_9026409'
openr,lunr,fname,/get_lun,/share
finfo=fstat(lunr)
nrec=finfo.size/finfo.rec_len
r=assoc(lunr,fla_fms_sky_st(nrec))
srec=strcompress(string(nrec))
print,'Reading'+srec+' HIFA Sky Records'
f2=r(0)
close,lunr & free_lun,lunr
;

; Initialize Arrays
; -----------------
num1=fltarr(6144) & num2=num1 & combined_chi_sq=num1
res1=fltarr(160) & res2=res1
s12=fltarr(160) & denom=s12
;

; Pixel Numbers and Weights for HISL
; ----------------------------------
pix1=f1.attitude.pixel_no
n1x=f1.coad_spec_head.comb_num_ifgs
num1(pix1(*)) = n1x(*)
;

; Pixel Numbers and Weights for HIFA
; ----------------------------------
pix2=f2.attitude.pixel_no
n2x=f2.coad_spec_head.comb_num_ifgs
num2(pix2(*)) = n2x(*)
;

; Find the Sky Pixels Common to Both Skymaps
; ------------------------------------------
FLA_MATCH,pix1,pix2,idx1,idx2
;

; Array of Matched Pixels
; -----------------------
pixel = pix1(idx1)
;

; Weights, Variances, and Real_Spectra for Matching Pixels
; --------------------------------------------------------
n1 = n1x(idx1)  &  n2 = n2x(idx2)
s1 = float(f1(idx1).spec_data.spec(4:163))
cvar1 = f1(idx1).spec_data.cvs_var(4:163)
s2 = float(f2(idx2).spec_data.spec(4:163))
cvar2 = f2(idx2).spec_data.cvs_var(4:163)
;

; For Each Matched Pixel, Compute Combined_Chi_Squared
; ----------------------------------------------------
cv1 = fltarr(160)  &  cv2 = cv1  &  s1k = cv1  &  s2k = cv1
nmatch = n_elements(idx1)
FOR kk = 0,nmatch-1 DO BEGIN
;
 denom(*) = (n1(kk)*w1(*))+(n2(kk)*w2(*))   ; N1*W1 + N2*W2
;
 s1k(*)=s1(*,kk)*g1(*)-o1(*)             ; Corrected HISL Spectrum
 s2k(*)=s2(*,kk)*g2(*)-o2(*)             ; Corrected HIFA Spectrum
 s12=((n1(kk)*w1*s1k)+(n2(kk)*w2*s2k))/denom   ; Merged HIGH Spectrum
;
 res1 = s1k-s12   ;  HISL Residual
 res2 = s2k-s12   ;  HIFA Residual
;
 cv1(*) = cvar1(*,kk)   ; HISL CVS_Variance
 cv2(*) = cvar2(*,kk)   ; HIFA CVS_Variance
;
 csqx = (res1*res1/cv1/g1/g1)+(res2*res2/cv2/g2/g2) ; Chi_Squared(nu)
;
 combined_chi_sq(pixel(kk))=total(csqx(fidx)) ; Chi_Squared combined over nu
;
ENDFOR
;

; Create IDL Save Set
; -------------------
sname='csdr$firas_save:fla_csq_high.iss'
save,filename=sname,id,num1,num2,freq_combined,combined_chi_sq
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Define and fill in the FIP_CSQ record
; -------------------------------------
print,' '
print,'Filling in the FIP_CSQ record.'
print,' '
csq_rec = FLA_FIP_CSQ(6144)
;
; Loop over the pixels
; --------------------
FOR j=0,6143 DO BEGIN
  csq_rec(j).outmap = id
  csq_rec(j).inmap_1 = map1
  csq_rec(j).inmap_2 = map2
  csq_rec(j).pixel = j
  csq_rec(j).weight_1 = num1(j)
  csq_rec(j).weight_2 = num2(j)
  csq_rec(j).num_freq = num_combined
  csq_rec(j).comb_chi_square = combined_chi_sq(j)
ENDFOR
;

; Write the FIP_CSQ file
; ----------------------
print,' '
print,'Writing the FIP_CSQ_HIGH file.'
print,' '
;
outfile1 = 'csdr$firas_out:fip_csq_' + id + '.F16_93HYBRID'
rec_len = 32                 ; fixed length record size in bytes
;
OPENW,10,outfile1, rec_len, /fixed
FOR j=0,6143 DO WRITEU,10,csq_rec(j)
CLOSE,10
;
print,'File "'+strupcase(outfile1)+'" Written.'
print,' '
;

; Set Return Status
; -----------------
error = 0
print,'FLA_CSQ_HIGH : Successful Completion'
print,' '
;

RETURN
END
