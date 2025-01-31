;___________________________________________________________________________
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FEF_8
;
;DESCRIPTION:  
;     A driver to process FEF HISL-HIFA Comparison
;
;CALLING SEQUENCE:  
;     This main block of IDL code is invoked directly from IDL.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;   ERROR     (O) =  Error Status                          = 0 or 1
;
;WARNING:
;   The Logicals CSDR$FIRAS_IN and CSDR$FIRAS_SAVE must be previously defined.
;
;   CSDR$FIRAS_IN must contain the input HISL.ISS and HIFA.ISS IDL Save Sets.
;                 These must be created by FTT prior to executing this code.
;   CSDR$FIRAS_SAVE will contain the output FEF_8.ISS IDL Save Set.
;
;EXAMPLE:
;   $ UPRIDL
;   FEF_8,error
;
;#
;COMMON BLOCKS:
;   None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;   Calls  FEF_SPECTRA_XX.PRO
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;   None
;  
;MODIFICATION HISTORY:
;   Written by Ken Jensen, Hughes STX, Nov 28, 1994
; 
;-
;___________________________________________________________________
;
Pro FEF_8,error
;

; Set Error Status
; ----------------
error = 1
;

; Procedure Invoked Correctly ?
; -----------------------------
if N_Params() ne 1 then begin
 print,'FEF_8,error'
 return
endif
;

; Restore Save Sets
; -----------------
idx=indgen(34)
restore,'csdr$firas_in:hisl.iss'
n1 = n_hisl  &  s1 = s_hisl  &  c1 = c_hisl
;
restore,'csdr$firas_in:hifa.iss'
n2 = n_hifa  &  s2 = s_hifa  &  c2 = c_hifa
;

; Call FEF_SPECTRA_XX
; -------------------
FEF_SPECTRA_XX,s1,s2,n1,n2,c1,c2,diff=diff,ratio=ratio
;

; FEF IDL Save Set
; ----------------
id1='HISL' & id2='HIFA'
sname='csdr$firas_save:fef_8.iss'
save,filename=sname,id1,id2,diff,ratio
print,' '
print,'IDL Save Set "'+strupcase(sname)+'" Created.'
print,' '
;

; Re-Define Unused Restored Parameters
; ------------------------------------
b_hisl=0 & d_hisl=0 & b_hifa=0 & d_hifa=0
;

; Return Status
; -------------
error = 0
;

;
return
end
