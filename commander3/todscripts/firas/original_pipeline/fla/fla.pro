Pro FLA,chanscan,error

;
error = 1
;

if N_Params() ne 2 then begin
 PRINT,'FLA,chanscan,error'
 RETURN
endif
;

chanscan=strupcase(chanscan)
if(chanscan eq 'RHSF')then chanscan='RHLF'
if(chanscan eq 'RLFL')then chanscan='RLSF'
if(chanscan eq 'LHSF')then chanscan='LHLF'
if(chanscan eq 'LLFL')then chanscan='LLSF'
if(chanscan eq 'RHSS')then idd='RHS'
if(chanscan eq 'RHLF')then idd='RHF'
if(chanscan eq 'RLSS')then idd='RLS'
if(chanscan eq 'RLSF')then idd='RSF'
if(chanscan eq 'RLLF')then idd='RLF'
if(chanscan eq 'LHSS')then idd='LHS'
if(chanscan eq 'LHLF')then idd='LHF'
if(chanscan eq 'LLSS')then idd='LLS'
if(chanscan eq 'LLSF')then idd='LSF'
if(chanscan eq 'LLLF')then idd='LLF'
;

n_freq = 0
ftype = strmid(idd,1,2)
if((ftype eq 'HS')or(ftype eq 'HF'))then n_freq = 167
if((ftype eq 'LS')or(ftype eq 'SF'))then n_freq = 34
if(ftype eq 'LF')then n_freq = 148
;

; Call FLA_FIT driver
; -------------------
PRINT,'Calling FLA_FIT Driver'
IF (n_freq eq 34) THEN FLA_FIT_LLR,chanscan,error
IF (n_freq eq 148) THEN FLA_FIT_LHR,chanscan,error
IF (n_freq eq 167) THEN FLA_FIT_HI,chanscan,error
;

if(error ne 0)then begin
 PRINT,'Error Returned from FLA_FIT !'
 RETURN
endif
;

; Call FLA_TCMBR (low frequency channels only)
; --------------------------------------------
IF (n_freq eq 167) THEN RETURN
;
PRINT,'Calling FLA_TCMBR'
FLA_TCMBR,chanscan,error
;
if(error ne 0)then PRINT,'Error Returned from FLA_TCMBR !'
;

RETURN
END
