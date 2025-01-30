Pro FMD_MATCH,a,b,suba,subb,long_flag,error
;+
; NAME:
;   FMD_MATCH
; PURPOSE:
;   Routine to match values in two vectors.
; CALLING SEQUENCE:
;	match,a,b,suba,subb,long_flag,error
; INPUTS:
;	a,b - two vectors to match elements
;       long_flag : "Y(es)" to make longword arrays
; OUTPUTS:
;	suba - subscripts of elements in vector a with a match
;		in vector b
;	subb - subscripts of the positions of the elements in
;		vector b with matchs in vector a.
;       error - Return status
;
;	suba and subb are ordered such that a(suba) equals b(subb)
; SIDE EFFECTS:
;	!ERR is set to the number of matches.
; RESTRICTIONS:
;	a and b should not have duplicate values within them.
;	You can use rem_dup function to remove duplicate values
;	in a vector
; HISTORY:
;	D. Lindler  Mar. 1986. (MATCH.PRO)
;       K. Jensen, Hughes STX, added a flag to make indices longword
;                              instead of integer.
;       K. Jensen, Hughes STX, added a return status flag, renamed MATCH2.
;       K. Jensen, Hughes STX, if matched indices contain duplicate
;                              values, trim them.
;       K. Jensen, Hughes STX, 07-MAY-1997, renamed FMD_MATCH.
;-
;-------------------------------------------------------------------------
;
error = 1  ;  Return status
;
if N_Params() ne 6 then begin
 print,'FMD_MATCH,a,b,suba,subb,long_flag,error'
 return
endif
;
long_flag = strupcase(strmid(long_flag,0,1))  ;  Longword option
;
na=n_elements(a)		;number of elements in a
nb=n_elements(b)		;number of elements in b
if(long_flag ne 'Y')then begin
 inda=indgen(na)			;indices in a
 indb=indgen(nb)			;indices in b
 vec=[intarr(na),replicate(1,nb)]       ;flag of which vector in combined list
endif
if(long_flag eq 'Y')then begin
 inda=lindgen(na)			;indices in a
 indb=lindgen(nb)			;indices in b
 vec=[lonarr(na),replicate(1,nb)]       ;flag of which vector in combined list
endif
c=[a,b]				;combined list of a and b
ind=[inda,indb]			;combined list of indices
				;	0 - a   1 - b
;
; sort combined list
;
sub=sort(c)
c=c(sub)
ind=ind(sub)
vec=vec(sub)
;
; find duplicates in sorted combined list
;
n=na+nb				;total elements in c
firstdup=where( (c eq shift(c,-1)) and (vec ne shift(vec,-1)))
if !err lt 1 then begin		;any found?
	suba=intarr(1)-1
	subb=intarr(1)-1
	return
end
dup=lonarr(!err*2)		;both duplicate values
if(long_flag ne 'Y')then even=indgen(n_elements(firstdup))*2
if(long_flag eq 'Y')then even=lindgen(n_elements(firstdup))*2
dup(even)=firstdup
dup(even+1)=firstdup+1
ind=ind(dup)			;indices of duplicates
vec=vec(dup)			;vector id of duplicates
suba=ind(where(vec eq 0))	;a subscripts
subb=ind(where(vec eq 1))	;b subscripts
;

; If duplicated indices, trim them
;
; Index 1
na = N_ELEMENTS(a)
nsa = N_ELEMENTS(suba)
IF (nsa gt na) THEN BEGIN
 subb = subb(SORT(suba))
 suba = suba(SORT(suba))
 flag2 = intarr(nsa)
 FOR i=long(1),long(nsa)-1 DO IF(suba(i) le suba(i-1))THEN flag2(i) = 1
 nx = WHERE(flag2 eq 0,cx)
 IF (cx le 0) THEN BEGIN
  PRINT,'MATCH : Error in Trimming Index 1 !'
  RETURN
 ENDIF
 suba = suba(nx)  &  subb = subb(nx)
ENDIF
;
; Index 2
nb = N_ELEMENTS(b)
nsb = N_ELEMENTS(subb)
IF (nsb gt nb) THEN BEGIN
 suba = suba(SORT(subb))
 subb = subb(SORT(subb))
 flag2 = intarr(nsb)
 FOR i=long(1),long(nsb)-1 DO IF(subb(i) le subb(i-1))THEN flag2(i) = 1
 nx = WHERE(flag2 eq 0,cx)
 IF (cx le 0) THEN BEGIN
  PRINT,'MATCH : Error in Trimming Index 2 !'
  RETURN
 ENDIF
 suba = suba(nx)  &  subb = subb(nx)
ENDIF
 
error = 0
;

RETURN
END
