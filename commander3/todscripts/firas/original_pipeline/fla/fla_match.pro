Pro FLA_MATCH,a,b,suba,subb
;+
; NAME:
;   FLA_MATCH
; PURPOSE:
;   Routine to match values in two vectors.
; CALLING SEQUENCE:
;	match,a,b,suba,subb
; INPUTS:
;	a,b - two vectors to match elements
; OUTPUTS:
;	suba - subscripts of elements in vector a with a match
;		in vector b
;	subb - subscripts of the positions of the elements in
;		vector b with matchs in vector a.
;
;	suba and subb are ordered such that a(suba) equals b(subb)
; SIDE EFFECTS:
;	!ERR is set to the number of matches.
; RESTRICTIONS:
;	a and b should not have duplicate values within them.
;	You can use rem_dup function to remove duplicate values
;	in a vector
; HISTORY:
;	Written by D. Lindler as procedure MATCH Mar. 1986.
;       Renamed FLA_MATCH by Ken Jensen, Hughes STX, 19 January 1995.
;       
;-
;-------------------------------------------------------------------------
na=n_elements(a)		;number of elements in a
nb=n_elements(b)		;number of elements in b
inda=indgen(na)			;indices in a
indb=indgen(nb)			;indices in b
c=[a,b]				;combined list of a and b
ind=[inda,indb]			;combined list of indices
vec=[intarr(na),replicate(1,nb)];flag of which vector in combined list
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
even=indgen(n_elements(firstdup))*2
dup(even)=firstdup
dup(even+1)=firstdup+1
ind=ind(dup)			;indices of duplicates
vec=vec(dup)			;vector id of duplicates
suba=ind(where(vec eq 0))	;a subscripts
subb=ind(where(vec eq 1))	;b subscripts
return
end

