   function FMD_choleskey,a
;
;	FDS_CHOLESKEY written by rashafer 94 aug 17
;	This is my form of the cholesky decomposition into a format that is
;	more useful for subsequent matrix manipulation in IDL.  Also it
;	doesn't have the worrysome aspect of mucking with the elements
;	that shouldn't change, but which NR_CHOLDC apparently does modify
;	(at small levels) anyway.
;		A - a real or double symmeteric positive definite matrix.
;	Returns
;		MY_CHOLDC - a lower triangular matrix, L, such that LL^T
;			equals A.
;
;       Renamed FMD_CHOLESKEY, K.Jensen, Hughes STX, 21-Jun-96
;
;
	sizestuff = size(a)
	ndim = sizestuff(0)
	L = make_array(size=sizestuff)
	one = L(0,0)+1.

	n = sizestuff(1)-1
	L(0,0) = sqrt(a(0,0))
	linv = one/l(0,0)
	l(1:n,0) = linv*a(0,1:n)
	for i = 1, n do begin
	    l(i,i) = sqrt(a(i,i)-total(l(i,0:i-1)*l(i,0:i-1)))
	    linv = one/l(i,i)
	    for j = i, n do begin
		l(j,i) = linv*(a(i,j)-total(l(i,0:i-1)*l(j,0:i-1)))
	      endfor
	  endfor
	return,L
    end
