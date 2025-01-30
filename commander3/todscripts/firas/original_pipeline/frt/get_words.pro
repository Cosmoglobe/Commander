function get_words, text

; separate text string(s) into array of words and return it,
;  text = string(s) with words delimited by blanks or commas.
; Frank Varosi NASA/Goddard 1989

	if N_elements( text ) LE 0 then return,""

	if N_elements( text ) GT 1 then begin
		Lens = strlen( text )
		Lenmax = max( Lens, maxi )
		text(maxi) = text(maxi) + " "
	   endif

	textb = byte( text )
	Len = N_elements( textb )
	if (Len LE 1) then return, string( textb )

	bb = byte( " " )
	bb = bb(0)

	cb = byte( "," )
	cb = cb(0)
	compos = where( textb EQ cb, ncomma )
	if (ncomma GT 0) then textb(compos) = bb

	tb = byte( 9 )
	tabpos = where( textb EQ tb, ntab )
	if (ntab GT 0) then textb(tabpos) = bb

	zpos = where( textb EQ 0, nzero )
	if (nzero GT 0) then textb(zpos) = bb

	bpos = where( textb EQ bb, nblank )
	if (nblank LE 0) then return, string( textb )

	words = strarr( nblank+1 )
	tend = [ bpos, Len-1 ]
	bpos = [ 0, bpos ]

	for i=0,nblank do  words(i) = string( textb( bpos(i):tend(i) ) )

	words = strtrim( words, 2 )
	w = where( words NE "" )

return, words(w)
end
