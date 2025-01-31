function conv_sec_jdt, seconds		;convert seconds to J.D. integers 
					;Frank Varosi NASA/GSFC (STX) 1990

; seconds = array of seconds relative to 1/1/1989 (double prec. float. point)

; function result = matrix of integers ( 6, Ndates ), where the
;   format of each row of matrix: [ year, day, hour, min, sec, millisec ] .

	if N_elements( Start_Year ) NE 1 then Start_Year = 89

	Nt = N_elements( seconds )
	jdt = intarr( Nt, 6 )

	hour = double( 3600. ) 
	conv = [ hour*24*365 , hour*24 , hour , 60. , 1. , 1.d-3 ]
	fsec = seconds

	for i=0,5 do begin

		isec = fix( fsec/conv(i) )
		jdt(0,i) = isec
		fsec = fsec - isec*conv(i)
	  endfor

	jdt(0,0) = jdt(*,0) + Start_Year    	;time relative to Start_Year

return, transpose( jdt )
end
