function conv_gmt_sec, gmt		;convert GMT byte arrays to seconds
						;Frank Varosi NASA/GSFC 1990

; gmt = byte array ( 14, ntimes) of ASCII character codes for
;				 format YYDDDHHMMSSTTT  (TTT is milliseconds).

; function result = array of seconds, relative to 1/1/1989.

	jdtb = transpose( gmt - 48B )

	s = size( jdtb )
	nchar = s(2)

	if (nchar LT 14) then  jdtb = [ [jdtb] , [bytarr( s(1), 14-nchar )] ]

	pt2 = fix( [10,1] )
	pt3 = fix( [100,10,1] )    ;fix to keep everything in 16-bit integers.

	jdt = [ [ fix( jdtb(*,0:1) # pt2 ) ] ,	$	; years 
		[ fix( jdtb(*,2:4) # pt3 ) ] ,	$	; days
		[ fix( jdtb(*,5:6) # pt2 ) ] ,	$	; hours
		[ fix( jdtb(*,7:8) # pt2 ) ] ,	$	; minutes
		[ fix( jdtb(*,9:10) # pt2 ) ],	$       ; seconds
		[ fix( jdtb(*,11:13) # pt3 ) ]	]	; milliseconds

	hour = double( 3600. ) 
	conv = [ hour*24*365 , hour*24 , hour , 60. , 1. , 1.d-3 ]

	Start_Year = 89

	jdt(0,0) = jdt(*,0) - Start_Year
	seconds = jdt # conv

	if N_elements( seconds ) EQ 1 then seconds = seconds(0)

return, seconds
end
