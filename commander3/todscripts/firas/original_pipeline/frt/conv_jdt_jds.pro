function conv_jdt_jds, jdt	;convert Julian Date-Time integers to strings
						;Frank Varosi NASA/GSFC 1990

; jdt = matrix of integers ( 6, Ndates ), (or just one row), where the
;   format of each row of matrix: [ year, day, hour, min, sec, millisec ] .

; function result = array of string(s): YYDDDHHMMSSTTT .

	time = transpose( jdt )
	T1dig = byte( time MOD 10 )
	T2dig = byte( (time/10) MOD 10 )

	jdtb = [ [ T2dig(*,0) ], [ T1dig(*,0) ], $			;year
		[byte( time(*,1)/100 )], [T2dig(*,1)], [T1dig(*,1)], $	;day
		[ T2dig(*,2) ], [ T1dig(*,2) ], $			;hour
		[ T2dig(*,3) ], [ T1dig(*,3) ], $			;minut.
		[ T2dig(*,4) ], [ T1dig(*,4) ], $			;secs.
		 [byte( time(*,5)/100 )], [T2dig(*,5)], [T1dig(*,5)] ]	;milli.

return, string( transpose( jdtb + 48B ) )
end
