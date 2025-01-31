function conv_adt_sec, adt, REF_ADT=refadt

	if N_elements( refadt ) NE 1 then refadt=9.56059992315d6 ;ref to 1989.
                                                
	st32 = (2L^30)*4.d-7		;# secs in 32-bit 100-nanosec counter.

	shi = (adt(1,*) - refadt) * st32
	slo = ishft( adt(0,*), -1) * 1.d-6 / 5

return, transpose( shi + slo )
end
