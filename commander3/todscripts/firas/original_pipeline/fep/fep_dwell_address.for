	subroutine fep_dwell_address (dw_byte, iadrs)
c
c  Decode dwell flag and address from status byte. Bits 0-4 are the dwell
c  address (0-31), bit 7 is the dwell flag (0=off, 1=on). Technique is to
c  convert the byte to an integer*2 number and check its value, using
c  integer division and multiplication to mask bits.
c
	byte dw_byte, null(2)
	integer*2 iadrs, equiv, onoff
	equivalence (null(1), equiv)
	data null/0,0/
	null(1) = dw_byte
c
c  Byte is now an integer. Bit 7 corresponds to a value of 128, so dividing
c  by 128 gives a value of 0 or one depending on whether the dwell bit is
c  off or on. If off, tell the main pgm by returing a negative dwell address.
c
	onoff = equiv/128
	if (onoff .eq. 0) then
	   iadrs = -64
	   return
	endif
c
c  Dwell bit is on. Mask bits 0-4 to get the dwell address 0-31 by subtraction.
c
	onoff = equiv/32
	iadrs = equiv - 32*onoff
	return
	end
