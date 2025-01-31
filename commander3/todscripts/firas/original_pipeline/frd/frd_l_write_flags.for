
	Integer*4 Function FRD_L_Write_Flags ( Limflags, Lun )

!	Program Name : FRD_L_Write_Flags
!
!	Programmer: Shirley M. Read, STX, January 1988

!	Program Description:
!	    This subroutine writes the record containing the enable/
!	    disable flags for the checking of limits for 110 quality flags 
!	    into the Fex_Limflags.Dat file.
!
!----------------------------------------------------------------------------

	Implicit None

!	Passed Parameters.

	Dictionary	'Fex_Limflags'

	Record		/Fex_Limflags/limflags

	Integer*4 lun   ! Logical unit number

!	Include files

	Include '($SSdef)'

!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 systatus      ! I/O status
	Integer*4 success / 1 /, error / 2 /

!	Set status to success.

	retstat = success

!       Open the output file, Fex_Limflags.Dat 

	    Open (unit=lun, file='Fex_Limflags.Dat',
	1 	status='new', form='unformatted', access='sequential',
	2	organization='sequential', iostat=systatus,
	3	recordsize=128, recordtype='fixed')

	    If (systatus .ne. 0) then
!	      Call Lib$Signal(fdq_openerr,%val(2),%val(systatus),
!	1	'Fex_Limflags.Dat')
	      Write(6,100) systatus
	      retstat = error
	    Endif
 100	format(1x,'Error: Opening Fex_Limflags.Dat. Status= ',
	1	z8.8)

!	Write the FEX_LIMFLAGS record containing the enable/disable flags
!	for limit checking from LIMFLAGS, the record structured buffer.

	If ( retstat.eq.success ) then
	  Write ( unit=lun, iostat=systatus ) limflags    
  	  if ( systatus.ne.0) then
!	    Call Lib$Signal(fdq_writerr,%val(2),%val(systatus),
!	1	'FEX_Limflags.Dat')
	    Write(6,200) systatus
	    retstat = error
	  endif
	endif
 200	format(1x,'Error: Writing Fex_Limflags.Dat. Status= ',
	1	z8.8)

!	Close Fex_Limflags.Dat

	If ( retstat .eq. success ) then
	  Close (unit=lun, iostat=systatus)
 	  if ( systatus.ne.0 ) then
!	    Call Lib$Signal(fdq_closerr,%val(2),%val(systatus),
!	1	'FEX_Limflags.Dat')
	    Write(6,400) systatus
	    retstat = error
	  else
	    Write(6,500)
	  endif
	Endif
 400	format(1x,'Error: Closing Fex_Limflags.Dat. Status= ',
	1	z8.8)
 500    format(1x,'Success: Record of Fex_Limflags.Dat',
	1	' written.') 

	FRD_L_Write_Flags = retstat

	Return
	End
