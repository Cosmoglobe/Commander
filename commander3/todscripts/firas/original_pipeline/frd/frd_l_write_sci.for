

	Integer*4 Function FRD_L_Write_Sci ( Scilim_Rec, Lun )

!	Program Name : FRD_L_Write_Sci
!
!	Programmer: Shirley M. Read, STX, January 1988

!	Program Description:
!	    This subroutine writes the record containing the red and yellow
!	  science limits into the Fex_Scilim.Dat file.
!
!----------------------------------------------------------------------------

	Implicit None

!	Passed Parameters.

	dictionary	'Fex_Scilim'

	record		/Fex_Scilim/scilim_rec

	integer*4 lun   ! Logical unit number

!	Include files

	Include '($SSdef)'

!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 systatus      ! I/O status
	Integer*4 success / 1 /, error / 2 /

!	Set status to success.

	retstat = success

!       Open the output file, Fex_Scilim.Dat 

	    Open (unit=lun, file='Fex_Scilim.Dat',
	1 	status='new', form='unformatted', access='sequential',
	2	organization='sequential', iostat=systatus,
	3	recordsize=256, recordtype='fixed')

	    If (systatus .ne. 0) then
!	      Call Lib$Signal(fdq_openerr,%val(2),%val(systatus),
!	1	'Fex_Scilim.Dat')
	      Write(6,100) systatus
	      retstat = error
	    Endif
 100	format(1x,'Error: Opening Fex_Scilim.Dat. Status= ',
	1	z8.8)

!	Write the FEX_SCILIM record containing the red limits and
!	the yellow limits from SCILIM_REC, the record structured buffer.

	If ( retstat.eq.success ) then
	  Write( unit=lun, iostat=systatus ) scilim_rec    
  	  if ( systatus.ne.0) then
!	    Call Lib$Signal(fdq_writerr,%val(2),%val(systatus),
!	1	'FEX_Scilim.Dat')
	    Write(6,300) systatus
	    retstat = error
	  endif
	Endif
 300	format(1x,'Error: Writing Fex_Scilim.Dat. Status= ',
	1	z8.8)

!	Close Fex_Scilim.Dat

	If ( retstat .eq. success ) then
	  Close (unit=lun, iostat=systatus)
 	  if ( systatus.ne.0 ) then
!	    Call Lib$Signal(fdq_closerr,%val(2),%val(systatus),
!	1	'FEX_Scilim.Dat')
	    Write(6,400) systatus
	    retstat = error
	  else
	    Write(6,500)
	  endif
	Endif
 400	format(1x,'Error: Closing Fex_Scilim.Dat. Status= ',
	1	z8.8)
 500    format(1x,'Success: Record written to ',
	1	'Fex_Scilim.Dat' )

	FRD_L_Write_Sci = retstat

	Return
	End
