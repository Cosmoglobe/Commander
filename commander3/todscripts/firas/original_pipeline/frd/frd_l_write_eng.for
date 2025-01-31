
	Integer*4 Function FRD_L_Write_Eng ( Englim_Rec, Lun )

!	Program Name : FRD_L_Write_Eng
!
!	Programmer: Shirley M. Read, STX, January 1988

!	Program Description:
!	    This subroutine writes the record containing the red and yellow
!	    engineering limits into the Fex_Englim.Dat file.
!
!----------------------------------------------------------------------------

	Implicit None

!	Passed Parameters.

	dictionary	'Fex_Englim'

	record		/Fex_Englim/englim_rec

	integer*4 lun   ! Logical unit number

!	Include files

	Include '($SSdef)'

!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 systatus      ! I/O status
	Integer*4 success / 1 /, error / 2 /

!	Set status to success.

	retstat = success

!       Open the output file, Fex_Englim.Dat 

	    Open (unit=lun, file='Fex_Englim.Dat',
	1 	status='new', form='unformatted', access='sequential',
	2	organization='sequential', iostat=systatus,
	3	recordsize=1024, recordtype='fixed')

	    If (systatus .ne. 0) then
!	      Call Lib$Signal(fdq_openerr,%val(2),%val(systatus),
!	1	'Fex_Englim.Dat')
	      Write(6,100) systatus
	      retstat = error
	    Endif
 100	format(1x,'Error: Opening Fex_Englim.Dat. Status= ',
	1	z8.8)

!	Write the FEX_ENGLIM record containing the red low limits,
!	the yellow low limits, yellow high limits and red high limits,
!	respectively, from ENGLIM_REC, the record structured buffer.

	If ( retstat.eq.success ) then
	  Write( unit=lun, iostat=systatus ) englim_rec    
  	  if ( systatus.ne.0) then
!	    Call Lib$Signal(fdq_writerr,%val(2),%val(systatus),
!	1	'FEX_Englim.Dat')
	    Write(6,200) systatus
	    retstat = error
	  endif
	Endif
 200	format(1x,'Error: Writing Fex_Englim.Dat. Status= ',
	1	z8.8)

!	Close Fex_Englim.Dat

	If ( retstat .eq. success ) then
	  Close (unit=lun, iostat=systatus)
 	  if ( systatus.ne.0 ) then
!	    Call Lib$Signal(fdq_closerr,%val(2),%val(systatus),
!	1	'FEX_Englim.Dat')
	    Write(6,400) systatus
	    retstat = error
	  else
	    Write(6,500)
	  endif
	Endif
 400	format(1x,'Error: Closing Fex_Englim.Dat. Status= ',        
	1	z8.8)
 500    format(1x,'Success: Record written to ', 
	1	'Fex_Englim.Dat')    

	FRD_L_Write_Eng = retstat

	return

	end
