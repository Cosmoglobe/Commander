	Integer*4 Function FRD_L_Store_Flags( FRLun, Flag, Flaglen, 
	1	  Number, Parm, Flagbuf )

!	Program Name : FRD_L_Store_Flags
!
!	Programmer: Shirley M. Read, STX, May 1988
!
!	Program Description:
!	    This subroutine stores the values of the limit checking flags
!	    into an array to be moved to the Fex_Limflag record structure.
!
!----------------------------------------------------------------------------

	Implicit None

!	Include files.

	Include '(FUT_Qualflag_Names)'

!	Passed Parameters.

	Character*15	Flag          ! Input flag field parameter
	Integer*2       Flaglen       ! Input length of flag name
	Integer*2       Number        ! Input number of values to be stored
	Byte    	Parm(20)      ! Input set of values of limit check flags
	Byte		Flagbuf(320)  ! Output array of limit flags
	Integer*4	FRLun	      ! Field retriever LUN

!	Functions.

	Integer*4 	FUT_Field_Attributes	! Get record field attribute

!	Local variables

	Character*64    fieldname     ! Record and field name
	Integer*2       offset        ! Offset of field in record
	Integer*2	length        ! Length of field
	Integer*4 	retstat	      ! Program Processing status
	Integer*4	status	      ! Return status
	Integer*4 	success / 1 /, error / 2 /
	Integer*4 	ix            ! Index
	Logical*1       found	      ! Flag for match found

!	Set status to success.

	retstat = success

!       Check for invalid flag name.

	found = .false.
	Do ix = 1, 110
	  If ( flag(1:flaglen) .eq. qualflag_names(ix)(1:flaglen)) then
	      found = .true.
	  Endif
	Enddo

!	Find the offset of the flag array or field for the output 
!	limflags record.

	If ( found ) then
	  fieldname = 'FEX_Limflags.Lim_Flags.' // Flag	   
	  status = FUT_Field_Attributes( frlun, fieldname, length,
     .					 offset)
	  if (.not. status) then
	    retstat = error
	    Write(6,150) status
	  endif
	Else	! Not found
	    retstat = error
	    Write(6,200) flag
	Endif

!	Check for valid number of flags to store for this field.

!	If ( number .gt. length ) then
!	  Write(6,250) number, flag(1:flaglen), length
!	  retstat = error
!	Endif

!	If status is success, store the flag values in the Limflags record.

	If ( retstat .eq. success ) then
	    Write(6,300) Flag, number
	    Do ix = 1, number
	      offset = offset + 1
	      Call Lib$Movc3(1, parm(ix), flagbuf(offset))
	      Write(6,350) ix, flagbuf(offset)
	    Enddo
	Endif

!	Formats

 150    format(1x,'Error from FUT_Field_Attributes. Status= ',z8.8)
 200    format(1x,'Error: Flag field name ',a,' is invalid.')
 250    format(1x,'Error: Input number of values',i3,' for flag ',a,
	1     ' is greater than field length of ',i4,' in record.') 
 300    format(1x,'Flag field ',a,' has limits stored in ',i3,
	1     ' positions of the flag array buffer.')
 350    format(10x,'Flag array position ',i3,' has hex value ',z2.2)

	FRD_L_Store_Flags = retstat

	Return
	End
