
	Integer*4 Function FRD_L_Store_Sci ( Parm, Lenparm, Limits,
	1	  Scilim_Rec )

!	Program Name : FRD_L_Store_Sci
!
!	Programmer: Shirley M. Read, STX, January 1988

!	Program Description:
!	    This subroutine stores the red and yellow science limits into 
!	  the scilim_rec structured buffer if the parameter name matches 
!	  one of the RDL field names for the Fex_Scilim record.
!	  An assumption is made that the limits for angles are input as degrees.
!
!	Modified:   Shirley M. Read, STX, May 1988
!	   The record structure for science limits was modified to include
!          all four limits in one record. The new COBETRIEVE 'Get_Config'
!          will thus be able to access one record for each call instead of 
!          four. This routine was modified to use the new data dictionary
!          definition for the science limits.
!
!----------------------------------------------------------------------------

	Implicit None

!	Passed Parameters.

	Character*20	Parm	! First parameter
	Integer*4	Lenparm   ! Length of first character parm
	Integer*2	Limits(2) ! Red and yellow limits

	Dictionary	'Fex_Scilim'

	Record		/Fex_Scilim/Scilim_Rec

!	Local variables

	Integer*4 retstat	! Program Processing status
	Integer*4 success / 1 /, error / 2 /
	Integer*4 ix            ! Index
	Integer*4 rec		! Index for record structure
	Integer*4 match         ! Var to store index of match
	Character*20 Match_Table(10)
	Data Match_Table(1) /'SC_HEAD20           '/
	Data Match_Table(2) /'GLIRAT_RH           '/
	Data Match_Table(3) /'GLIRAT_RL           '/
	Data Match_Table(4) /'GLIRAT_LH           '/
	Data Match_Table(5) /'GLIRAT_LL           '/
	Data Match_Table(6) /'EARTH_LIMB          '/
	Data Match_Table(7) /'SUN_ANGLE           '/
	Data Match_Table(8) /'MOON_ANGLE          '/
	Data Match_Table(9) /'ALTITUDE            '/
	Data Match_Table(10) /'DATA_READY          '/
	Real*4 Pi
	Parameter ( Pi = 3.14159265359 )
	Real*4 D180
	Parameter ( D180 = 180.0 )
	Real*4 Radangle		! Angle in radians
	Integer*2 Iradangle     ! Integer value of angle in units of 
				! 10**-4 radians

!	Set status to success.

	retstat = success

!       Check the input parameter and set the match index if found.

	Match = 0
	Do ix = 1, 10
	  If ( Parm(1:lenparm) .eq. Match_Table(ix)(1:lenparm) ) then 
	    Match = ix
	  Endif
	Enddo

!	If match found, store the limits in the Scilim_Rec.

	If ( Match .ne. 0 ) then
          If ( Match .eq. 1 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Sci_Head.SC_Head20 = limits(rec)
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 2 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Glitch_rate(1) = limits(rec)/100.
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 3 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Glitch_rate(2) = limits(rec)/100.
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 4 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Glitch_rate(3) = limits(rec)/100.
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 5 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Glitch_rate(4) = limits(rec)/100.
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 6 ) then
	    Do rec = 1, 2
	      radangle = ( Pi * limits(rec) ) / D180
	      iradangle = nint ( radangle * 10000.0 )
	      Scilim_Rec.Lim(rec).Attitude.Earth_Limb = iradangle
	      Write (6,200) parm(1:lenparm), iradangle
	    Enddo
	  Elseif ( Match .eq. 7 ) then
	    Do rec = 1, 2
	      radangle = ( Pi * limits(rec) ) / D180
	      iradangle = nint ( radangle * 10000.0 )
	      Scilim_Rec.Lim(rec).Attitude.Sun_Angle = iradangle
	      Write (6,200) parm(1:lenparm), iradangle
	    Enddo
	  Elseif ( Match .eq. 8 ) then
	    Do rec = 1, 2
	      radangle = ( Pi * limits(rec) ) / D180
	      iradangle = nint ( radangle * 10000.0 )
	      Scilim_Rec.Lim(rec).Attitude.Moon_Angle = iradangle
	      Write (6,200) parm(1:lenparm), iradangle
	    Enddo
	  Elseif ( Match .eq. 9 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Attitude.Altitude = limits(rec)
	      Write (6,300) parm(1:lenparm), limits(rec)
	    Enddo
	  Elseif ( Match .eq. 10 ) then
	    Do rec = 1, 2
	      Scilim_Rec.Lim(rec).Sci_Head.Data_Ready(1) = limits(rec)
	      Write (6,100) parm(1:lenparm), limits(rec)
	    Enddo
	  Endif
	Else
	  Write(6,400) Parm(1:lenparm)
	Endif

!	Formats
 100    format(1x,'Parameter ',a,' stored in Scilim_Rec. Limit = ',i10) 
 200    format(1x,'Parameter ',a,' stored in Scilim_Rec. Limit',
	1	'(Units of 10**-4 Radians) = ',i10)
 300    format(1x,'Parameter ',a,' stored in Scilim_Rec. Limit',
	1	'(Units of 0.1 Km) = ',i10) 
 400	format(1x,'Parameter ',a,' not found in Match Table.') 
	
	FRD_L_Store_Sci = retstat

	Return
	End
