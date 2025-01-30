	Integer*4 Function FUT_Combine_HiLo (rdghi, rdglo, trantemp, tranhwid,
	2                                    wt)

C-----------------------------------------------------------------------------
C      PURPOSE:	Subroutine to combine high- and low-current readings of an
C		individual GRT temperature.  Was originally FIT_Combine_Temps.
C
C      AUTHOR:  Fred Shuman
C               STX
C               1989 October 20
C    ...following the method used in FUT_Temperature_List, which, unlike this
C    routine, combines all the GRTs at once, and does so for each Major Frame.
C    ...1991 Dec 3:  It now follows the method of FUT_Temp_List.
C
C      CALLING SEQUENCE:
C	    status = FUT_Combine_HiLo (rdghi, rdglo, trantemp, tranhwid, wt)
C
C      CALLING ARGUMENTS:
C	      NAME        TYPE  I/O    DESCRIPTION
C	    rdghi          R*4   I     raw temperature reading (high current)
C	    rdglo          R*4   I     raw temperature reading (low current)
C	    trantemp       R*4   I     hi/lo current temp. transition midpt
C	    tranhwid       R*4   I     hi/lo current temperature transition
C	                                  half-width
C	    wt             R*4   O     weight for combining temperatures
C
C      INCLUDE FILES:
C	    CCT_Status_Record
C	    CCT_Get_Config
C	    $SSDef
C
C      SUBROUTINES CALLED:
C	    CT_GMT_to_Binary
C	    LIB$Signal
C=============================================================================
C   Modification History:
C
C   1990 Jan 25, D. Bouler, STX. SPR 5533
C       Output statement overflows record. An internal write was used to
C       convert a real number to a character string. Since the character
C       string was only 5 characters long, and the format for conversion
C       was F10.2 , there was a record statement overflow. The character
C       strings were re-defined to 7 characters, and format to 6.2 .
C       Also changed definition of message FIT_BADTEMPS in FIT_MSG because
C       it had one parameter that was not being used.
C
C    1991 Sep 25, F. Shuman, STX, SPRs 9078, 9079.
C       Rename routine from FIT_Combine_Temps to FUT_Combine_HiLo and move into
C       FUT to be used by both FIT and FEP.
C
C    1991 Dec 3, F. Shuman, Hughes STX, SPR 9327.
C       Make the method for combining hi/lo current rdgs in FUT_Combine_HiLo
C       consistent with that in FUT_Temp_List.
C-----------------------------------------------------------------------------

	Implicit None

	Include '(CCT_Status_Record)'
	Include '(CCT_Get_Config)'
	Include '($SSDEF)'

C  Calling arguments:
	Real      * 4      rdghi,rdglo	!temp to combine--high, low current rdgs
	Real      * 4      trantemp	! hi/lo-current temp. transition midpt
	Real      * 4      tranhwid	! hi/lo-current temperature transition
					!  interval radius (half-width)
	Real      * 4      wt		!weighting factor returned for combining
	                                ! temps

	Real      * 4      tranlo, tranhi	!low, high transition endpoints

	External           FUT_Normal

C------------------------------------------------------------------------------C
C  End of declarations, start of code..........................................C
C------------------------------------------------------------------------------C

	FUT_Combine_HiLo = %Loc(FUT_Normal)

	tranlo = trantemp - tranhwid
	tranhi = trantemp + tranhwid

	If (rdghi .Lt. -999.) Then
	   If (rdglo .Lt. -999.) Then
	      wt = -9999.
	   Else
	      wt = 0.
	   End If
	Else If (rdglo .Lt. -999.) Then
	   wt = 1.
	Else If (rdglo .Lt. tranlo) Then
	   wt = 0.
	Else If (rdglo .Gt. tranhi) Then
	   wt = 1.
	Else
	   wt = 0.5*((rdglo - trantemp)/tranhwid + 1.0)
	End If

	Return
	End
