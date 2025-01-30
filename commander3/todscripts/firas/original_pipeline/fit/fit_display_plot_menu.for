	Integer*4 Function FIT_Display_Plot_Menu ( page, menu,
     .			menu_num, menu_page, menu_element,
     .			nfields )

C--------------------------------------------------------------------------
C    PURPOSE: Display field menu and select from the menu.
C	      
C
C       REQUIREMENTS REFERENCE NUMBERS:
C
C       INPUT DATA:
C
C       OUTPUT DATA:
C
C    AUTHOR: Fred Shuman, STX	(adapted from FEP_Display_Plot_Menu
C            October 1, 1987	 by Rob Kummerer, STX, March 26, 1987)
C            version: Nov 2, 1987
C
C    INVOCATION: STATUS = FIT_DISPLAY_PLOT_MENU ( PAGE, MENU,
C     						  MENU_NUM, MENU_PAGE,
C						  MENU_ELEMENT,
C						  NFIELDS )
C
C    INPUT PARAMETERS:
C	PAGE			I*2		Menu page number.
C	MENU(*)			C*32		Choices on selected sub-menu.
C	MENU_NUM		I*4		Number of menu elements.
C
C    OUTPUT PARAMETERS: 
C	MENU_PAGE(*)		I*2		Selected menu page.		
C	MENU_ELEMENT(*)		I*2		Selected menu element.
C	NFIELDS			I*2		Number of fields.
C
C    SUBROUTINES CALLED: 
C	OTS$Cvt_TI_L
C	STR$UpCase
C	SMG$Create_Virtual_Display
C	SMG$Create_Virtual_Keyboard
C	SMG$Delete_Virtual_Display
C	SMG$Delete_Virtual_Keyboard
C	SMG$Paste_Virtual_Display
C	SMG$Draw_Rectangle
C	SMG$Draw_Line
C	SMG$Erase_Chars
C	SMG$Change_Rendition
C	SMG$Ring_Bell
C	SMG$Put_Line
C	SMG$Flush_Buffer
C	SMG$Read_String
C	SMG$Set_Cursor_Abs
C	SMG$Erase_PasteBoard
C
C    COMMON VARIABLES USED: None
C
C    INCLUDE FILES: 
C	FUT_Error
C	FIT_Menu
C	$SMGDEF
C	$TRMDEF
C	$SSDEF
C---------------------------------------------------------------------
C      Revise:
C         Spr 5655, Add the hard coded menu text file.
C             Version 5.8, H. Wang, STX, Mar 8, 1990.  
C
C         Spr 7694, 7706, Rename the calibration resistors.
C             Version 7.1, H. Wang, STX, Nov. 13, 1990.  
C----------------------------------------------------------------------

	Implicit	None

	Include		'(FUT_Error)'
	Include		'(FIT_Menu)'
	Include 	'($SMGDEF)'
	Include 	'($TRMDEF)'
	Include 	'($SSDEF)'

	Integer 	*4	i
	Integer 	*4	j
	Integer 	*4	k
	Integer 	*4	l
	Integer 	*4	m
	Character	*4	inp
	Integer		*4	off
	Integer		*4	irow
	Integer		*4	icol
	Integer		*4	pbid
	Integer		*4	vdid
	Integer		*4	kbid
	Integer		*4	status
	Character	*80	menu_line
	Logical		*1	display_menu
	Integer		*4	rend/2/
	Logical		*1	select
        Logical         *1      not_avail

	Integer		*2	page
	Character	*64	title(6)
	Character	*32	menu(*)
	Character	*32	selprmt
	Integer		*4	menu_num
	Integer		*2	nfields
	Integer		*2	menu_page(*)
	Integer		*2	menu_element(*)

	Integer		*4	OTS$Cvt_TI_L
	Integer		*4	STR$UpCase
	Integer		*4	SMG$Create_Virtual_Display
	Integer		*4	SMG$Create_Virtual_Keyboard
	Integer		*4	SMG$Delete_Virtual_Display
	Integer		*4	SMG$Delete_Virtual_Keyboard
	Integer		*4	SMG$Paste_Virtual_Display
	Integer		*4	SMG$Draw_Rectangle
	Integer		*4	SMG$Draw_Line
	Integer		*4	SMG$Erase_Chars
	Integer		*4	SMG$Change_Rendition
	Integer		*4	SMG$Ring_Bell
	Integer		*4	SMG$Put_Line
	Integer		*4	SMG$Flush_Buffer
	Integer		*4	SMG$Read_String
	Integer		*4	SMG$Set_Cursor_Abs
	Integer		*4	SMG$Erase_PasteBoard

	External	FIT_Normal

	Common /PB/	pbid, irow, icol

	Data title /'A Side GRT Temperatures Menu',
     .		    'B Side GRT Temperatures Menu',
     .		    'Combined GRT Temperatures Menu',
     .		    'IPDU Temperatures Menu',
     .		    'Voltages Menu',
     .		    'Currents Menu'/

	Save rend

c
c Initialize.
c
	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Display ( irow, icol, vdid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Create_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

c
c Put out the title.
c
	status = SMG$Draw_Line ( vdid, 3, 1, 3, icol )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	i = Len(title(page))
	Do While (title(page)(i:i) .Eq. ' ')
	   i = i - 1
	End Do

	status = SMG$Set_Cursor_Abs ( vdid, 2, 40-i/2 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid, title(page) )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

c
c Put out the menu.
c
	off = menu_num / 2

	Do m=1,off

	   Write (menu_line, 30) m, menu(m), m+off, menu(m+off)
30	   Format ( i3, 2x, a, 5x, i3, 2x, a )

	   status = SMG$Set_Cursor_Abs ( vdid, m+3, 3 )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Put_Line ( vdid, menu_line )
	   If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	End Do

	status = SMG$Set_Cursor_Abs ( vdid, m+4, 5 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Put_Line ( vdid, 'E  Exit to main menu' //
     .			'   R  Reset menu' )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Draw_Rectangle ( vdid, 1, 1, irow, icol )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

c
c Display the menu.
c
	status = SMG$Paste_Virtual_Display ( vdid, pbid, 1, 1 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Flush_Buffer ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Set_Cursor_Abs ( vdid, irow-1, 3 )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	display_menu = .True.

c
c Flag the fields already selected from the current menu.
c
	Do i=1,nfields
	   If (menu_page(i) .Eq. page) Then
	      j = menu_element(i)
	      If (j .Le. off) Then
		 status = SMG$Change_Rendition ( vdid, j+3, 4,
     .							1, 2, 2 )
	         If (status .Ne. SS$_Normal) Then
	   	    FIT_Display_Plot_Menu = status
	   	    Call LIB$Signal(%Val(status))
		    Return
		 End If
	      Else
		 k = j - off
		 status = SMG$Change_Rendition ( vdid,  k+3, 46,
     .							1, 2, 2 )
	         If (status .Ne. SS$_Normal) Then
	   	    FIT_Display_Plot_Menu = status
	   	    Call LIB$Signal(%Val(status))
		    Return
		 End If
              End If
	   End If
	End Do

c
c Select from the menu.
c
	Do While (display_menu)

	   selprmt = 'Select a menu item:'
	   status = SMG$Read_String ( kbid, inp, selprmt,,,,,,,vdid )
           If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Erase_Chars ( vdid, icol-2, irow-1, 2 )
           If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = SMG$Set_Cursor_Abs ( vdid, irow-1, 3 )
           If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   status = STR$UpCase ( inp, inp )
           If (status .Ne. SS$_Normal) Then
	      FIT_Display_Plot_Menu = status
	      Call LIB$Signal(%Val(status))
	      Return
	   End If

	   If (inp(1:1) .Eq. 'E') Then
c
c Exit from the field select menu.
c
	      display_menu = .False.

	   Else If (inp(1:1) .Eq. 'R') Then
c
c Reset the menu.
c
	      j = 0
	      Do i=1,nfields
		 If (menu_page(i) .Ne. page) Then
		    j = j + 1
		    menu_page(j) = menu_page(i)
		    menu_element(j) = menu_element(i)
		 Else
		    l = menu_element(i)
		    If (l .Le. off) Then
		       status = SMG$Change_Rendition ( vdid, l+3, 2,
     .						1, 4, 0 )
                       If (status .Ne. SS$_Normal) Then
	         	  FIT_Display_Plot_Menu = status
	         	  Call LIB$Signal(%Val(status))
	         	  Return
	               End If
		    Else
		       k = l - off
		       status = SMG$Change_Rendition ( vdid,  k+3,
     .						44, 1, 4, 0 )
                       If (status .Ne. SS$_Normal) Then
	         	  FIT_Display_Plot_Menu = status
	         	  Call LIB$Signal(%Val(status))
	         	  Return
	               End If
		    End If
		 End If
	      End Do

	      nfields = j
	      Do i=nfields+1,total_menu_num
		 menu_page(i) = 0
		 menu_element(i) = 0
	      End Do

	   Else

c
c Field select from menu.
c
	      status = OTS$Cvt_TI_L(inp,j,%val(4),%val(17))
              If (status .Ne. SS$_Normal) Then
	         FIT_Display_Plot_Menu = status
	         Call LIB$Signal(%Val(status))
	         Return
	      End If

c
c	Check the selection.
c
              not_avail = .false.
              if (page .eq. grta_page) then
                if (j .eq. 1) not_avail = .true.
                if (j .eq. 2) not_avail = .true.
                if (j .eq. 3) not_avail = .true.
                if (j .eq. 4) not_avail = .true.
              endif 
              if (page .eq. grtb_page) then
                if (j .eq. 1) not_avail = .true.
                if (j .eq. 2) not_avail = .true.
                if (j .eq. 3) not_avail = .true.
                if (j .eq. 4) not_avail = .true.
                if (j .eq. 9) not_avail = .true.
                if (j .eq. 15) not_avail = .true.
                if (j .eq. 25) not_avail = .true.
                if (j .eq. 31) not_avail = .true.
              endif 
              if (page .eq. comb_page) then
                if (j .eq. 25) not_avail = .true.
                if (j .eq. 31) not_avail = .true.
                if (j .eq. 5) not_avail = .true.
                if (j .eq. 6) not_avail = .true.
                if (j .eq. 7) not_avail = .true.
                if (j .eq. 8) not_avail = .true.
                if (j .eq. 21) not_avail = .true.
                if (j .eq. 22) not_avail = .true.
                if (j .eq. 23) not_avail = .true.
                if (j .eq. 24) not_avail = .true.
              endif 
   	      If (j .Ge. 1 .And. j .Le. menu_num) Then
		 If (nfields .Lt. 4 .and. (.not. not_avail)) Then
		    If (j .Le. off) Then
		       status = SMG$Change_Rendition ( vdid, j+3, 4,
     .						1, 2, rend )
                       If (status .Ne. SS$_Normal) Then
	         	  FIT_Display_Plot_Menu = status
	         	  Call LIB$Signal(%Val(status))
	          	  Return
	               End If
		    Else
		       k = j - off
		       status = SMG$Change_Rendition ( vdid,  k+3,
     .						46, 1, 2, rend )
                       If (status .Ne. SS$_Normal) Then
	         	  FIT_Display_Plot_Menu = status
	         	  Call LIB$Signal(%Val(status))
	          	  Return
	               End If
		    End If

c
c	Check if element already selected.
c
		    select = .True.
		    Do i=1,nfields
		       If (menu_page(i) .Eq. page .And.
     .			   menu_element(i) .Eq. j) Then
			  select = .False.
		       End If
		    End Do

c
c	If not, save it...
c
		    If (select) Then
		       nfields = nfields + 1
		       menu_page(nfields) = page
		       menu_element(nfields) = j
		    Else
		       status = SMG$Ring_Bell ( vdid, 3 )
                       If (status .Ne. SS$_Normal) Then
	                  FIT_Display_Plot_Menu = status
	                  Call LIB$Signal(%Val(status))
	                  Return
	               End If
		    End If

		 Else
		    status = SMG$Ring_Bell ( vdid, 3 )
                    If (status .Ne. SS$_Normal) Then
	               FIT_Display_Plot_Menu = status
	               Call LIB$Signal(%Val(status))
	               Return
	            End If
		 End If

	      Else
		 status = SMG$Ring_Bell ( vdid, 3 )
                 If (status .Ne. SS$_Normal) Then
	            FIT_Display_Plot_Menu = status
	            Call LIB$Signal(%Val(status))
	            Return
	         End If
	      End If

	   End If

	End Do

c
c Clean up and exit.
c
	status = SMG$Erase_PasteBoard ( pbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Keyboard ( kbid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If

	status = SMG$Delete_Virtual_Display ( vdid )
	If (status .Ne. SS$_Normal) Then
	   FIT_Display_Plot_Menu = status
	   Call LIB$Signal(%Val(status))
	   Return
	End If


	FIT_Display_Plot_Menu = %Loc(FIT_Normal)

	Return
	End
