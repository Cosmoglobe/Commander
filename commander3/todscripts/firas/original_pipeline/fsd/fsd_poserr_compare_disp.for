C-------------------------------------------------------------------------------

	Integer*4 Function FSD_Poserr_Compare_Disp ( 
	1	  Smooth_Ifgs, Template_Ifg, Sigmas )

C-------------------------------------------------------------------------------
C
C	Purpose:
C
C       To compare smoothed IFGs to a template IFG and plot template IFG,
C	difference between individual IFG and template, and derivative of 
C       template.
C
C	Author: Shirley M. Read
C		STX, July 1988
C
C	Invocation:
C		Status = FSD_Poserr_Compare_Disp ( Smooth_Ifgs, 
C			 Template_Ifg, Sigmas )
C
C	Modificaton History:
C
C	  Author	    Date	  Modification
C	  ----------------------------------------------------------------------
C	R. Kummerer	Aug  8, 1988	SPR 2061, Calculate and display IFG
C					noise.
C	R. Kummerer	Sep 27, 1989	SPR 4635, Correct "error handling"
C					such that POSERR issues an abort
C					message only when an error occurs.
C
C	S. Alexander	Mar  7, 1990	SER 5728, Add capability to pass
C					in a PLT command file.
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name	       Dimension  Type	  Description
C 	  ----------------------------------------------------------------------
C	  Smooth_Ifgs  (512,Num)  R*4     Smoothed IFGS
C	  Template_Ifg (512)      R*4     Template IFG
C	  Sigmas       (Num)      R*4     RMS noise from IFGs - template 
C	
C	Output Parameters:
C	  Name		Type		Description
C 	  ----------------------------------------------------------------------
C	
C	Subroutines Called:
C         
C
C	Common Variables Used:
C
C	  Name		Type	Use	Description
C	  ----------------------------------------------------------------------
C
C	Include Files:
C		Fut_Params.txt
C	        Fsd_Poserr.txt
C
C------------------------------------------------------------------------------

	Implicit None

C	Include files.

	Include '(fut_params)'
	Include '(fsd_poserr)'

C	Passed Parameters.

	Real*4 smooth_ifgs(512,num)
	Real*4 template_ifg(512)
	Real*4 sigmas(num)
	Real*4 ensemble_sigma(512)	! deviation of ifg from ensemble

C	Local Declarations.

	real*4 deriv(512)	! Derivative of the template at each point
	real*4 diff(512)        ! Difference between the IFG and the template
	real*4 y(512,5)         ! Vectors to be plotted
	real*4 ymin             ! Minimum Y value for plot
	real*4 ymax             ! Maximum Y value for plot
	integer*4 mx1 / 512 /	! Actual first dimension of Y array
	integer*4 npts / 512 /  ! Number of points to plot
	integer*4 nveci / 5 /   ! Number of vectors to be plotted	
	integer*4 iery(5) /0,0,0,0,0/ ! Flags for plot errors:
				    ! -1 = plot error as sqrt y
				    !  0 = no errors
				    ! +1 = explicit errors
	integer*4 ier               ! Error flag
	integer*4 ncmd              ! Number of commands

	integer*4 status        ! Processing Status
	integer*4 success / 1 /, error / 2 /
	integer*4 ifg_num	! Selected IFG number in series
	real*4 sum		! summation variable
	real*4 ensemble_avg	! IFG average over freq channel
	integer*4 i,j,k		! Indices

	character*80 cmd(20)    ! Plot commands
	data cmd(1) /  !Outer top text line for plot
     .  'LA OT TEMPLATE, IFG - TEMPLATE, TEMPLATE DERIVATIVE, IFG NOISE' /
	data cmd(3) / 'LA X POSITION' /
	data cmd(4) / 'LA Y IFG DATA IN COUNTS' /
	data cmd(5) / 'PLOT VERTICAL' /
	character*5 label_t/ 'LA T ' /
	character*4 lifg / 'IFG ' /
	character*6 lgmt / ' GMT= ' /
	character*7 lchan / ' CHAN= ' /
	character*12 lscan / ' SCAN MODE= ' /
	character*3 lnum      ! Number of IFG in char

	external    fsd_normal
	external    fsd_abort
C
C	Set processing status to success.
C
	status = success
C
C	Compute the derivative of the template IFG at each point.
C	If three points of the IFG are fit to a parabola and the derivative 
C	is computed at the center point Xi (the points on the x axis are equally
C	spaced), the derivative may be computed by the formula:
C		 DY/DX = 2 a Xi + b = ( Yi+1 - Yi-1 ) / ( Xi+1 - Xi-1 )
C		       in which ( Xi+1 - Xi-1 ) is a constant value of 2.
C
	do i = 2, 511

	  deriv(i) = ( template_ifg(i+1) - template_ifg(i-1)) / 2

	enddo
C
C	Compute the ensemble sigmas.
C
	do i = 1, 512

	  sum = 0.
	  do j = 1, num
	    sum = sum + smooth_ifgs(i,j)
	  end do

	  ensemble_avg = sum / num

	  sum = 0.
	  do j = 1, num
	    sum = sum + (smooth_ifgs(i,j) - ensemble_avg) ** 2
	  end do

	  ensemble_sigma(i) = sqrt(sum/(num-1))

	end do

C	Estimate the slope at the two end points.

	deriv(1) = deriv(2)
	deriv(512) = deriv(511) 
C
C	Load the plot vectors with values constant to all plots.

	call lib$movc3(512*4,x(1),y(1,1))	! X values
	call lib$movc3(512*4,template_ifg(1),y(1,2))	! Y1 values
	call lib$movc3(512*4,deriv(1),y(1,4))	! Y4 values
	call lib$movc3(512*4,ensemble_sigma(1),y(1,5))	! Y5 values

C	Write commands for min and max values on each template Y plot.	

	do i = 1, 4
	  if (i .ne. 2) then
	    j = i + 1
	    ymin = y(1,j)
	    ymax = y(1,j)
	    do k = 2, 512
	      if ( y(k,j) .lt. ymin ) ymin = y(k,j)
	      if ( y(k,j) .gt. ymax ) ymax = y(k,j)
	    enddo
	    write ( cmd(i+5), '(a3,i1,2(1x,e14.4))')
     .          'R Y', j, ymin, ymax
	  end if
        enddo

	ncmd = 9

	print *, 'Display template, IFG - template and derivative:'
C
	ifg_num = 1

	do while ( ifg_num .ne. 0 )
	  type 30
30	  format(' Enter serial number of IFG (0 to quit) > ', $)
	  accept *, ifg_num
C
	  if ((ifg_num.gt.0).and. (ifg_num.le.NUM)) then

C	Compute the difference between the individual IFG and the template IFG.

	     do i = 1, 512
		diff(i) = smooth_ifgs(i,ifg_num) - template_ifg(i)
	     enddo

             call lib$movc3(512*4,diff(1),y(1,3))	! Y2 values

C	Set up the top label.

	     write ( lnum,'(i3)') ifg_num
 	     cmd(2) = label_t // lifg // lnum // lgmt // 
     .        gmt(ifg_num) // lchan // fac_channel_ids(chan_id) // lscan// scan

C	Set the min and max of the differences.

	     ymin = y(1,3)
	     ymax = y(1,3)
	     do k = 2, 512
	       if ( y(k,3) .lt. ymin ) ymin = y(k,3)
	       if ( y(k,3) .gt. ymax ) ymax = y(k,3)
	     enddo
	     write ( cmd(7), '(a3,i1,2(1x,e14.4))')
     .          'R Y', 3, ymin, ymax

	     if (plt_com .eq. fac_present) then
	        ncmd = ncmd + 1
	        cmd(ncmd) = '@' // plt_com_file
	     end if

	     call PLT(y,iery,mx1,npts,nveci,cmd,ncmd,ier)

	  elseif (ifg_num .ne. 0) then

	     print *, 'Improper choise of spectum number. Try again.'

	  endif ! Ifg_num
	enddo ! Ifg_num

C	Set return status according to processing. Since no functions are
C	invoked for Build 4.1, the Fut_Params will be used for status.

	if ( status .eq. success ) then
	  fsd_poserr_compare_disp = %Loc(FSD_Normal)
	else
	  fsd_poserr_compare_disp = %Loc(FSD_Abort)
	endif

	return
	end
