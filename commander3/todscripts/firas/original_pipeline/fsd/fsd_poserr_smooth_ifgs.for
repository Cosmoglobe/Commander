C-------------------------------------------------------------------------------

	Integer*4 Function FSD_Poserr_Smooth_Ifgs ( Fgs, Np, Smooth_Ifgs )

C-------------------------------------------------------------------------------
C
C	Purpose:
C
C       To smooth the individual IFGs using a boxcar function as a low pass
C	filter and a user selected number of points for averaging at each
C	position in the IFG.
C
C	Author: Shirley M. Read
C		STX, July 1988
C
C	Invocation:
C		Status = FSD_Poserr_Smooth_Ifgs ( Fgs, Np, Smooth_Ifgs ) 
C
C	Modificaton History:
C
C	  Author	    Date	  Modification
C	  ----------------------------------------------------------------------
C	R. Kummerer	Sep 27, 1989	SPR 4635, Correct "error handling"
C					such that POSERR issues an abort
C					message only when an error occurs.
C
C	Input Files:
C
C	Output Files:
C
C	Input Parameters:
C	  Name	       Dimension  Type	  Description
C 	  ----------------------------------------------------------------------
C	  Fgs          (512,Num)  R*4     Raw IFGS
C	  Np                      I*4     Number of points for averaging
C	
C	Output Parameters:
C	  Name	       Dimension  Type	  Description
C 	  ----------------------------------------------------------------------
C	  Smooth_Ifgs  (512,Num)  R*4     Smoothed IFGS
C	
C	Subroutines Called:
C         
C
C	Common Variables Used:
C
C	  Name	       Dimension  Type	Use              Description
C	  ----------------------------------------------------------------------
C	  Num                     I*4   Array dimension  Number of IFS in set
C	  Nsw          (100)      I*2   Divisor          Number of mirror sweeps
C
C	Include Files:
C		Fut_Params.txt
C	        Fsd_Poserr.txt
C
C	Processing Method:
C	Set the return status to success.
C	Normalize each IFG to number of mirror sweeps.
C	Using a boxcar function with user selected number of points, smooth
C	each IFG and output the smoothed IFGs.
C	Set the return status to success or failure according to processing.
C
C------------------------------------------------------------------------------

	Implicit None

C	Include files.

	Include '(fut_params)'
	Include '(fsd_poserr)'

C	Passed Parameters.

	Real*4 Fgs(512,num)
	Integer*4 Np
	Real*4 Smooth_Ifgs(512,num)

C	Local Declarations.

	integer*4 status             ! Processing status
	integer*4 success / 1 /, error / 2 /
	real*4    norm_ifgs(512,100) ! IFGs divided by number of mirror sweeps
	integer*2 numsweep(100)      ! Number of mirror sweeps in each IFG
	real*4    acc_ifg            ! Accumulator for IFG values
	real*4    weight             ! Weight total for averaging
	integer*2 zero / 0 /
	integer*2 ix, jx, kx         ! Indices

C	Externals

	External  FSD_Normal
	External  FSD_Abort

C	Initialize processing to success.

	status = success

C 	Normalize IFGs to number of mirror sweeps.

	do ix = 1, num
	  if ( nsw(ix) .le. zero ) then
	    numsweep(ix) = 1
	  else
	    numsweep(ix) = nsw(ix)
	  endif
	enddo

	do ix = 1, num
	  do jx = 1, 512
	     norm_ifgs(jx,ix) = fgs(jx,ix) / numsweep(ix)
	  enddo
	enddo

c	Compute the averaged points of the smoothed IFGs using the 
c	boxcar function and the user selected number of points.

	weight = 2 * np + 1
	do ix = 1, num
	  do jx = np + 1, 512 - np
	    acc_ifg = 0.0
	    do kx = jx - np, jx + np
	      acc_ifg = acc_ifg + norm_ifgs(kx,ix)
	    enddo
	    smooth_ifgs(jx,ix) = acc_ifg / weight
	  enddo	  
	enddo

c	Estimate the smoothed points at the end points of the IFGs.

	do ix = 1, num
	  do jx = 1, np
	    weight = np + jx
	    acc_ifg = 0.0
	    do kx = 1, jx + np
	      acc_ifg = acc_ifg + norm_ifgs(kx,ix)
	    enddo
	    smooth_ifgs(jx,ix) = acc_ifg / weight
          enddo
	  do jx = 513 - np, 512
	    weight = np + 513 - jx
	    acc_ifg = 0.0
	    do kx = jx - np, 512
	      acc_ifg = acc_ifg + norm_ifgs(kx,ix)
	    enddo
	    smooth_ifgs(jx,ix) = acc_ifg / weight
          enddo
	enddo

c	Set return status to success or failure. For Build 4.1 there are no
c	functions invoked by FSD_Poserr_Smooth_Ifgs and no condition handler
c	for FSD_Poserr. Thus the parameters in FUT_Params are used for the
c	return status. Future enhancements may require code modifications.

	if (status .eq. success) then
	  fsd_poserr_smooth_ifgs = %Loc(FSD_Normal)
	else
	  fsd_poserr_smooth_ifgs = %Loc(FSD_Abort)
	endif	  

	return
	end
