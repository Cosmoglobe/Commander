	integer * 4 function  fsl_read_reference (fakeit, upmode)

c-------------------------------------------------------------------------------
c
c	Function FSL_READ_REFERENCE
c
c       This function reads the appropriate apodization function for the scan 
c       mode (1 - 6) from the FEX_APODL reference dataset.  It then reads the 
c       appropriate electronics transfer function from the FEX_ETFL reference 
c       dataset. All data sets to be calibrated by FSL from the 4 channel and 
c       6 scan mode combinations have an apodization function and an 
c       electronics transfer function in the reference data sets.
c
c       The peak position in the FIL coadd record is for the non-linearized 
c       ifg. Requirements for FSL have changed from using the peak position
c       for the linearized ifg to the non-linearized ifg. Thus the peak 
c       position is now taken from the FIL coadd record.
c
c	Author:  
c                FCF_Read_Reference
c                Gene Eplee
c		 General Sciences Corp.
c		 9 October 1992
c
c                FSL_Read_Reference
c                Shirley M. Read
c                Hughes STX Corporation
c                August 1995
c
c-------------------------------------------------------------------------------
c
c	Input:
c		fakeit		integer * 4		fakeit flag
c		upmode		integer * 4		microprocessor mode
c
c	Output:
c		none
c
c	Subroutines called:
c		dfftrf
c		fut_apod_recnuml  
c		fut_get_recnum
c		lib$movc5
c		lib$signal
c
c	Include files:
c		fsl_config.txt
c		fsl_invoc.txt
c		fut_params.txt
c
c-------------------------------------------------------------------------------
c
c	Changes for FCF:
c
c	Modifications to recover low frequency short fast data.
c	Gene Eplee, GSC, 25 October 1993
c	SER 11395
c
c       Changes for FSL:
c
c       Shirley M. Read, Hughes STX Corporation, August 7, 1995 
c       Modified FCF_Read_Reference to FSL_Read_Reference for the new FIRAS 
c       pipeline which will process long spectra to get improved frequency 
c       resolution.
c           1. Changed status, include files, and function names for FSL.
c           2. Removed the FFT of the apodization function and its application
c              to the variances. The spectral variances now include 
c              apodization which was performed in FIL.
c           3. Removed all special processing for low frequency channel FS
c              and FL scan modes. The processsing is already performed in FIL.
c              All 6 scan modes now have corresponding apodization functions 
c              and electronics transfer functions. The FS and FL adds per 
c              group must be reset to 2 for the call to FUT_Get_Recnum.
c           4. Removed the call to FUT_Default_Peak. Due to changed requirements
c              for FSL to use the non-linearized peak position, the peak 
c              position can now be taken from the FIL coadd record.
c           5. Increased the array lengths for FFT and related functions.
c           6. Called the new FUT routine for apodization record number.
c
c-------------------------------------------------------------------------------

	implicit none

	include '(fut_params)'
	include '(fsl_config)'
	include '(fsl_invoc)'

	integer * 4	arecno		!  apodization function record number
	integer * 4	erecno		!  etf record number
	integer * 4	fakeit		!  fakeit flag
	integer * 4	io_stat		!  I/O return status
	integer * 4	j		!  a counter
	integer * 4	k		!  a counter
	integer * 4	linearized/0/	!  linearization flag: 0=non-linearized
                                        !                      1=linearized
	integer * 4	rstatus		!  return status
	integer * 4	status		!  return status
	integer * 4	upmode 		!  microprocessor mode
	integer * 4     scan_length     !  scan length (0=short, 1=long)
	integer * 4     ngroup          !  adds per group

	integer * 4	fut_apod_recnuml

	external	fsl_normal
	external	fsl_readapod
	external	fsl_readetf

C
C  Initialize function status.
C
	status = %loc(fsl_normal)
C
C  Read the appropriate apodization function.
C

c
c  Get the record number. Call the new fut_apod_recnuml which handles all
c  6 scan modes. The calling sequence is different from the previous function.
c
	rstatus = fut_apod_recnuml (fcc_chan, fcc_smode, fakeit, upmode,
     .				   fcc_ngroup, linearized, arecno)

c
c  Read the apodization function into a real*8 512 point array. The array is 
c  stored in the FSL_Config include file.
c
	call lib$movc5 (0,,0,4096,apod_fcn)
	read (dir_lun(1)arecno, iostat=io_stat) apod_fcn

	if (io_stat .ne. 0) then
	   status = %loc(fsl_readapod)
	   call lib$signal (fsl_readapod, %val(2), %val(arecno), %val(io_stat))
	endif

	if (status .eq. %loc(fsl_normal)) then
C
C  Read the appropriate electronics transfer function.
C

c
c  Get the record number. For FS and FL scan modes (5 and 6), the adds per 
c  group in the call to FUT_Get_Recnum must be set to 2.
c
	   if ((fcc_smode .eq. 5) .or. (fcc_smode .eq. 6)) then
	       ngroup = 2	  
	       call fut_get_recnum (fakeit, fcc_speed, fcc_chan, upmode, 
     .                              ngroup, erecno)
	   else
	       call fut_get_recnum (fakeit, fcc_speed, fcc_chan, upmode, 
     .                              fcc_ngroup, erecno)
	   endif
c
c  Read the ETF. The ztrans complex array is stored in the FSL_Config 
c  include file.
c
	   read (dir_lun(2)erecno, iostat=io_stat) ztrans

	   if (io_stat .ne. 0) then
	      status = %loc(fsl_readetf)
	      call lib$signal (fsl_readetf, %val(1), %val(erecno),
     .					    %val(io_stat))
	   endif

	endif	!  status from apodization function read


	fsl_read_reference = status

	return
	end
