	INTEGER*4 FUNCTION FDQ_SAVE_CAT_INFO ( ID, NENTRIES, 
	1	MAX_ENTRIES, CATALOG_REC, NSEG, FILEX, DELFLG )
C/
C/	PROGRAM NAME:
C/	  FDQ_SAVE_CAT_INFO
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine takes the information in the Catalog_Rec gotten by 
C/	  COBETRIEVE query catalog routines and stores it into the input
C/	  storage arrays for use in grouping the dataset files for the new
C/	  Fortran user open in COBETRIEVE. The routine extracts and stores
C/	  the file extension names and the deletion flags for each file,
C/	  calls FDQ_Check_Entry to see if the segment has been modified and
C/	  stores the value of the modify flag bit for each segment. If the
C/	  status is OK, the number of segments for this dataset is stored 
C/	  in NSEG. If the number of segments exceeds the maximum allowed 
C/	  number of entries, a message is signalled but processing continues.
C/
C/	AUTHOR:
C/	  Edwin H. Fung
C/	  GSFC
C/	  April 22, 1987
C/
C/	MODIFIED BY:
C/
C/	  Shirley M. Read
C/	  STX
C/	  January 15, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  March 1988
C/        REASON:       Function had to be rewritten to save information
C/	                required for the new Fortran User Open in COBETRIEVE.
C/
C/	CALLING SEQUENCE:
C/	  Return Status = FDQ_SAVE_CAT_INFO (ID, NENTRIES, MAX_ENTRIES, 
C/		          CATALOG_REC, NSEG, FILEX, DELFLG )
C/
C/	INPUT PARAMETERS:
C/	  ID	        I*2     Dataset index ID: Channels 1 to 4, HKP 5
C/	  NENTRIES	I*4	Number of entries found from CT_CAT_INFO
C/	  MAX_ENTRIES	I*4	Maximum no. of entries that can be processed
C/				by a single run of DATA QUALIFY.  This 
C/				affects the array sizes of catalog records.
C/	  CATALOG_REC	B*1	Catalog entries as returned by COBETRIEVE query
C/
C/	OUTPUT PARAMETERS:
C/	  NSEG(5)	         I*2  Number of segments for each science or  
C/				      HKP dataset. Must be <= to MAX_ENTRIES.
C/	  FILEX(MAX_ENTRIES,5)   C*20 File extension names for each entry.
C/	                	      COBETRIEVE catalog record field.
C/	  DELFLG(MAX_ENTRIES,5)  L*1  Delete flag for each file entry.
C/		               
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILES USED:
C/	  TBD
C/
C/	SUBROUTINES CALLED:
C/
C/	METHOD USED:
C/	  Trivial
C/
C/---------------------------------------------------------------------------
CH CHANGE LOG:
CH	August 24, 1988, R. Kummerer, SPR 2409, Always passed first element
CH	of CATALOG_REC array to FDQ_CHECK_ENTRY, instead of current element
CH	being examined.
CH
CH      Version 4.1.1 10/20/88, SPR 2665, Shirley M. Read, STX
CH    		FDQ needs the L3 FDQ_Save_Cat_Info to get common file versions.
CH	   	Processing of FDQ requires that the housekeeping file and
CH		four science files have the same file extension and version
CH	        numbers. This constraint occured when the new COBETRIEVE 
CH	        allowed multiple version numbers for the same data.
CH
CH	Version 4.2.1 02/09/89, SPR 2700, Shirley M. Read, STX
CH		FDQ falls over for out of order and bad time ranges.The FIRAS
CH	 	Stripper now writes the science records using the telemetry
CH		minor frame time at transmit for the primary key in order to
CH		ensure a valid time tag. A new FIRAS preprocessor reads the raw
CH		science records, computes the midpoint of collect time and flags
CH		any records with bad midpoint times. FDQ must now access the
CH		computed collect time from the new field in the science record
CH	        in order to group the four science channels with the bracketing
CH		housekeeping frames and set the data quality summary flag when
Ch		the badtime flag is set. The bad record is written back to the
CH		archive immediately so that the next record may be considered
CH	        for the group. The average midpoint of collect times is still
CH		used as the time tag of the FDQ engineering record and is thus
CH		stored in the science record as the cross reference. However,
CH		the transmit times of the science records are stored in the 
CH		FDQ engineering record. FDQ must now enter the science gain into
CH		each science record from the bracketing housekeeping frames. The
CH		two LS bits must now be checked in the certification mask when
CH		running FDQ. The bit must be set for FPP but not for FDQ.
CH
CH	Version 4.4.1 07/22/89, SER 4168, R. Kummerer, STX
CH		There have been problems during test operations for FIRAS
CH		processing due to the required clean-up of the archives after
CH		an FPP or FDQ abort. The FPR tracking system compounds the
CH		problems. Files with non-matching version numbers seem often
CH		to result from improper clean-up. Bad record times cause
CH		SEGCTL to abort and mess up the tracking system. It was
CH		decided to change the modify of the science records in FPP
CH		and FDQ to a simple COBETRIEVE read of the existing records
CH		from a dataset and write a modifed dataset with the same
CH		information which was entered on the modify. Two new science
CH		data sets will be required: a science dataset of raw science
CH		data plus FPP input and a science dataset with FPP science
CH		data plus FDQ input. These datasets will be FPP_SDF_xx, where
CH		xx is the channel id (RH, RL, LH or LL) and FDQ_SDF_xx, where
CH		xx is the channel id. The new datasets must be opened and
CH		processed in FPP and FDQ. 
CH
CH	Version 4.4.1 08/21/89, SER 4210, R. Kummerer, STX
CH		Prevent overlaps in raw science segments related changes.
CH
C/---------------------------------------------------------------------------

	IMPLICIT	NONE

C	Passed Parameters

	integer*2       id		! Id for channel index or HKP index
	integer*4       nentries	! Number of catalog entries
	integer*4	max_entries     ! Maximum allowed number of cat entries
	integer*2	nseg(5)		! Number of segments for each dataset
	character*20    filex(max_entries,5)   ! Extension names of files 
					       ! to be processed
	logical*1	delflg(max_entries,5)  ! Delete flag for each file

	integer*4	iext		!number of characters in file extension

C	Include files and external parameters.

	dictionary 'CCM_CME_Catalog_Entry'	!Info retrieved from CT
	record     /CCM_CME_Catalog_Entry/Catalog_Rec(max_entries)

	external   	FDQ_NORMAL
	external	FDQ_ABERR
	external        FDQ_BADATIDX
	external	FDQ_CHECKENTRY
	external        FDQ_NUMCATINFO

C	Functions

C	Local Declarations

	integer*2       ntrunc		! Number of segments exceeding maximum
	integer*4	cur_info	! Index for catalog records
	integer*4 	retstat		! Return status
	integer*4 	success / 1 /, error / 2 /  ! Values for status
	integer*4       zero / 0 /      ! Zero value
	integer*4	status		! Dummy status variable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C	Set return status to success.

	retstat = success

C	Check for a valid dataset index to arrays.

	if ( id .lt. 1 .or. id .gt. 5 ) then
	  call lib$signal(FDQ_BADATIDX, %val(1), %val(id))
	  retstat = error
	endif

C	Initialize counter for current information file segments .

	cur_info = 1

	do while (cur_info .le. nentries .and. cur_info .le. max_entries)

C	 Check status before proceeding.

         if ( retstat .eq. success ) then

C        Move catalog record fields into storage arrays.

	    iext = index(catalog_rec(cur_info).filename_extension, ' ') - 1

	    filex(cur_info,id) = catalog_rec(cur_info).filename_extension(1:iext)

	    delflg(cur_info,id) = catalog_rec(cur_info).deletion_flag

	 endif		! Return status is success

C        Increment current info counter for loop.

	 cur_info = cur_info + 1

	enddo

C	Check status before proceeding. If OK, store the actual number of 
C	file segments for the dataset. Also check if the number of catalog 
C	records from CT query exceeded the maximum allowed. If so, signal 
C	a message and continue processing.

        if ( retstat .eq. success ) then
	  nseg(id) = cur_info - 1
	  ntrunc = nentries - (cur_info - 1)
	  if ( ntrunc .gt. 0 ) then
	      call Lib$Signal( FDQ_NUMCATINFO, %val(3),
	1	   %val(Max_Entries), %val(id), %val(ntrunc))
	  endif
	endif  	! Return status is success

C	Set function to return status

	if (retstat.eq.success) then
	  FDQ_SAVE_CAT_INFO = %loc(FDQ_NORMAL)
	else
	  FDQ_SAVE_CAT_INFO = %loc(FDQ_ABERR)
	endif

	return
	end
