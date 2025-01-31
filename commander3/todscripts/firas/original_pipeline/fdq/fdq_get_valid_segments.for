	INTEGER*4 FUNCTION FDQ_GET_VALID_SEGMENTS(MAX_ENTRIES, FILEX,
	1	DELFLG, NSEG, MORE_SEGMENTS, FILENAMES)
C/
C/	PROGRAM NAME:
C/	  FDQ_GET_VALID_SEGMENTS
C/
C/	PROGRAM DESCRIPTION:
C/	  This routine takes as input the filled arrays FILEX and DELFLG
C/	  (filled by FDQ_SAVE_CAT_INFO), NSEG (the no. of segments for 
C/	  each channel), and returns the next "valid" group of file segments
C/	  (one for each channel to be processed if it exists) and the 
C/	  corresponding HKP.  A "valid" file segment is one which has the
C/	  same file extension name as the HKP file and does not have the 
C/	  deletion flag set. The channel file must not have the run flag
C/	  set already. 
C/
C/	AUTHOR:
C/	  EDWIN H. FUNG
C/	  GSFC
C/	  APRIL 23, 1987
C/
C/	MODIFIED BY:
C/	  
C/	  Shirley M. Read
C/	  STX
C/	  January 12, 1988
C/	  REASON: 	Converted from subroutine to function for interface
C/			with Fut_Error condition handler. Added error checking
C/		        and calls to Lib$Signal.
C/
C/	  Shirley M. Read
C/	  STX
C/	  March 1988
C/	  REASON:       The new Cobetrieve query catalog routines can now be
C/	                used to obtain information about each dataset. 
C/	                The input arrays created from the catalog records
C/	                may be used to build the set of filenames for the
C/ 			new Cobetrieve FORTRAN user open. This routine was
C/	 		rewritten to make use of the catalog information.
C/	                
CH	CHANGE LOG:
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
C/
C/	CALLING SEQUENCE:
C/	  Return_Status =  FDQ_GET_VALID_SEGMENTS(MAX_ENTRIES, FILEX, DELFLG,
C/			   NSEG, MORE_SEGMENTS, FILENAMES)
C/
C/	INPUT PARAMETERS:
C/	  MAX_ENTRIES	I*4	Maximum number of segments for processing in one
C/				run of DATA QUALIFY.  This is defined as a parameter
C/				in calling routine FDQ_GET_NXT_SEGMENT, but since
C/				it is used in the dimension of CENOS and DDTS (which
C/				are passed input parameters), this value has to be
C/				passed to this routine also.
C/	  FILEX(max_entries,5)  C*20  File extension names for the segments.
C/	  DELFLG(max_entries,5) L*1   Deletion flags for each segment.
C/	  NSEG(5)	        I*2   Number of catalog entries to be processed
C/				      (max. is MAX_ENTRIES). Each dataset has 
C/				      its own separate count (that's why NSEG 
C/				      is an array of 5 elements).
C/	  MORE_SEGMENTS(4)      L*1   Flag (one for each channel) to indicate
C/				      whether there are any more catalog 
C/				      entries for a channel. This is also an 
C/				      output parameter.
C/
C/	OUTPUT PARAMETERS:
C/	  MORE_SEGMENTS(4)      L*1   Flag (one for each channel) to indicate 
C/			              whether there are any more catalog 
C/				      entries for a channel. This is also an 
C/				      input parameter.
C/	  FILENAMES(5)          C*39  Filenames for the group of corresponding
C/				      science channel and HKP files.
C/
C/	INPUT/OUTPUT FILES:
C/	  NONE
C/
C/	INCLUDE FILE USED:
C/	  NONE
C/
C/	SUBROUTINES CALLED:
C/	  NONE
C/
C/	METHOD USED:
C/	  < PDL > :
C/
C/	INITIALIZE filenames to blank.
C/	INITIALIZE found, deleted to false.
C/	DO for channel 1 to 4
C/	   IF chan_done(channel) THEN
C/	      SET more_segments(channel) to false
C/	   ENDIF
C/	ENDDO
C/	DO WHILE not found and current HKP segment < number HKP segments
C/	   INCREMENT current HKP segment
C/         SET the current HKP file extension to current name in FILEX array
C/	   DO for channel 1 to 4
C/	      IF more_segments(channel) THEN
C/	         DO for current channel segment to number channel segments
C/	            IF channel file extension matches HPK file extension THEN
C/		       RESET current channel segment 
C/	               IF channel segment run flag OK and not deleted THEN
C/		 	  SET found to true
C/			  STORE the channel filename in the FILENAMES array
C/	               ELSE
C/			  SET deleted flag
C/		       ENDIF
C/		       IF current channel segment equals last segment for 
C/			  channel THEN
C/		          SET chan_done(channel) flag to true
C/		       ENDIF
C/		    ENDIF
C/		 ENDDO
C/	      ENDIF
C/	   ENDDO
C/	   IF found THEN
C/	      STORE the HKP filename in the FILENAMES array
C/	   ELSE
C/	      SIGNAL information : deleted
C/	   ENDIF
C/	ENDDO
C/	IF not found and last HKP file has been checked for match THEN
C/	   SET more_segments for all channels to false
C/	ENDIF
C/      RETURN with function status set
C/

	IMPLICIT	NONE

!	Passed Parameters

	integer*4	max_entries,i
	character*20    filex(max_entries,5)	! File extension names
	logical*1       delflg(max_entries,5)   ! Deletion flags
	integer*2	nseg(5)		        ! Number of segments
	logical*1	more_segments(4)        ! Flag for more science files
	character*39    filenames(5)            ! Filenames for segment group

!	Include Files and External Parameters

	include		'($ssdef)'

	external   	FDQ_NORMAL
	external	FDQ_ABERR
	external        FDQ_SEGDELETE
	external        FDQ_FPPRUN
	external        FDQ_FDQRUN
	external        FDQ_NOSCIMATCH
	
!!!!!!!!!!!!!!!!!!!!!!!!!
!			!
!     Local storage     !
!			!
!!!!!!!!!!!!!!!!!!!!!!!!!

	integer*4 	retstat		! Return status
	integer*4 	success / 1 /, error / 2 /  ! Values for status
	integer*4       zero  / 0 /
	integer*2       one   / 1 /
	integer*2       two   / 2 /
	character*7     hkp_rdl / 'NFS_HKP' /   ! HKP RDL name
	character*10    sci_rdl(4) / 'FPP_SDF_RH', 'FPP_SDF_RL',   ! Sci RDLs
	1			     'FPP_SDF_LH', 'FPP_SDF_LL' /
	integer*4       fpp_cur_seg / 0 /	  ! Current FPP segment
	integer*4       HKP_cur_seg / 0 /	  ! Current HKP segment
	integer*4       sci_cur_seg(4) / 4 * 1 /  ! Current science segments
	character*20    Fpp_file_ext		  ! FPP file extension
	character*20    hkp_file_ext		  ! Housekeeping file extension
	character*28    hkp_filename 		  ! Housekeeping filename
	logical*1	found			  ! Matching science found
        logical*1       done/.false./
	logical*1       deleted			  ! Deleted flag
	logical*1	chan_done(4) / 4 * .false. /  ! Completed channels
	integer*2	ix, jx		          ! indices
	character*39    blank			  ! Blank characters
	character*1     blanks(39) / 39 * ' ' /
	equivalence     (blank, blanks(1))
        character*20    nblank
        Character*1     nblanks(20)/20*' '/
	equivalence     (nblank, nblanks(1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			   !
!     Code begins here     !
!			   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Set return status to success.

	retstat = success

!	Initialize filenames for current call.

	do ix = 1, 5
	  filenames(ix) = blank
	enddo
	
!	Initialize flags for current call.

	found = .false.
	deleted = .false.

!	Update more_segments for channels which are complete.

	do ix = 1, 4				! Update more_segments
	   if ( chan_done(ix) ) then
		more_segments(ix) = .false.
	   endif
	enddo

        If ((.not. more_segments(1)) .and.
	1   (.not. more_segments(2)) .and.
	1	(.not. more_segments(3)) .and.
	1	(.not. more_segments(4)) ) done = .true.

!	DO loop to get a next group of filenames for a new segment.
!	Continue looping until at least one matching science file is found
!	to go with the housekeeping data or the list of housekeeping file
!	extensions is exhausted.
        if (.not. done) then
	DO WHILE ( (.not. found) .and. (FPP_cur_seg .lt. nseg(1)) )

	  FPP_cur_seg = FPP_cur_seg + 1
	  HKP_cur_seg = HKP_cur_seg + 1
          Do i=1,4
           If (filex(FPP_cur_seg,i) .ne. nblank) then 
    	     FPP_file_ext = filex(FPP_cur_seg,i)	! FPP file extension
           Endif
          Enddo
	  HKP_file_ext = filex(HKP_cur_seg,5)	! HKP file extension

	  do ix = 1, 4
	    if ( more_segments(ix) ) then
	      do jx = sci_cur_seg(ix), nseg(ix)
	        if ( filex(jx,ix) .eq. FPP_file_ext ) then
		  sci_cur_seg(ix) = jx		! Update sci current segment
	          if ( delflg(jx,ix) ) then	! Delete flag set in catalog
		    deleted = .true.
	          else				! Successful match found
		    found = .true.
	            filenames(ix) = sci_rdl(ix) // '.' // filex(jx,ix)
	          endif
	          if ( sci_cur_seg(ix) .eq. nseg(ix)) chan_done(ix) = .true.
                endif
	      enddo				! sci_cur_seg to nseg
	    endif				! more_segments(ix)
	  enddo					! channel 1 to 4
	  hkp_filename = hkp_rdl // '.' // hkp_file_ext
	  if ( found ) then			! Store HKP filename
	    filenames(5) = hkp_filename
	  elseif (deleted) then
	    call lib$signal( FDQ_SEGDELETE, %val(1), hkp_filename )
	    deleted = .false.
C	  else
C	    call lib$signal( FDQ_NOSCIMATCH, %val(1), hkp_filename)
	  endif
	ENDDO			! while not found and hkp_cur_seg < nseg(5)

!	If all HKP processed and no science to match, reset flags.

	If ((.not. found) .and. (nseg(5) .eq. hkp_cur_seg )) then
	  do ix = 1, 4
	    more_segments(ix) = .false.
	  enddo
         Endif
!	Set function return status
        endif 
	if (retstat.eq.success) then
	  FDQ_GET_VALID_SEGMENTS = %loc(FDQ_NORMAL)
	else
	  FDQ_GET_VALID_SEGMENTS = %loc(FDQ_ABERR)
	endif

	return
	end
