      integer*4  function  fss_match_group  (pixel,  rec_num, sky_buff,
     &                                       num_recs, indices, num_match,
     &                                       bins, num_bins, dbins, 
     &                                       num_dbins, temptol)
C---------------------------------------------------------------------------
C    PURPOSE: Search through the sorted and merged records to find indices
C             of records that match scan mode, ical bin, and dihedral bin
C             for the specified pixel.
C
C    AUTHOR: L. Rosen, HSTX, 7 February 1995
C
C    INVOCATION:  status = fss_match_group (pixel, rec_num, sky_buff, num_recs,
C                                           indices, num_match, bins, num_bins,
C                                           dbins, num_dbins, temptol)
C
C    INPUT PARAMETERS:       I*4   pixel        The pixel to match
C                            I*4   rec_num      Record number of record to match
C                            array sky_buff     The records to be written
C                            I*4   num_recs     The number of records
C                            R*4   bins(100)    ICAL temperature bins
C                            I*4   num_bins     Number of ICAL bins
C                            R*4   dbins(100)   dihedral temperature bins
C                            I*4   num_dbins    Number of dihedral bins
C                            R*4   temptol      ICAL temp. tolerance in K.
C    OUTPUT PARAMETERS:
C                            I*4   indices(10)  Indices of matching records
C                            I*2   num_match    Number of matching records
C
C    FUNCTIONS CALLED:
C                            FSS_MATCH_PIXEL - return index of 1st record
C                            with matching pixel number.
C
C
C    SPR 12197 - Modifications to implement /DBINS qualifier in place
C                of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C    PDL:
C
C   Get index of first coadd group with matching pixel.  (Call FSS_MATCH_PIXEL)
C   Loop through short science records, starting at index.
C      If pixel doesn't match:
C         Exit record loop.
C      Else
C         If scan mode doesn't match:
C            Advance record counter to next record with non-zero transitionflag.
C         Else
C            If ICAL bin doesn't match:
C               Advance record counter to next rec with non-zero transitionflag.
C         Else
C            If Dihedral bin doesn't match:
C               Advance record counter to next rec with non-zero transitionflag.
C            Else
C               Store record counter in array of matching indices.
C            EndIf
C         EndIf
C      EndIf
C   End Loop
C   Return array of indices of matching records.
C   End.
C-----------------------------------------------------------------------------
      implicit none
c
c Include files.
c
      include       '(fss_include)'

c Passed Parameters.

      integer*4      pixel                 ! Pixel number to match
      integer*4      rec_num               ! Number of record to match
      integer*4      num_recs              ! Number of records
      real*4         bins(100)             !ICAL temperature bins
      integer*4      num_bins              !Number of ICAL bins
      real*4         dbins(100)            !Dihedral temperature bins
      integer*4      num_dbins             !Number of Dihedral bins
      real*4         temptol               !ICAL temp. tolerance in K.
      integer*4      indices(fss_max_groups) !Indices of matching records
      integer*2      num_match             !Number of matching records
c
c Dictionary and record declarations.
c
      dictionary 'fss_sssky'
      record     /fss_sssky/ sky_buff (fss_max_buff)
c
c Functions called.
c
      integer*4      fss_match_pixel
c
c Local variables.
c
      integer*4      zero (fss_max_groups) / fss_max_groups * 0 /
      integer*2      i
      integer*4      index             ! Index of record being checked
      logical*1      done              ! Flag indicating done searching
      logical*1      next              ! Flag to do a loop again
      integer*2      bin_count         ! ICAL bin index
      real*4         bin_max           ! ICAL bin temp + tolerance
      integer*2      bin1, bin2        ! bin to match and bin being checked
c
c******  Begin code *******************************************
c
c Initialize stuff to 0
c
      fss_match_group = 0
      call lib$movc3 (fss_max_groups*4, zero, indices)
      num_match = 0
      done = .false.

C Get index of first coadd group with matching pixel.  (Call FSS_MATCH_PIXEL)

      index = fss_match_pixel (pixel, sky_buff, num_recs)
      if (index .eq. 0) done = .true.

C Loop through short science records, starting at index.

      do while (index .le. num_recs .and. .not. done)

C If pixel doesnt match,  Exit record loop.

         if (sky_buff(index).pixel_no .ne. pixel) then
            done = .true.
         else

C Check scan mode.  If it doesnt match then advance record counter until the
C next record with non-zero transition flag.

            if (sky_buff(index).mtm_scan_speed .ne.
     &          sky_buff(rec_num).mtm_scan_speed .or.
     &          sky_buff(index).mtm_scan_length .ne.
     &          sky_buff(rec_num).mtm_scan_length) then

               if (index .eq. num_recs) then
                  done = .true.
               else
                  next = .true.
                  do while (index .lt. num_recs .and. next)
                     index = index + 1
                     next = sky_buff(index).transition_flag .eq. 0
                  enddo
               endif
            else

C Check ICAL bin.  This is harder to check since its values are floats and
C not exact integers.  First get ICAL bins by comparing ICAL temp with bin +
C temperature tolerance.  If ICAL bins dont match then advance record counter
C to next rec with non-zero transition flag.

               next = .true.
               bin_count = 1
               do while (bin_count .le. num_bins .and. next)
                  bin_max = bins(bin_count) + temptol
                  if (sky_buff(rec_num).ical_temp .le. bin_max) then
                     bin1 = bin_count
                     next = .false.                    ! ICAL bin found
                  else
                     bin_count = bin_count + 1
                  endif
               enddo
               next = .true.
               bin_count = 1
               do while (bin_count .le. num_bins .and. next)
                  bin_max = bins(bin_count) + temptol
                  if (sky_buff(index).ical_temp .le. bin_max) then
                     bin2 = bin_count
                     next = .false.                    ! ICAL bin found
                  else
                     bin_count = bin_count + 1
                  endif
               enddo
               if (bin1 .ne. bin2) then
                  if (index .eq. num_recs) then
                     done = .true.
                  else
                     next = .true.
                     do while (index .lt. num_recs .and. next)
                        index = index + 1
                        next = sky_buff(index).transition_flag .eq. 0
                     enddo
                  endif
               else

C Check Dihedral bin. If dihedral bins dont match then advance record counter
C to next rec with non-zero transition flag.

                  next = .true.
                  bin_count = 1
                  do while (bin_count .le. num_dbins .and. next)
                     bin_max = dbins(bin_count)
                     if (sky_buff(rec_num).dihedral_temp .le. bin_max) then
                        bin1 = bin_count
                        next = .false.                    ! Dihedral bin found
                     else
                        bin_count = bin_count + 1
                     endif
                  enddo
                  next = .true.
                  bin_count = 1
                  do while (bin_count .le. num_dbins .and. next)
                     bin_max = dbins(bin_count)
                     if (sky_buff(index).dihedral_temp .le. bin_max) then
                        bin2 = bin_count
                        next = .false.                    ! Dihedral bin found
                     else
                        bin_count = bin_count + 1
                     endif
                  enddo
                  if (bin1 .ne. bin2) then
                     if (index .eq. num_recs) then
                        done = .true.
                     else
                        next = .true.
                        do while (index .lt. num_recs .and. next)
                           index = index + 1
                           next = sky_buff(index).transition_flag .eq. 0
                        enddo
                     endif
                  else

C Found a matching coadd in the pixel.  Store record counter in array of
C matching group indices.  Then advance index to next coadd group.

                     next = .true.
                     do while (index .lt. num_recs .and. next)
                        num_match = num_match + 1
                        indices (num_match) = index
                        index = index + 1
                        next = sky_buff(index).transition_flag .eq. 0
                     enddo
                     if (index .eq. num_recs .and. next) then
                        num_match = num_match + 1
                        indices (num_match) = index
                        done = .true.
                     endif
                  endif   ! matching dihedral bin
               endif      ! matching ical bin
            endif         ! matching scan mode
         endif            ! matching pixel
      enddo
      return
      end
