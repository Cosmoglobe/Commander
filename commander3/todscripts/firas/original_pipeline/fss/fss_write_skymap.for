      integer*4 function fss_write_skymap (lu_output, sky_buff, num_recs,
     &                                     bins, num_bins, dbins, num_dbins, 
     &                                     temptol, split_array, size_array)

C------------------------------------------------------------------------
C    PURPOSE: Write the sorted and merged records to an output skymap.
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION:  status = fss_write_skymaps (lu_output, sky_buff, num_recs,
C                          bins, num_bins, dbins, num_dbins, temptol, 
C                          split_array, size_array)
C
C    INPUT PARAMETERS:       I*4   lu_output    The output skymap LU
C                            array sky_buff     The records to be written
C                            I*4   num_recs     The number of records to write
C                            R*4   bins(*)      ICAL temperature bins
C                            I*4   num_bins     Number of ICAL bins
C                            R*4   dbins(*)     dihedral temperature bins
C                            I*4   num_dbins    Number of dihedral bins
C                            R*4   temptol      ICAL temp. tolerance in K.
C                            I*2   split_array()  Split flags of coadd groups
C                            I*4   size_array()  Size of coadd groups
C
C    OUTPUT PARAMETERS:      None.  
C 
C    SUBROUTINES CALLED:     csa_field_offset_values
C                            csa_write_pixels
C                            fss_get_neighbors
C                            lib$movc3
C
C    COMMON VARIABLES USED:  None.
C 
C    INCLUDE FILES:          fss_include
C
C----------------------------------------------------------------------
C Changes:
C    1/30/95  L. Rosen.  FSS must now write out coadds from neighboring
C             pixels of the same scan mode and ical bin, for each coadd.
C             The first transition flag of the neighbor coadds will be
C             set to 3.  This will allow FIC to use neighbor ifgs to
C             form the coadd template for deglitching.  FSS_WRITE_SKYMAP
C             now loops through coadds, writing to output, instead of
C             doing just the one write of the entire skymap.
C    New PDL:
C
C   If there are records to write:
C      Call CSA to set offset values.
C   Loop through records:                   (loop always begins new coadd group)
C      Write the first record of the coadd group.        (with csa_write_pixels)
C      While in the group:                     (next record transition flag = 0)
C         Write the record.                              (with csa_write_pixels)
C      Get array of indices of sorted neighbor records. (FSS_GET_SORT_NEIGHBORS)
C      If there are neighbors:
C         Copy the first record and set its transition flag to 3.
C         Write the copied/modified record.           (with csa_write_pixels)
C         For each remaining neighbor record:
C            Set the neighbor record pixel number = the main coadd group pixel
C            Set the transition flat to 0.
C            Write the copied/modified record.     (with csa_write_pixels)
C         End For each neighbor record
C      End If there are any neighbors.
C   End loop over records.
C   End.
C
C   SPR 12197 - Modifications to implement /DBINS qualifier in place
C               of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      include       '(fss_include)'
c
c Function declarations
c
      integer*4      csa_field_offset_values
      integer*4      csa_write_pixels
      integer*4      fss_get_sort_neighbors

c Passed parameters

      integer*4      num_recs              !Number of records to write
      real*4         bins(*)               !ICAL temperature bins
      integer*4      num_bins              !Number of ICAL bins
      real*4         dbins(*)              !dihedral temperature bins
      integer*4      num_dbins             !Number of dihedral bins
      real*4         temptol               !ICAL temp. tolerance in K.
      integer*4      lu_output             !Output skymap LU
      integer*2      split_array(fss_max_groups) ! Array of split group flags
      integer*4      size_array(fss_max_groups) ! Array of sizes of groups
c
c Dictionary and record declarations
c                                             
      dictionary 'fss_sssky'
      record     /fss_sssky/ sky_buff(fss_max_buff), rec_copy

c local variables

      integer*4      status                !General return status
      integer*4      i,j,k,n               !Loop counters
      integer*4      pixel                 !Pixel number
      integer*4      rec_num               !Record number of 1st of group
      logical*1      more                  !More records to write
      integer*4      group_num             !Coadd group number
      integer*4      neighbs(fss_max_groups)  !Indices sorted neighbor records
      integer*2      num_neighbs           !Number of neighbor records
c
c Declare externals
c
      external       fss_normal
      external       fss_noskyrecs
      external       fss_csawriterr
      external       csa_normal
c*
c******  Begin code *******************************************
c*
      fss_write_skymap = %loc(fss_normal)         

      if (num_recs .ne. 0) then

          status = csa_field_offset_values(pix_offset,time_offset,
     .                                     end_offset,lu_output)
          if (.not. status) then
               call lib$signal(%val(status))
          else

C Loop through coadd groups.
C    1st record is start of coadd group.  Save pixel and record number.
C    Write the record.  Copy records until group ends.  j is record counter.

             group_num = 0
             j = 1
             do while (j .le. num_recs .and. status)
                pixel = sky_buff(j).pixel_no
                rec_num = j
                group_num = group_num + 1
                status = csa_write_pixels (lu_output, sky_buff(j),1, blk_cnt)
                more = j .lt. num_recs         ! true if more records to write
                if (more) then
                   do while (status .eq. %loc (csa_normal) .and. more)
                      j = j + 1
                      if (sky_buff(j).transition_flag .eq. 0) then
                         status = csa_write_pixels (lu_output, sky_buff(j), 1,
     &                                              blk_cnt)
                         more = j .lt. num_recs
                      else
                         more = .false.
                      endif
                   enddo
                else
                   j = j + 1                             ! Ends main the loop
                endif
                if (status .ne. %loc (csa_normal)) then
                   call lib$signal(fss_csawriterr, %val(1), %val(status))
                else

c Finished writing coadd group.  Get array of indices of sorted neighbor
c records.  Call FSS_GET_SORT_NEIGHBORS.

                   status = fss_get_sort_neighbors (rec_num, sky_buff, num_recs,
     &                      bins, num_bins, dbins, num_dbins, temptol,
     &                      split_array(group_num), size_array(group_num),
     &                      neighbs, num_neighbs)

C If there are neighbors:
C    Copy the first record and set its transition flag to 3.
C    Write the copied/modified record.           (with csa_write_pixels)

                   if (num_neighbs .gt. 0) then
                      call lib$movc3 (fss_byte_recsize, sky_buff(neighbs(1)),
     &                                rec_copy)
                      rec_copy.transition_flag = 3
                      rec_copy.pixel_no = pixel
                      status = csa_write_pixels (lu_output, rec_copy,1,blk_cnt)
                      n = 1
                      more = n .lt. num_neighbs                 ! true if more
                      do while (status .eq. %loc(csa_normal).and. more)
                         n = n + 1
                         call lib$movc3 (fss_byte_recsize,
     &                                   sky_buff(neighbs(n)), rec_copy)
                         rec_copy.pixel_no = pixel
                         rec_copy.transition_flag = 0
                         status = csa_write_pixels (lu_output, rec_copy,
     &                                              1, blk_cnt)
                         more = n .lt. num_neighbs
                      enddo
                      if (status .ne. %loc (csa_normal)) then
                         call lib$signal (fss_csawriterr, %val(1), %val(status))
                      endif
                   endif            ! If there are any neighbors
                endif               ! Status check on writing coadd group
             enddo                  ! for all short sky records
          endif                     ! if offset values set ok.
      else                          ! if there are records to write
          call lib$signal(fss_noskyrecs)
          status = %loc(fss_noskyrecs)
      endif
      fss_write_skymap = status
      return
      end
