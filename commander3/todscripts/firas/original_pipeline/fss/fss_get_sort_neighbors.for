      integer*4  function  fss_get_sort_neighbors  (rec_num, sky_buff, num_recs,
     &                                              bins, num_bins, dbins,
     &                                              num_dbins, temptol, 
     &                                              split, size, neighbs,
     &                                              num_neighbs)

C------------------------------------------------------------------------
C    PURPOSE: Get array of indices of neighbor pixel records with matching
C       scan mode, ICAL bin, and dihedral bin, sorted in order of increasing
C       absolute difference of galactic latitude between the record and the
C       mean of the group for which we are finding neighbors.
C
C    AUTHOR: L. Rosen, HSTX, 13 February 1995
C
C    INVOCATION:  status = fss_get_sort_neighbors  (rec_num, sky_buff, num_recs,
C    &                     bins, num_bins, dbins, num_dbins, temptol,
C    &                     split, size, neighbs, num_neighbs)
C
C    INPUT PARAMETERS:
C     integer*4      rec_num                   ! Record number of 1st of group
C     record/fss_sssky/ sky_buff(fss_max_buff) ! Array of short sky recs
C     integer*4      num_recs                  ! Number of records in ss array
C     real*4         bins(100)                 ! ICAL temperature bins
C     integer*4      num_bins                  ! Number of ICAL bins
C     real*4         dbins(100)                ! Dihedral temperature bins
C     integer*4      num_dbins                 ! Number of dihedral bins
C     real*4         temptol                   ! ICAL temp. tolerance in K.
C     integer*2      split                     ! Coadd group split flag
C     integer*4      size                      ! Size of coadd group
C    OUTPUT PARAMETERS:
C     integer*4      neighbs(fss_max_groups)   !Indices sorted neighbor records
C     integer*2      num_neighbs               ! Number of neighbor records
C 
C    SUBROUTINES CALLED:
C                    fss_avg_gallat
C                    fss_match_group
C                    upx_8_neighbors
C                    sor$begin_sort
C                    sor$release_rec
C                    sor$merge_sort
C                    sor$return_rec
C                    sor$end_sort
C                    lib$movc3
C
C    INCLUDE FILES:  fss_include, '($SSDef)'
C
C
C    SPR 12197 - Modifications to implement /DBINS qualifier in place
C                of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FSS_GET_SORT_NEIGHBORS PDL:
C
C   Calculate the average galactic latitude of the main coadd group.
C                                                            (FSS_AVG_GALLAT)
C   Get neighboring pixel numbers.                     (Call UPX_8_NEIGHBORS)
C   For each neighbor pixel:
C      Find coadd group(s) with the same pixel,             (FSS_MATCH_GROUP)
C         scan speed, scan length, ICAL bin, and dihedral bin.
C      If there are any matching groups:
C         If the main coadd group had been split:
C            Store record indices and absolute galactic latitude differences
C               from the matching groups that are within the time range.  See
C               FSS_FORM_GROUPS for definition of array of split flags.
C         Else
C            Store record indices and absolute galactic latitude differences
C               from the matching groups into array of neighbor record indices.
C         EndIf
C      End If there are any matching groups.
C   End For each neighbor pixel.
C   If there are any indices in the array of neighbor record indices:
C      Use SOR$ routines to sort indices of neighboring records in increasing
C         order of absolute galactic latitude difference from the average
C         latitude of the main group.
C   End If
C   Return (array of sorted indices of neighbors)
C   End.
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      include       '($ssdef)'
      include       '(fss_include)'
c
c Function declarations
c
      integer*4      fss_avg_gallat
      integer*4      fss_match_group
      logical*1      time_le
      integer*4      sor$begin_sort
      integer*4      sor$release_rec
      integer*4      sor$return_rec
      integer*4      sor$end_sort
      integer*4      sor$sort_merge
c
c Passed parameters
c
      integer*4      rec_num                   ! Record number of 1st of group
      dictionary 'fss_sssky'
      record/fss_sssky/ sky_buff(fss_max_buff) ! Array of short sky recs
      integer*4      num_recs                  ! Number of records in ss array
      real*4         bins(100)                 ! ICAL temperature bins
      integer*4      num_bins                  ! Number of ICAL bins
      real*4         dbins(100)                ! dihedral temperature bins
      integer*4      num_dbins                 ! Number of dihedral bins
      real*4         temptol                   ! ICAL temp. tolerance in K.
      integer*2      split                     ! Coadd group split flag
      integer*4      size                      ! Size of coadd group
      integer*4      neighbs(fss_max_groups)   !Indices sorted neighbor records
      integer*2      num_neighbs               ! Number of neighbor records
C
C Local variables
C
      integer*4      status                ! General return status
      real*4         avg_gallat            ! Average galactic latitude
      integer*4      pixel                 ! Pixel number
      integer*4      last_pixel / -9999 /  ! Pixel number of last group
      integer*4      neighb_pix (8)        ! Pixel numbers of neighbors
      integer*4      num                   ! Number of neighbors (7 or 8)
      integer*2      pix                   ! Neighbor pixel counter
      integer*2      rec                   ! Matching record counter
      integer*2      i                     ! counter
      integer*4      indices(fss_max_groups) ! Indices of matching records
      integer*2      num_match             ! Number of matched recs in neighbor
      structure /ind_lat_st /
         integer*4   index                 ! Index of record in sky_buff array
         real*4      gallatdif             ! Abs difference of galactic latitude
      end structure
      integer*4      irec_size             ! Size of ind_lat_st in bytes
      parameter      (irec_size = 8) 
      record /ind_lat_st/ ind_lat_recs (fss_max_groups)  ! Array to sort
      character*8    one_rec               ! "Record" for sorting ind_lat_recs
      integer*4      last                  ! Record number of last in group
      integer*4      last_adt (2)          ! Time of last record of group
      integer*4      first_adt (2)         ! Time of first record of group
      integer*2      sort_keys(5)          ! Provide information on sorting 
c
      save last_pixel                      ! Save last pixel number
      save neighb_pix
      save num
c
c Declare externals
c
      external       fss_normal
      external       DSC$K_DType_F               !SOR type definition (real*4)
c*
c******  Begin code *******************************************
c*
      fss_get_sort_neighbors = %loc(fss_normal)         

C Calculate the record number of the last of the coadd group

      last = rec_num + size - 1

C Calculate the average galactic latitude of the main coadd group.

      status = fss_avg_gallat (sky_buff, rec_num, num_recs, avg_gallat)

c Now find neighbors.  If pixel number is the same as it was for the last
c coadd group, then we save time by not having to get the neighbor pixel
c numbers again.

      pixel = sky_buff (rec_num).pixel_no
      if (pixel .ne. last_pixel) then
         last_pixel = pixel
         do pix=1,8
            neighb_pix(pix) = -1
         enddo
         call upx_8_neighbors (pixel, 6, neighb_pix, num)
      endif                                  ! check whether pixel = last_pixel

C Find coadd groups in neighbor pixels with matching scan mode, ical bin,
C and dihedral bin.
C Use function FSS_Match_Group to find indices of all records of groups
C that match (ususally 0 or 1 group).  If none are found, num_match = 0.

      num_neighbs = 0
      do pix=1, num
         status = fss_match_group (neighb_pix(pix), rec_num, sky_buff, num_recs,
     &                             indices, num_match, bins, num_bins, dbins, 
     &                             num_dbins, temptol)

C If the main coadd group is not split then:
C    Store record indices and absolute galactic latitude differences
C    from the matching groups into array of neighbor record indices.
C Else
C    Store record indices and absolute galactic latitude differences
C    from the matching groups that are within the time range.  Split value:
C         0 - group not split             1 - 1st split group
C         2 - middle split groups         3 - last split group
C    If first split, only use records with time <= end time of group,
C    If last split, only use records with time >= end time of group,
C    If middle split, only use records with time in between.
C EndIf

         if (split .eq. 0) then
            do rec = 1, num_match
               ind_lat_recs (num_neighbs + rec).index = indices(rec)
               ind_lat_recs (num_neighbs + rec).gallatdif = abs ( real (
     &            sky_buff(indices(rec)).galactic_latitude) - avg_gallat )
            enddo
            num_neighbs = num_neighbs + num_match
         elseif (split .eq. 1) then
            last_adt(1) = sky_buff(last).time(1)
            last_adt(2) = sky_buff(last).time(2)
            do rec = 1, num_match
               if (time_le (sky_buff(indices(rec)).time, last_adt)) then
                  num_neighbs = num_neighbs + 1
                  ind_lat_recs (num_neighbs).index = indices(rec)
                  ind_lat_recs (num_neighbs).gallatdif = abs ( real (
     &               sky_buff(indices(rec)).galactic_latitude) - avg_gallat )
               endif
            enddo
         elseif (split .eq. 3) then
            first_adt(1) = sky_buff(rec_num).time(1)
            first_adt(2) = sky_buff(rec_num).time(2)
            do rec = 1, num_match
               if (time_le (first_adt, sky_buff(indices(rec)).time)) then
                  num_neighbs = num_neighbs + 1
                  ind_lat_recs (num_neighbs).index = indices(rec)
                  ind_lat_recs (num_neighbs).gallatdif = abs ( real (
     &               sky_buff(indices(rec)).galactic_latitude) - avg_gallat )
               endif
            enddo
         else
            first_adt(1) = sky_buff(rec_num).time(1)
            first_adt(2) = sky_buff(rec_num).time(2)
            last_adt(1) = sky_buff(last).time(1)
            last_adt(2) = sky_buff(last).time(2)
            do rec = 1, num_match
               if (time_le (first_adt, sky_buff(indices(rec)).time) .and.
     &             time_le (sky_buff(indices(rec)).time, last_adt)) then
                  num_neighbs = num_neighbs + 1
                  ind_lat_recs (num_neighbs).index = indices(rec)
                  ind_lat_recs (num_neighbs).gallatdif = abs ( real (
     &               sky_buff(indices(rec)).galactic_latitude) - avg_gallat )
               endif
            enddo
         endif
      enddo                                      ! do for each neighbor pixel

C   If there are any indices in the array of neighbor record indices:
C      Use SOR$ routines to sort indices of neighboring records in increasing
C         order of absolute galactic latitude difference from the average
C         latitude of the main group.
C   End If

      if (num_neighbs .gt. 0) then
c*
c*    Intialize Sort_Keys for the SOR$ routines
c*
         Sort_Keys(1) = 1                     !  Total of one key
         Sort_Keys(2) = %loc(DSC$K_Dtype_F)   !  Key 1 = abs diff gal lat
         Sort_Keys(3) = 0                     !  Ascending sort
         Sort_Keys(4) = 4                     !  Start byte of key
         Sort_Keys(5) = 4                     !  Length in bytes of key
c
c Initialize the sort process
c
         status = SOR$Begin_Sort (Sort_Keys, irec_size)
         if (status .ne. ss$_normal) then
            call lib$signal(%val(status))
         else
c
c Release records of sky_buff array to sort process
c
            rec = 1
            do while (status .eq. ss$_normal .and. rec .le. num_neighbs)
               call lib$movc3 (irec_size, ind_lat_recs(rec), %ref(one_rec))
               status = SOR$Release_Rec ( one_rec )
               if (.not. status) call lib$signal(%val(status))
               rec = rec + 1
            enddo
c
c Perform sort
c
            if (status .eq. ss$_normal) then
               status = SOR$Sort_Merge ()
               if (status .ne. ss$_normal) then
                  call lib$signal(%val(status))
               endif
            endif
c
c Extract sorted records
c    lib$movc3 will move first 4 bytes (index) into neighbs array
c
            rec = 1
            do while (status .eq. ss$_normal .and. rec .le. num_neighbs)
               status = SOR$Return_Rec ( one_rec )
               if (status .ne. ss$_normal) then
                  call lib$signal(%val(status))
               else
                  call lib$movc3 (4, %ref(one_rec), neighbs(rec))
                  rec = rec + 1
               endif
            enddo
c
c End the sort process
c
            if (status .eq. ss$_normal) then
               status = SOR$End_Sort ()
               if (status .ne. ss$_normal) call lib$signal(%val(status))
            endif
         endif
      endif
      if (status .ne. ss$_normal .and. num_neighbs .gt. 0) then
         fss_get_sort_neighbors = status
      endif
      return
      end
