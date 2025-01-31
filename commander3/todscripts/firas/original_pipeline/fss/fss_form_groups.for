      integer*4 function fss_form_groups(sky_buff, tot_recs, chan,
     .                                   num_groups, bad_groups, good_ifgs,
     .                                   bad_ifgs, bins, num_bins, dbins,
     .                                   num_dbins, temptol, minimum,
     .                                   pix_array, mode_array, bin_array, 
     .                                   dbin_array, size_array, split_array)

C------------------------------------------------------------------------
C    PURPOSE: Form groups from the records in the sky array. 
C             Set transition flags for groups. Pass numbers for 
C             coadd groups, number of ifgs used in coadd groups, 
C             and number of ifgs not in a coaddable group. 
C
C    AUTHOR: D. Bouler, STX, Jan, 1990
C
C    INVOCATION:   STATUS = fss_form_groups(sky_buff, tot_recs, chan,
C                                           num_groups, bad_groups, good_ifgs,
C                                           bad_ifgs, bins, num_bins, dbins,
C                                           num_dbins, temptol, minimum, 
C                                           pix_array, mode_array, bin_array, 
C                                           dbin_array, size_array, split_array)
C
C
C    INPUT PARAMETERS:     ARRAY sky_buff          short sci records 
C                          I*4   tot_recs          number of records in array
C                          R*4   bins(*)           ICAL Bin temps 
C                          I*4   num_bins          number of ICAL bins
C                          R*4   dbins(*)          dihedral Bin temps 
C                          I*4   num_dbins         number of dihedral bins
C                          R*4   temptol           ICAL Temperature tolerance
C                          I*4   minimum(4)        min number in group
C                          I*4   chan              current channel
C
C    OUTPUT PARAMETERS:    ARRAY sky_buff            grouped array 
C                          I*4   num_groups          number of coadd groups
C                          I*4   bad_groups          coadd groups < minimum
C                          I*4   good_ifgs           Number ifgs in coadd groups
C                          I*4   bad_ifgs            Number ifgs not in coadds
C                          I*4   pix_array (max_groups) pixel numbers for coadds
C                          I*4   mode_array(max_groups) modes for coadds
C                          I*4   size_array(max_groups) size of coadds
C                          I*4   bin_array(max_groups)  ICAL temp bins of coadds
C                          I*4   dbin_array(max_groups) dihedral temp bins of coadds
C                          I*2   split_array(max_groups) split group flags
C
C    SUBROUTINES CALLED:    None.
C
C    COMMON VARIABLES USED: None.
C 
C    INCLUDE FILES:         fss_include
C
c
c----------------------------------------------------------------------
C Changes:
C      SPR 9539,  FSS report is incorrect for one record passed/one coadd group
C         H. Wang, Hughes STX, March 4, 1992
C
C      SPR 9579,  FSS to limit numbers of IFGS per coadd group
C         H. Wang, Hughes STX, March 11, 1992
C
C      10 February 1995, FSS_FORM_GROUPS must set a flag indicating if coadd
C         groups are split (because group was larger than fac_max_num (100))
C         and if so whether the group is the first, last, or a middle of a
C         split.  This info is used for neighbor acceptability.
C         0 - group not split             1 - 1st split group
C         2 - middle split groups         3 - last split group
C
C      SPR 12197 - Modifications to implement /DBINS qualifier in place
C                  of /MAX_DIHED,  K.Jensen, HSTX, 23-May-1995
C
C----------------------------------------------------------------------
C                         PDL for FSS_FORM_GROUPS
c                         D. Bouler, Jan, 1990
c
c  initialize total number of coadd groups to one
c  initialize number of bad   coadd groups to zero
c  initialize number of good  coadd groups to zero
c  initialize pixel, mode, and temp arrays to zero
c  initialize group size array for first group to one
c  initialize new_group flag to false
c
c  set group pixel     to pixel of first record
c  set group ical temp to max of temp bin of first record 
c
c  do i = 2, tot_recs
c
c     pixel = sky_rec(i).pixel_no
c     temp  = sky_rec(i).ical_temp
c     mode  = 2 * sky_rec(i).mtm_speed + sky_rec(i).mtm_length
c     set new_group flag
c
c     if (not new group) then
c         increment group size array by 1
c     endif
c
c     if (new group or (last record)) then
c
c         pix_array (num_groups) = grp_pixel
c         mode_array(num_groups) = grp_mode
c         bin_array (num_groups) = bptr
c
c         if (new_group .and. (i .le. tot_recs)) then
c             set record to mark = i - size_array(num_groups)
c         elseif ((.not. new_group) .and. (i .eq. tot_recs)) then
c             set record to mark = i - size_array(num_groups) + 1
c         endif
c
c          if (number in old group .ge. minimum(chan)) then
c              set transition flag of first record in old group to 2
c              increment number of coadd groups by 1
c              increment number of good ifgs by number of ifgs in old group
c          else
c              set transition flag of first record in old group to 1
c              increment number of bad ifgs by number of ifgs in old group
c          endif
c
c          if (i .lt. tot_recs) num_groups = num_groups + 1
c
c          if ((new group) .and. (last record)) then
c               set transition flag of last record to 1
c               increment number of coadd groups
c               increment number of bad ifgs by 1
c               set size of last group to one
c          endif
c
c      endif
c    
c      if (not last record) then
c          reset number in group to 1
c          reset group pixel, scan mode, temp.
c      endif
c
c  enddo
c
c  fss_form_groups
c
c  return
c  end (pdl)
C----------------------------------------------------------------------

      implicit none
c
c Include files
c
      include   '(fss_include)'
      include   '(fut_params)'
C
C Variable declarations
C
      integer*4 chan                     ! channel to process
      integer*4 minimum(4)               ! Min number in group 
      integer*4 num_bins                 ! number of ICAL temp bins
      integer*4 num_dbins                ! number of dihedral temp bins
      integer*4 grp_pixel                ! current pixel number for group
      integer*4 grp_mode                 ! current scan mode for group
      integer*4 pixel                    ! pixel for current record 
      integer*4 speed                    ! mtm speed
      integer*4 length                   ! mtm length
      integer*4 mode                     ! scan mode for current record 
      integer*4 tot_recs                 ! number of sky_buff records
      integer*4 num_groups               ! number of coadd groups 
      integer*4 bad_groups               ! coadd groups < minimum
      integer*4 pix_array (fss_max_groups) ! array of pixels
      integer*4 mode_array(fss_max_groups) ! array of modes
      integer*4 size_array(fss_max_groups) ! array of group size
      integer*4 bin_array (fss_max_groups) ! Array of ICAL temp bin for each group
      integer*4 dbin_array (fss_max_groups)! Array of dihedral temp bin for each group
      integer*2 split_array(fss_max_groups) ! Array of split group flags
      integer*4 good_ifgs                ! number of ifgs in coadd groups
      integer*4 bad_ifgs                 ! number of ifgs not in coadd groups
      integer*4 i,j,k,l                  ! loop counters
      integer*4 grp_size                 ! Number of records in current group
      integer*4 status                   ! general return status
      integer*4 bptr                     ! Pointer to current ICAL temp bin
      integer*4 dbptr                    ! Pointer to current dihedral temp bin
      integer*4 mark_rec                 ! Record to mark transition flag
      integer*4 temp_num, mark_temp,tot_num,last_num, num_div,ix
      real*4    bins(*)                  ! ICAL Temp bins
      real*4    dbins(*)                 ! Dihedral Temp bins
      real*4    grp_temp                 ! ICAL Temp max for current group
      real*4    grp_dtemp                ! dihedral Temp max for current group
      real*4    temp                     ! ical temp for current record 
      real*4    dtemp                    ! dihedral temp for current record 
      real*4    temptol                  ! ICAL Temp tolerance

      logical*4 in_bin                   ! Is temp in current bin?
      logical*4 new_group                ! Is record start of new group ?
c
c Declare externals
c
      external fss_normal
c
c Dictionary and record declarations
c                                             
      dictionary 'fss_sssky'
      record /fss_sssky/ sky_buff(fss_max_buff)
c*
c******  Begin code *******************************************
c*
      fss_form_groups = %loc(fss_normal)         
      status          = %loc(fss_normal)         
c
c Initialize variables.
c
      num_groups = 1
      bad_groups = 0
      good_ifgs  = 0
      bad_ifgs   = 0
c      do i = 1, fss_max_groups
      do i = 1, 160
         pix_array(i)  = 0
         mode_array(i) = 0
         bin_array(i)  = 0
         dbin_array(i)  = 0
         size_array(i) = 0
         split_array(i) = 0
      enddo
      size_array(1) = 1
      new_group     = .false.
c
c Set parameters for current group equal to the first record.
c Use mode = (2 * length) + (speed); values are 0,1,2,3 (SS,SF,LS,LF).
c      
      bptr   = 1
      temp   = sky_buff(1).ical_temp
      in_bin = (abs(temp - bins(bptr)) .le. temptol)

      do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
         bptr = bptr + 1
         in_bin = (abs(temp - bins(bptr)) .le. temptol)
      enddo

      dbptr   = 1
      dtemp   = sky_buff(1).dihedral_temp
      in_bin = (dtemp .le. dbins(dbptr))

      do while ( (.not. in_bin) .and. (dbptr .lt. num_dbins))
         dbptr = dbptr + 1
         in_bin = (dtemp .le. dbins(dbptr))
      enddo

      grp_pixel  = sky_buff(1).pixel_no
      length     = 2 * sky_buff(1).mtm_scan_length  
      speed      = sky_buff(1).mtm_scan_speed
      grp_mode   = length + speed
      grp_temp   = bins(bptr) + temptol
      grp_dtemp  = dbins(dbptr)
c
c Search through records, setting transition flag when old group ends.
c Group break will occur when (1) Pixel number changes,
c                             (2) Scan mode changes,
c                             (3) ical temp crosses temp bin boundary.
c                             (4) dihedral temp crosses temp bin boundary.
c
               if (tot_recs .eq. 1) then
                   num_groups             =  1
                   mode    = speed + length
                   temp    = sky_buff(num_groups).ical_temp
                   dtemp   = sky_buff(num_groups).dihedral_temp
                   pix_array (num_groups) =grp_pixel
                   mode_array(num_groups) = mode
                   size_array(num_groups) = 1
                   split_array(num_groups) = 0
                   sky_buff(num_groups).transition_flag = 1
                   bad_groups = bad_groups + 1
                   bad_ifgs   = bad_ifgs   + 1
             
                   bptr = 1
                   in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
                      bptr = bptr + 1
                      in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   enddo
                   bin_array (num_groups) = bptr

                   dbptr = 1
                   in_bin = (dtemp .le. dbins(dbptr))
                   do while ( (.not. in_bin) .and. (dbptr .lt. num_dbins))
                      dbptr = dbptr + 1
                      in_bin = (dtemp .le. dbins(dbptr))
                   enddo

                   dbin_array (num_groups) = dbptr

       Else
        do i  = 2, tot_recs
          
         length  = 2 * sky_buff(i).mtm_scan_length
         speed   = sky_buff(i).mtm_scan_speed 
         mode    = speed + length
         temp    = sky_buff(i).ical_temp
         dtemp   = sky_buff(i).dihedral_temp
         pixel   = sky_buff(i).pixel_no

         new_group = (pixel .ne. grp_pixel) .or. (mode .ne. grp_mode)  .or. 
     .               (temp  .gt. grp_temp) .or. (dtemp .gt. grp_dtemp)  

         if ((.not. new_group) .and. (i .le. tot_recs)) then
              size_array(num_groups) = size_array(num_groups) + 1
         endif
c
c We found the beginning of a new group; increment the number of groups.
c If group is coaddable, set transition flag of old group to two
c                        and increment number of good ifgs.
c If group is not coaddable, set transition flag of old group to one
c                        and increment number of bad ifgs.
c Save the group pixel number, mode, size, and mean temp.
c If the last record is a group by itself, pick up that group too.
c
         if (new_group .or. (i .eq. tot_recs)) then
              pix_array (num_groups) = grp_pixel
              mode_array(num_groups) = grp_mode
              bin_array (num_groups) = bptr
              dbin_array (num_groups) = dbptr

              if (new_group .and. (i .le. tot_recs)) then
                  mark_rec = i - size_array(num_groups)
              elseif ((.not. new_group) .and. (i .eq. tot_recs)) then
                  mark_rec = i - size_array(num_groups) + 1
              endif

C Set a flag indicating if coadd groups are split (because group was larger
C than fac_max_num (100)) and if so whether the group is the first, last, or a
C middle of a split.  This info is used for neighbor acceptability.
C    0 - group not split             1 - 1st split group
C    2 - middle split groups         3 - last split group

              if (size_array(num_groups) .ge. minimum(chan)) then
                  sky_buff(mark_rec).transition_flag = 2
                  if (size_array(num_groups) .gt. fac_max_num ) then
                    tot_num = size_array(num_groups)
                    num_div= (size_array(num_groups)- 1)/fac_max_num
                    temp_num = size_array(num_groups)/(num_div + 1)
                    size_array(num_groups) = temp_num
                    good_ifgs = good_ifgs + size_array(num_groups)
                    split_array(num_groups) = 1
                    do ix =1, num_div-1
                        num_groups = num_groups + 1
                        mark_temp = mark_rec + (ix)*temp_num
                        size_array(num_groups) = temp_num
                        pix_array (num_groups) = grp_pixel
                        mode_array(num_groups) = grp_mode
                        bin_array (num_groups) = bptr
                        dbin_array (num_groups) = dbptr
                        sky_buff(mark_temp).transition_flag = 2
                        good_ifgs = good_ifgs + size_array(num_groups)
                        split_array(num_groups) = 2
                     enddo
                        last_num = tot_num - num_div*temp_num 
                        mark_temp = mark_rec + num_div *temp_num
                        num_groups = num_groups + 1
                        size_array(num_groups) = last_num
                        good_ifgs = good_ifgs + size_array(num_groups)
                        pix_array (num_groups) = grp_pixel
                        mode_array(num_groups) = grp_mode
                        bin_array (num_groups) = bptr
                        dbin_array (num_groups) = dbptr
                        sky_buff(mark_temp).transition_flag = 2
                        split_array(num_groups) = 3
                   else
                     split_array(num_groups) = 0
                     good_ifgs = good_ifgs + size_array(num_groups)
                   endif 
               else
                  sky_buff(mark_rec).transition_flag = 1
                  bad_groups = bad_groups + 1
                  bad_ifgs   = bad_ifgs   + size_array(num_groups)
                  split_array(num_groups) = 0
               endif

               if (i .lt. tot_recs) num_groups = num_groups + 1

               if (new_group .and. (i .eq. tot_recs)) then
                   num_groups             = num_groups + 1
                   pix_array (num_groups) = pixel
                   mode_array(num_groups) = mode
                   size_array(num_groups) = 1
                   sky_buff(i).transition_flag = 1
                   bad_groups = bad_groups + 1
                   bad_ifgs   = bad_ifgs   + 1
             
                   bptr = 1
                   in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
                      bptr = bptr + 1
                      in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   enddo

                   bin_array (num_groups) = bptr

                   dbptr = 1
                   in_bin = (dtemp .le. dbins(dbptr))
                   do while ( (.not. in_bin) .and. (dbptr .lt. num_dbins))
                      dbptr = dbptr + 1
                      in_bin = (dtemp .le. dbins(dbptr))
                   enddo

                   dbin_array (num_groups) = dbptr

                endif   !    (new_group .and. (i .eq. (tot_recs))
c
c  Find the temperature bins to start on and reset the group parameters.
c
               if (i .lt. tot_recs) then

                   bptr   = 1
                   in_bin = (abs(temp - bins(bptr)) .le. temptol)

                   do while ( (.not. in_bin) .and. (bptr .lt. num_bins))
                      bptr = bptr + 1
                      in_bin = (abs(temp - bins(bptr)) .le. temptol)
                   enddo

                   dbptr   = 1
                   in_bin = (dtemp .le. dbins(dbptr))

                   do while ( (.not. in_bin) .and. (dbptr .lt. num_dbins))
                      dbptr = dbptr + 1
                      in_bin = (dtemp .le. dbins(dbptr))
                   enddo

                   size_array(num_groups) = 1
                   grp_pixel = sky_buff(i).pixel_no
                   length    = 2 * sky_buff(i).mtm_scan_length 
                   speed     = sky_buff(i).mtm_scan_speed
                   grp_mode  = length + speed
                   grp_temp  = bins(bptr) + temptol
                   grp_dtemp = dbins(dbptr)
                endif  

           endif    !   (pixel .ne. grp_pixel) ...
      enddo         !   i = 2, tot_recs
      Endif  
      fss_form_groups = status

      return
      end
