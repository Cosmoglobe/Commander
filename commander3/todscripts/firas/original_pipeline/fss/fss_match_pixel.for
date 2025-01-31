      integer*4  function  fss_match_pixel  (pixel, sky_buff, num_recs)
C---------------------------------------------------------------------------
C    PURPOSE: Search through the sorted and merged records to find the
C             first occurence of a record with a matching pixel.
C
C    AUTHOR: L. Rosen, HSTX, 3 February 1995
C
C    INVOCATION:  status = fss_match_pixel (pixel, sky_buff, num_recs)
C
C    INPUT PARAMETERS:       I*4   pixel        The pixel to match
C                            array sky_buff     The records to be written
C                            I*4   num_recs      The number of records
C
C    OUTPUT PARAMETERS:      None.  Return value of function is index of
C                            first record with matching pixel number.  If
C                            pixel not found return 0.
C
C    Method:  Binary search records to get a record with matching pixel.  Then
C             work backwards to find the first record with that pixel. 
C-----------------------------------------------------------------------------
      implicit none
c
c Include files.
c
      include       '(fss_include)'

c Passed Parameters.

      integer*4      pixel                 ! Pixel number to match
      integer*4      num_recs              ! Number of records to write
c
c Dictionary and record declarations.
c
      dictionary 'fss_sssky'
      record     /fss_sssky/ sky_buff (fss_max_buff)
c
c Local variables.
c
      integer*4      pix               ! pixel value at current index
      integer*4      index             ! current index being checked
      integer*4      delta
      integer*4      low, high         ! indices of range ends during search
c
c******  Begin code *******************************************
c
      fss_match_pixel = 0

      low = 1
      high = num_recs
      delta = int ((high - low) / 2)
      index = low + delta
      pix = sky_buff (index).pixel_no
      do while (pix .ne. pixel .and. delta .gt. 0)
         if (pixel .lt. pix) then
            high = index
         else                      ! pixel .gt. pix
            low = index
         endif
         delta = int ((high - low) / 2)
         index = low + delta
         pix = sky_buff (index).pixel_no
      enddo
      if (pix .ne. pixel .and. high .ne. low) then  ! try index = high
         index = high
         pix = sky_buff (index).pixel_no
      endif
      if (pix .ne. pixel) then
         fss_match_pixel = 0
      else

c Found a matching pixel.  Now work back to first occurence.

         do while (pix .eq. pixel .and. index .gt. 1)
            index = index - 1
            pix = sky_buff (index).pixel_no
         enddo
         if (pix .ne. pixel) then
            fss_match_pixel = index + 1
         else
            fss_match_pixel = index
         endif
      endif
      return
      end
