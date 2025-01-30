      integer*4 function frd_gal_cut_flv(pixel,galexc_val)
c-------------------------------------------------------------------------------
c
c     Purpose: Return the first pixel number greater than or equal to the input
c              pixel number whose galactic center is outside of the input 
c              galactic latitude cutoff on either side of the galactic plane.
c
c     Author: S. Brodd, HSTX, 10/95.
c
c     Input: pixel         integer*4         Input pixel number.
c            galexc_val    real*4            Galactic latitude exclusion.
c
c     Output: pixel         integer*4         Output pixel number.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Include files.
c
      include '(fut_params)'
c
c     Return statuses.
c
      external frd_normal,frd_eof
c
c     Input parameters.
c
      real*4 galexc_val
c
c     Input/output parameters.
c
      integer*4 pixel
c
c     Local variables.
c
      real*4 glat,evec(3),gvec(3)
c
c     Determine the galactic latitude of the input pixel number.
c
      call firas_cenpix(pixel,evec)
      call xcc_e_to_g(evec,fac_epoch,gvec)
      glat = atan2d(gvec(3),sqrt(gvec(1)**2 + gvec(2)**2))
c
c     Find the first pixel outside the excluded galactic latitude range.
c
      do while ((abs(glat) .lt. galexc_val) .and. (pixel .lt. 6143))
         pixel = pixel + 1
         call firas_cenpix(pixel,evec)
         call xcc_e_to_g(evec,fac_epoch,gvec)
         glat = atan2d(gvec(3),sqrt(gvec(1)**2 + gvec(2)**2))
      end do
c
c     Check the validity of the output pixel number.
c
      if ((pixel .lt. 6143)  .or.
     .    ((pixel .eq. 6143) .and. (abs(glat) .ge. galexc_val))) then
         frd_gal_cut_flv = %loc(frd_normal)
      else
         frd_gal_cut_flv = %loc(frd_eof)
      end if

      return
      end
