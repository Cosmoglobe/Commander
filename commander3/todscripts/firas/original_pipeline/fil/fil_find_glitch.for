      integer*4 function fil_find_glitch(aifg,noise,glitch_pos,glitch_amp)
c-------------------------------------------------------------------------------
c
c     Purpose: Find the largest peak in the ifg and return its position and its
c              amplitude scaled by a factor dependent on the ratio of the 
c              amplitude to the noise.
c
c     Author: R. Isaacman, GSC, 1/91
c
c     Input: aifg               r*4(512)  Real-valued interferogram.
c            noise              r*4  Noise of the interferogram.
c
c     Output: glitch_pos         i*4  Peak (glitch) position.
c             glitch_amp         r*4  Weighted peak (glitch) amplitude.
c
c     Modifications: K. Jensen, Hughes STX, 12/94 : New threshold and gain 
c                    values.
c                    K. Jensen, HSTX, 8/24/95 : Thresholds re-defined, based
c                    on LLSS minical analysis.
c                    K. Jensen, HSTX, 9/25/95, SPR 12256 : Thresholds 
c                    re-defined for positive/negative symmetry.
c
c-------------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fil_normal
c
c     Input parameters.
c
      real*4 aifg(512),noise
c
c     Output parameters.
c
      integer*4 glitch_pos
      real*4 glitch_amp
c
c     Local variables.
c
      integer*4 ifg_pos
      real*4 thresh

      fil_find_glitch = %loc(fil_normal)

      glitch_amp = 0.0
c
c     Find the highest peak of the interferogram.
c
      do ifg_pos = 1,512
         if (abs(aifg(ifg_pos)) .gt. abs(glitch_amp)) then
            glitch_amp = aifg(ifg_pos)
            glitch_pos = ifg_pos
         end if
      end do
c
c     If the peak amplitude exceeds the threshold, set a new
c     glitch amplitude. Positive peak thresholds are larger.
c
      thresh = glitch_amp/noise
      if ((thresh .ge. 5.5) .or. (thresh .le. -5.5)) then
         glitch_amp = glitch_amp * 0.2
      else if ((thresh .ge. 3.7) .or. (thresh .le. -3.7)) then
         glitch_amp = glitch_amp * 0.7
      else
         glitch_amp = 0.0
      end if

      return
      end
