      integer*4 function fic_find_glitch(aifg,noise,glitch_pos,glitch_amp)
c-------------------------------------------------------------------------------
c
c     Purpose: Find the largest peak in the ifg and return its position and its
c              amplitude scaled by a factor dependent on the ratio of the 
c              amplitude to the noise. If the peak amplitude is smaller than 
c              3.5 times the noise in the interferogram then it is assigned a 
c              value of zero.
c
c     Author: R. Isaacman, GSC, 1/91
c
c     Input: aifg               r*4(512)  Real-valued interferogram.
c            noise              r*4  Noise of the interferogram.
c
c     Output: glitch_pos         i*4  Peak (glitch) position.
c             glitch_amp         r*4  Weighted peak (glitch) amplitude.
c
c     Modifications:
c
c     SER 7985, S. Alexander, STX, 9/1/91: Implement new fic requirements.
c
c-------------------------------------------------------------------------------

      implicit none
c
c     Return statuses.
c
      external fic_normal
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

      fic_find_glitch = %loc(fic_normal)

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
c     Compare the peak amplitude to the 3.5 noise threshold.
c
      thresh = abs(glitch_amp/noise)
      if (thresh .ge. 5.0) then
         glitch_amp = glitch_amp * 0.2
      else if (thresh .ge. 4.0) then
         glitch_amp = glitch_amp * 0.6
      else if (thresh .ge. 3.5) then
         glitch_amp = glitch_amp * 0.8
      else
         glitch_amp = 0.0
      end if

      return
      end
