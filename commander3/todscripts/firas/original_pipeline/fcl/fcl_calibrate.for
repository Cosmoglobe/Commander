      integer*4 function fcl_calibrate(gain,n,n2,eq,da)
c-------------------------------------------------------------------------------
c
c     Purpose: Apply calibration to four intermediate matrices.
c
c     Author: S. Brodd, HSTX, 12/95, SPR 12291
c
c     Input: gain  c*16(361)      Gain function
c            n     c*16(361,361)  Intermediate matrix N
c            n2    c*16(361,361)  Intermediate matrix N2
c            eq    c*16(8,361)    Intermediate matrix EQ
c            da    c*16(257,361)  Intermediate matrix DA
c
c     Output: n     c*16(361,361)  Intermediate matrix N
c             n2    c*16(361,361)  Intermediate matrix N2
c             eq    c*16(8,361)    Intermediate matrix EQ
c             da    c*16(257,361)  Intermediate matrix DA
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Return status.
c
      external fcl_normal
c
c     Input parameters.
c
      complex*16 gain(361)
c
c     Input/output parameters.
c
      complex*16 n(361,361),n2(361,361),eq(8,361),da(257,361)
c
c     Local variables.
c
      integer*4 pos,row,col
      complex*16 gain_conjg(361)

      fcl_calibrate = %loc(fcl_normal)
c
c     Calculate the conjugate of the gain function.
c
      do pos = 1,361
         gain_conjg(pos) = conjg(gain(pos))
      end do
c
c     Multiply N by the gain function on the left and by its conjugate on the 
c     right, and multiply N2 by the gain function on the left and on the right.
c
      do row = 1,361
         do col = 1,361
            n(row,col) = n(row,col) * gain(row) * gain_conjg(col)
            n2(row,col) = n2(row,col) * gain(row) * gain(col)
         end do
      end do
c
c     Multiply EQ by the conjugate of the gain function on the right.
c
      do row = 1,8
         do col = 1,361
            eq(row,col) = eq(row,col) * gain_conjg(col)
         end do  
      end do
c
c     Multiply DA by the gain function on the right.
c
      do row = 1,257
         do col = 1,361
            da(row,col) = da(row,col) * gain(col)
         end do   
      end do

      return
      end
