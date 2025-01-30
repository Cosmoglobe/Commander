      integer*4 function fcc_calibrate(apod,etf,gain,peak,n,n2,eq,da)
c-------------------------------------------------------------------------------
c
c     Purpose: Apply calibration to four intermediate matrices.
c
c     Author: S. Alexander, HSTX, 7/93, SER 11189
c
c     Input: apod  r*8(512)       Apodization function
c            etf   c*16(512)      Electronics transfer function
c            gain  c*16(257)      Gain function
c            peak  i*4            Peak position
c            n     c*16(512,512)  Intermediate matrix N
c            n2    c*16(512,512)  Intermediate matrix N2
c            eq    c*16(8,512)    Intermediate matrix EQ
c            da    c*16(257,512)  Intermediate matrix DA
c
c     Output: n     c*16(512,512)  Intermediate matrix N
c             n2    c*16(512,512)  Intermediate matrix N2
c             eq    c*16(8,512)    Intermediate matrix EQ
c             da    c*16(257,512)  Intermediate matrix DA
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Return status.
c
      external fcc_normal
c
c     Input parameters.
c
      real*8 apod(512)
      complex*16 etf(512),gain(257)
      integer*4 peak
c
c     Input/output parameters.
c
      complex*16 n(512,512),n2(512,512),eq(8,512),da(257,512)
c
c     Local variables.
c
      integer*4 pos,row,col,dim
      complex*16 etf_conjg(512),gain_conjg(257),dspec(512)

      fcc_calibrate = %loc(fcc_normal)
c
c     Set variable DIM to 512 to avoid literals in calls to IMSL subroutines.
c
      dim = 512
c
c     Calculate the conjugate of the electronics transfer function and of the
c     gain function.
c
      do pos = 1,512
         etf_conjg(pos) = conjg(etf(pos))
      end do
      do pos = 1,257
         gain_conjg(pos) = conjg(gain(pos))
      end do
c
c     Multiply the intermediate matrix N by the electronics transfer function 
c     on the left and by its conjugate on the right. Multiply the intermediate
c     matrix N2 by the electronics transfer function on the left and on the
c     right. Multiply the intermediate matrix EQ by the conjugate of the 
c     electronics transfer function on the right, and multiply the intermediate
c     matrix DA by the electronics transfer function on the right.
c
      do row = 1,512
         do col = 1,512
            n(row,col) = n(row,col)*etf(row)*etf_conjg(col)
            n2(row,col) = n2(row,col)*etf(row)*etf(col)
         end do   
      end do
      do row = 1,8
         do col = 1,512
            eq(row,col) = eq(row,col)*etf_conjg(col)
         end do   
      end do
      do row = 1,257
         do col = 1,512
            da(row,col) = da(row,col)*etf(col)
         end do   
      end do
c
c     Replace each column of N by its inverse Fourier transform.
c
      do col = 1,512
         do row = 1,512
            dspec(row) = n(row,col)
         end do
         call dfftcb(dim,dspec,dspec)
         do row = 1,512
            n(row,col) = dspec(row)
         end do
      end do
c
c     Replace N2 by its inverse two-sided Fourier transform.
c
      call dfft2b(dim,dim,n2,dim,n2,dim)
c
c     Replace each row of N and each row of EQ by its forward Fourier
c     transform.
c
      do row = 1,512
         do col = 1,512
            dspec(col) = n(row,col)
         end do
         call dfftcf(dim,dspec,dspec)
         do col = 1,512
            n(row,col) = dspec(col)
         end do
      end do
      do row = 1,8
         do col = 1,512
            dspec(col) = eq(row,col)
         end do
         call dfftcf(dim,dspec,dspec)
         do col = 1,512
            eq(row,col) = dspec(col)
         end do
      end do
c
c     Replace each row of DA by its inverse Fourier transform.
c
      do row = 1,257
         do col = 1,512
            dspec(col) = da(row,col)
         end do
         call dfftcb(dim,dspec,dspec)
         do col = 1,512
            da(row,col) = dspec(col)
         end do
      end do
c
c     Multiply N and N2 by the apodization function on the left and right;
c     multiply EQ and DA by the apodization function on the right.
c
      do row = 1,512
         do col = 1,512
            n(row,col) = n(row,col)*apod(row)*apod(col)
            n2(row,col) = n2(row,col)*apod(row)*apod(col)
         end do   
      end do
      do row = 1,8
         do col = 1,512
            eq(row,col) = eq(row,col)*apod(col)
         end do   
      end do
      do row = 1,257
         do col = 1,512
            da(row,col) = da(row,col)*apod(col)
         end do   
      end do
c
c     Replace each row of the four intermediate matrices by the row rotated 
c     so that the interferogram peak position is in the first position.
c
      do row = 1,512
         do col = 1,512
            dspec(col) = n(row,col)
         end do
         do col = 1,512
            n(row,col) = dspec(1+mod(col+peak-2,512))
         end do
         do col = 1,512
            dspec(col) = n2(row,col)
         end do
         do col = 1,512
            n2(row,col) = dspec(1+mod(col+peak-2,512))
         end do
      end do
      do row = 1,8
         do col = 1,512
            dspec(col) = eq(row,col)
         end do
         do col = 1,512
            eq(row,col) = dspec(1+mod(col+peak-2,512))
         end do
      end do
      do row = 1,257
         do col = 1,512
            dspec(col) = da(row,col)
         end do
         do col = 1,512
            da(row,col) = dspec(1+mod(col+peak-2,512))
         end do
      end do
c
c     Replace each column of N and N2 by the column rotated so that the 
c     interferogram peak position is the first position.
c
      do col = 1,512
         do row = 1,512
            dspec(row) = n(row,col)
         end do
         do row = 1,512
            n(row,col) = dspec(1+mod(row+peak-2,512))
         end do
         do row = 1,512
            dspec(row) = n2(row,col)
         end do
         do row = 1,512
            n2(row,col) = dspec(1+mod(row+peak-2,512))
         end do
      end do
c
c     Replace each column of N by its forward Fourier transform.
c
      do col = 1,512
         do row = 1,512
            dspec(row) = n(row,col)
         end do
         call dfftcf(dim,dspec,dspec)
         do row = 1,512
            n(row,col) = dspec(row)
         end do
      end do
c
c     Replace N2 by its forward two-sided Fourier transform.
c
      call dfft2d(dim,dim,n2,dim,n2,dim)
c
c     Replace each row of N and EQ by its inverse Fourier transform.
c
      do row = 1,512
         do col = 1,512
            dspec(col) = n(row,col)
         end do
         call dfftcb(dim,dspec,dspec)
         do col = 1,512
            n(row,col) = dspec(col)
         end do
      end do
      do row = 1,8
         do col = 1,512
            dspec(col) = eq(row,col)
         end do
         call dfftcb(dim,dspec,dspec)
         do col = 1,512
            eq(row,col) = dspec(col)
         end do
      end do
c
c     Replace each row of DA by its forward Fourier transform.
c
      do row = 1,257
         do col = 1,512
            dspec(col) = da(row,col)
         end do
         call dfftcf(dim,dspec,dspec)
         do col = 1,512
            da(row,col) = dspec(col)
         end do
      end do
c
c     Divide N by the electronics transfer function on the left and by its 
c     conjugate on the right, multiply N by the gain function on the left and 
c     by its conjugate on the right, and normalize N by dividing by 512
c     squared. Divide N2 by the electronics transfer function on the left and on
c     the right, multiply N2 by the gain function on the left and on the right,
c     and normalize N2 by dividing by 512 squared.
c
      do row = 1,257
         if (etf(row) .ne. 0.0) then
            do col = 1,257
               if (etf(col) .ne. 0.0) then
                  n(row,col) = ((n(row,col) / (etf(row)*etf_conjg(col))) *
     .                          (gain(row) * gain_conjg(col))) / 262144.0
                  n2(row,col) = ((n2(row,col) / (etf(row)*etf(col))) *
     .                           (gain(row) * gain(col))) / 262144.0
               else
                  n(row,col) = ((n(row,col) / etf(row)) *
     .                          (gain(row) * gain_conjg(col))) / 262144.0
                  n2(row,col) = ((n2(row,col) / etf(row)) *
     .                           (gain(row) * gain(col))) / 262144.0
               end if
            end do
         else
            do col = 1,257
               if (etf(col) .ne. 0.0) then
                  n(row,col) = ((n(row,col) / etf_conjg(col)) *
     .                          (gain(row) * gain_conjg(col))) / 262144.0
                  n2(row,col) = ((n2(row,col) / etf(col)) *
     .                           (gain(row) * gain(col))) / 262144.0
               else
                  n(row,col) = (n(row,col) * (gain(row) * gain_conjg(col))) / 
     .                          262144.0
                  n2(row,col) = (n2(row,col) * (gain(row) * gain(col))) / 
     .                           262144.0
               end if
            end do
         end if
c
c     Divide DA by the electronics transfer function on the right, multiply DA 
c     by the gain function on the right, and normalize DA by dividing by 512.
c
         do col = 1,257
            if (etf(col) .ne. 0.0) then
               da(row,col) = ((da(row,col) / etf(col)) *
     .                         gain(col)) / 512.0
            else
               da(row,col) = (da(row,col) * gain(col)) / 512.0
            end if
         end do   
      end do
c
c     Divide EQ by the conjugate of the electronics transfer function on the 
c     right, multiply EQ by the conjugate of the gain function on the right,
c     and normalize EQ by dividing by 512.
c
      do row = 1,8
         do col = 1,257
            if (etf(col) .ne. 0.0) then
               eq(row,col) = ((eq(row,col) / etf_conjg(col)) *
     .                         gain_conjg(col)) / 512.0
            else
               eq(row,col) = (eq(row,col) * gain_conjg(col)) / 512.0
            end if
         end do  
      end do

      return
      end
