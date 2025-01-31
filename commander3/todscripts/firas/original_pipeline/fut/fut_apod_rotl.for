      integer*4 function fut_apod_rotl(aifg,apod,fft_len,peak,rifg)
c-------------------------------------------------------------------------------
c
c     Purpose: This function apodizes and rotates coadded interferograms 
c              using apodization functions read from the binary reference 
c              dataset FEX_APODL. 
c
c     Author: Alice Trenholme, GSC, 2/95, SER 12244
c             Created for new long spectrum pipeline routines 
c             (FIL, FFL, FSL, FCL).
c
c-------------------------------------------------------------------------------
c
c     Input: aifg(512)      r*4  Input real-valued interferogram.
c            apod(512)      r*8  Apodization function.
c            fft_len        i*4  Length of fast Fourier transform.
c            peak           i*4  Interferogram peak position.
c
c     Output: rifg(720)      r*8  Output apodized, rotated, double-precision
c                                 interferogram.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     External references.
c
      external fut_normal
c
c     Input parameters.
c
      real*4 aifg(512)
      real*8 apod(512)
      integer*4 fft_len,peak
c
c     Output parameter.
c
      real*8 rifg(720)
c
c     Local variables.
c
      integer*4 ifg_pos,sec_pos
      real*8 apifg(720)

      fut_apod_rotl = %loc(fut_normal)
c
c     Initialize the apodized interferogram and the apodized, rotated 
c     interferogram.
c
      call lib$movc5(0,,0,5760,apifg)
      call lib$movc5(0,,0,5760,rifg)
c
c     Apodize the interferogram, zero-padding from 513 to 720.
c
      do ifg_pos = 1, 512
         apifg(ifg_pos) = dble(aifg(ifg_pos))*apod(ifg_pos)
      end do
c
c     Rotate the interferogram so that its peak is at the first position.
c
      do ifg_pos = 1, fft_len + 1 - peak
         sec_pos = ifg_pos + peak - 1
         rifg(ifg_pos) = apifg(sec_pos)
      end do

      do ifg_pos = fft_len + 2 - peak, fft_len
         sec_pos = ifg_pos + peak - fft_len - 1
         rifg(ifg_pos) = apifg(sec_pos)
      end do

      return
      end
