      integer*4 function fut_apod_recnuml(channel,scan_mode,fake_it,sci_mode,
     .                                    adds_per_group,linearized,recnum)
c-------------------------------------------------------------------------------
c
c     Purpose: Return the record number for the appropriate apodization 
c              function from the direct-access apodization file for revised
c              pipeline facilities (FIL, FSL, FFL, FCL).
c
c     Authors: Rich Isaacman, ARC, 6/89
c              Alice Trenholme, GSC, 2/95, SER 12244
c
c     Input: channel            i*4  Value of channel, 1-4 = RH-LL.
c            scan_mode          i*4  Value of scan mode, 1-6 = SS-FL.
c            fake_it            i*4  Value of fake-it bit, 0-1.
c            sci_mode           i*4  Value of science mode, 0-4.
c            adds_per_group     i*4  Value of adds per group, 1-12.
c            linearized         i*4  Whether apodization function should be for
c                                    linearized interferograms, 0 = no, 1 = yes.
c
c     Output: recnum            i*4  Record number of appropriate
c                                    apodization function.
c
c     Modifications:
c
c-------------------------------------------------------------------------------
      implicit none
c
c     Return statuses.
c
      external fut_normal
c
c     Input parameters.
c
      integer*4 channel,scan_mode,fake_it
      integer*4 sci_mode,adds_per_group,linearized
c
c     Output parameters.
c
      integer*4 recnum
c
c     Local variables.
c
      integer*4 ch,sc

      fut_apod_recnuml = %loc(fut_normal)
c
c     Determine whether channel is high or low.
c
      if ((channel .eq. 1) .or. (channel .eq. 3)) then 
          ch = 1 
      else 
          ch = 0
      end if
c
c     Determine whether digital filters are on or off.
c
      if ((sci_mode .eq. 1) .or. (sci_mode .eq. 3)) then 
          sc = 0
      else 
          sc = 1
      end if
c
c     Apodization function for fake-it data is the last one (433); functions
c     for FS and FL data are after 384.
c
      if (fake_it .eq. 1) then
         recnum = 433
      else if (scan_mode .ge. 5) then
         recnum = 384 + 4*(adds_per_group-1) + 2*sc + linearized + 1
      else
         recnum = 32*(adds_per_group-1) + 16*(mod(scan_mode-1,2)) + 
     .            8*(int((scan_mode-1)/2)) + 4*ch + 2*sc + linearized + 1
      end if

      return
      end
