      Integer*4 Function FUT_Ave_Bin_Times(bin_time_array,num,
     .                                     weights,bin_time_ave)

c------------------------------------------------------------------------------
c     Purpose: To add and average quadword binary times (the adt).
c
c     Input:
c
c       1. Integer*4 bin_time_array(2) of adts
c       2. Integer*4 num of elements in the array
c       3. Integer*4 weights(num) of summed value
c
c     Output:
c
c       1. Integer*4 bin_time_ave(2) binary average of the adts
c
c     Requirement:
c       This algorithm assumes that the sum of the differences
c       of the binary values is < (2**31) - 1 [ie is representable
c       as a positive i*4 number].
c
c       Not a good assumption, as it turned out; this limits the sum of the
c       differences to ~3.5 min.  It is now represented in an i*4 array of
c       dimension 2.
c             D. Bouler, STX, Mar 09, 1990.
c
c
c     Author: J. Durachta
c             ST Systems Corporation
c             February 10,1987
c
c     Algorithm:
c
c       Store the first adt
c       Loop over the rest of the adts
c         Call lib$subx to find the difference
c         between the first adt and the others
c         Multiply the differences by their weight (1 or 0)
c         Sum the differences
c       End loop
c
c       Divide the differnce by the num to produce
c       the average difference
c
c       Add the average difference back to the first
c       adt to produce the average binary time
c
c       end
c------------------------------------------------------------------------------
c
c       Changes:
c
c       3/9/87    J. Durachta - Logic altered to prevent inclusion of channels
c                 with 0 IFGs.
c
c       3/26/87   J. Durachta - Logic altered to exclude first IFG from loop
c                 (ie. do i = j+1,num).
c
c       7/30/87   J. Durachta: Variable ave_bin_diff(2) initialized to 0 to
c                              prevent sign contamination from previous use.
c
c       2/19/88   R. Kummerer  Return success status FUT_NORMAL.
c
c       10/28/88  R. Kummerer  SPR 2702, Convert average time to CT "precision".
c
c       03/09/90  D. Bouler    spr 6440 Use i*4 array of dimension 2 to
c                              hold sum of differences. Use VAX extended
c                              arithmetic calls instead of bit shifting.
c
c	12/10/90  F. Shuman & L. Rosen   SPR 7821.  To multiply each time
c		difference by its weight, LIB$EMul was being used.  Its known
c		bug was being worked around, but in a way that failed for time
c		differences greater than 2^31 * 10^-7 sec, or about 3.5 min. 
c		LIB$MultF_Delta_Time could not handle delta times of both signs,
c		besides which, it would give only R*4 (roughly 3-byte) accuracy.
c		LIB$EMul, while buggy, proved to be patchable.  Secondly, to do
c		the final division of the sum of time differences by the sum of
c		the weights, the existing replacement for LIB$EDiv (using
c		Real*8, AUT_DFloat2ADT, and AUT_ADT2DFloat) didnt work on
c		delta times, which always have their sign bit set.  LIB$EDiv,
c		though even buggier than LIB$EMul, was patched into service.
c------------------------------------------------------------------------------

        Implicit None

        Integer*2 num                  !number of adts to be averaged
        Integer*4 bin_time_array(2,num)!array of binary times
        Integer*4 weights(num)         !statistical weight of summed value
        Integer*4 bin_time_ave(2)      !average of adts

        Logical*1 loop, negflag

        Integer*4 i,j                  !a counter
        Integer*4 div                  !divisor to produce average
        Integer*4 temp_diff(2)         !temp to hold most sig part of bin_diff
        Integer*4 bin_diff(2)          !sum of the diffs of the adts
        Integer*4 bin_num(2), temp(2)  !temporary holding vars
        Integer*4 diff_sum(2)          !sum of bin time diffs
        Integer*4 ave_bin_diff(2)      !bin ave of time diffs

        Character*14 gmt               !time conversion parameter

        Integer*4 wtd_diff(2)          !adt difference times the weight
        Integer*4 status               !return status

        Integer*4 quo, rem, top(2), offset(2)  !dummy variables, and...
        Integer*4 zero /0/, one(2) /1,0/       !constants... for Lib$EMul,
        Integer*4 half(2) /'80000000'X,0/      ! Lib$AddX, and Lib$SubX calls

        Integer*4 Lib$EDiv, Lib$EMul, Lib$AddX, Lib$SubX

        External FUT_Normal, FUT_OflWrnA, FUT_OflWrnB
        External SS$_Normal, SS$_IntOvf, SS$_IntDiv


        If (num .Gt. 100) Call Lib$signal(fut_oflwrna,%val(0),fut_oflwrnb)

!    Initialize diff sum and divisor

        ave_bin_diff(2) = 0
        diff_sum(1) = 0
        diff_sum(2) = 0
        j = 0
        loop = .true.

!       Find first non zero weight. This has no effect in average over
!       single channel, but excludes 0 IFG channel in grand average if
!       it is first.

        Do While(loop)
          j = j + 1
          If (weights(j) .Ne. 0) loop = .False.
        End Do

        div = weights(j)

!    Initialize bin ave time.

        bin_time_ave(1) = bin_time_array(1,j)
        bin_time_ave(2) = bin_time_array(2,j)

        If (num .Gt. 1) Then
          Do i = j+1,num

!           Exclude channels with 0 IFGs.

            If (weights(i) .Ne. 0) Then

              wtd_diff(1) = 0
              wtd_diff(2) = 0

              status = lib$subx(bin_time_array(1,i), bin_time_ave, bin_diff)

c  Get weighted difference = weight(i) * diff.  Calculation is done in two parts
c  with LIB$EMul.  Now, for some unfathomable reason, when the Extended-Multiply
c  routine, LIB$EMul, is given a multiplicand (the double-longword argument)
c  with the high-order bit of the low-order longword set, it renders a product
c  correct for a high-order longword one less than the one it is supplied.  We
c  must check the sign of the low-order component of bin_diff to correct this.
c  Then multiply the lower (1) order component by the weight; store in wtd_diff.
c  Multiply the higher (2) order component by the weight, adding the carry term
c  (wtd_diff(2)), then store in temp_diff, which becomes the higher order term
c  of wtd_diff.  Use LIB$AddX to accumulate the sum of weighted differences.

              If (bin_diff(1) .Lt. 0) bin_diff(2) = bin_diff(2) + 1
              status = lib$emul(bin_diff(1), weights(i), zero, wtd_diff)
              status = lib$emul(bin_diff(2), weights(i), wtd_diff(2), temp_diff)
              wtd_diff(2) = temp_diff(1)

              status = lib$addx(diff_sum, wtd_diff, diff_sum)

              div = div + weights(i)

            End If
          End Do
c
c   Lib$EDiv also needs a lot of help.  We want to use it here with dividends
c (numerators) of either sign, but negative numerators result in Integer Over-
c flow, with no division being done.  So we do a double-longword negate in this
c case, before and after the divide.  This consists of doing a 1's-complement of
c each longword (L=-1-L), then Extended Addition of a quadword "1".
c   Secondly, Lib$EDiv seems to suffer a deformity similar to that of Lib$EMul;
c it is supposedly designed to divide a Quadword dividend by a Longword divisor
c and return Longword quotient and remainder.  This seems to imply that the
c dividend can range up to (divisor*(2^32) - 1), the upper limit for a 32-bit
c quotient.  But the actual limit is that of a 31-bit quotient.  This causes a
c very annoying interruption in each Longword of a long-division process.
c Praise Allah we are only using it twice!

c Object:  to divide DIFF_SUM(2:1) by DIV and put result into AVE_BIN_DIFF(2:1)
c   If DIFF_SUM is negative, negate it and remember this, so we can negate the
c   result.  TOP(2:1), which gets the numerator, is picked apart in order to
c   submit it to EDiv twice, once for each Longword (LW) of the answer.
c
          If (diff_sum(2) .Lt. 0) Then
            negflag = .True.
            top(1) = -1 - diff_sum(1)
            top(2) = -1 - diff_sum(2)
            status = Lib$AddX(top, one, top)
          Else
            negflag = .False.
            top(1) = diff_sum(1)
            top(2) = diff_sum(2)
          End If
c
c First, set up division for the high-order LW.  For this to overflow, the sum
c   of diffs would have to be 2^63 * 10^-7 sec =~ 10^12 sec =~30000 yr
c
          temp(2) = 0
          temp(1) = top(2)
          status = LIB$EDiv ( div, temp, quo, rem )
          ave_bin_diff(2) = quo
c
c Next, set up division for the low-order LW.  If this overflows, offset the
c   dividend (TEMP) by subtracting the double-LW, [2^31(=HALF) * div], redo the
c   division, and add HALF back to the quotient.  Note that if div is even,
c   HALF*div will be [Int(div/2), 0]; if odd, this will be [Int(div/2), HALF].
c   (Given in the order [LW(2), LW(1)], since LW(2) is the high-order part.)
c
          temp(2) = rem
          temp(1) = top(1)
          status = LIB$EDiv ( div, temp, quo, rem )
          If (status .Eq. %Loc(SS$_IntOvf)) Then
             offset(2) = div/2
             If (div .Eq. 2*(div/2)) Then
                offset(1) = 0
             Else
                offset(1) = half(1)
             End If
             status = LIB$SubX ( temp, offset, temp)
             status = LIB$EDiv ( div, temp, quo, rem )
             ave_bin_diff(1) = quo
             status = LIB$AddX ( ave_bin_diff, half, ave_bin_diff)
          Else
             ave_bin_diff(1) = quo
          End If
c
c Now negate the result if we negated the incoming dividend.
c
          If (negflag) Then
             ave_bin_diff(1) = -1 - ave_bin_diff(1)
             ave_bin_diff(2) = -1 - ave_bin_diff(2)
             status = LIB$AddX (ave_bin_diff, one, ave_bin_diff)
          End If

c Add the average binary time diff back to initial time to get average of times.

          status = lib$addx(bin_time_ave, ave_bin_diff, bin_time_ave)

        endif

!     Convert the average time to ASCII format and back to binary format.
!     This ensures the low order bits of the average time will be the "same"
!     as if the average time came directly from a GMT ASCII string.  Apparently
!     low order bits in a binary time can differ and yet produce the ASCII
!     string with CT_BINARY_TO_GMT.  By ensuring the same low order bits,
!     CT utilities such as TIME_GE and TIME_LE will give consistent results.

        call ct_binary_to_gmt(bin_time_ave,gmt)
        call ct_gmt_to_binary(gmt,bin_time_ave)

        fut_ave_bin_times = %loc(fut_normal)

        return
        end
