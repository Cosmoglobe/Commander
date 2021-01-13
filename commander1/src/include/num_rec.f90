!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module num_rec
use healpix_types
use misc_utils
implicit none
private

public :: dsvbksb, dsvdcmp, isort, sort, iindexx, indexx, gammp, gammln

contains

! The routines PYTHAG, SVD and IPSORT were taken from SLATEC:
! http://www.netlib.org/slatec/src/<routine>.f
! The SLATEC routines are in the public domain.
!
! Changes by MR:
!  - processed with f2f90 sed script
!  - replaced REAL with REAL(DP)
!  - replaced xx.0E0 with xx.0_DP
!  - rewrote some loops to avoid shared CONTINUE statements and loop end
!    statements which are not CONTINUE
!  - removed calls to XERMSG; the error code is checked in the calling routine
!
! Changes by EH:
!  psort has been edited from slatec's IPSORT
!   IX --> XX
!   ITEMP --> XTEMP
!
!DECK PYTHAG
      FUNCTION PYTHAG (A, B)
!***BEGIN PROLOGUE  PYTHAG
!***SUBSIDIARY
!***PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PYTHAG
      REAL(DP) PYTHAG
      REAL(DP) A,B
!
      REAL(DP) P,Q,R,S,T
!***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END FUNCTION

!DECK SVD
      SUBROUTINE SVD (NM, M, N, A, W, MATU, U, MATV, V, IERR, RV1)
!***BEGIN PROLOGUE  SVD
!***SUBSIDIARY
!***PURPOSE  Perform the singular value decomposition of a rectangular
!            matrix.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVD-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure SVD,
!     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
!
!     This subroutine determines the singular value decomposition
!          T
!     A=USV  of a REAL M by N rectangular matrix.  Householder
!     bidiagonalization and a variant of the QR algorithm are used.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, U and V, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!          Note that NM must be at least as large as the maximum
!          of M and N.
!
!        M is the number of rows of A and U.
!
!        N is the number of columns of A and U and the order of V.
!
!        A contains the rectangular input matrix to be decomposed.  A is
!          a two-dimensional REAL array, dimensioned A(NM,N).
!
!        MATU should be set to .TRUE. if the U matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATU is a LOGICAL variable.
!
!        MATV should be set to .TRUE. if the V matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATV is a LOGICAL variable.
!
!     On Output
!
!        A is unaltered (unless overwritten by U or V).
!
!        W contains the N (non-negative) singular values of A (the
!          diagonal elements of S).  They are unordered.  If an
!          error exit is made, the singular values should be correct
!          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
!          REAL array, dimensioned W(N).
!
!        U contains the matrix U (orthogonal column vectors) of the
!          decomposition if MATU has been set to .TRUE.  Otherwise,
!          U is used as a temporary array.  U may coincide with A.
!          If an error exit is made, the columns of U corresponding
!          to indices of correct singular values should be correct.
!          U is a two-dimensional REAL array, dimensioned U(NM,N).
!
!        V contains the matrix V (orthogonal) of the decomposition if
!          MATV has been set to .TRUE.  Otherwise, V is not referenced.
!          V may also coincide with A if U does not.  If an error
!          exit is made, the columns of V corresponding to indices of
!          correct singular values should be correct.  V is a two-
!          dimensional REAL array, dimensioned V(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          K          if the K-th singular value has not been
!                     determined after 30 iterations.
!
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
!
!     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  SVD
!
      INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
      REAL(DP) A(NM,*),W(*),U(NM,*),V(NM,*),RV1(*)
      REAL(DP) C,F,G,H,S,X,Y,Z,SCALE,S1
!      REAL(DP) PYTHAG
      LOGICAL MATU,MATV
!
!***FIRST EXECUTABLE STATEMENT  SVD
      IERR = 0
!
      DO I = 1, M
!
         DO J = 1, N
            U(I,J) = A(I,J)
         END DO
      END DO
!     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0_DP
      SCALE = 0.0_DP
      S1 = 0.0_DP
!
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0_DP
         S = 0.0_DP
         SCALE = 0.0_DP
         IF (I .GT. M) GO TO 210
!
         DO K = I, M
           SCALE = SCALE + ABS(U(K,I))
         END DO
!
         IF (SCALE .EQ. 0.0_DP) GO TO 210
!
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
!
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
!
         DO 150 J = L, N
            S = 0.0_DP
!
            DO K = I, M
              S = S + U(K,I) * U(K,J)
            END DO
!
            F = S / H
!
            DO K = I, M
              U(K,J) = U(K,J) + F * U(K,I)
            END DO
  150    CONTINUE
!
  190    DO K = I, M
           U(K,I) = SCALE * U(K,I)
         END DO
!
  210    W(I) = SCALE * G
         G = 0.0_DP
         S = 0.0_DP
         SCALE = 0.0_DP
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
!
         DO K = L, N
           SCALE = SCALE + ABS(U(I,K))
         END DO
!
         IF (SCALE .EQ. 0.0_DP) GO TO 290
!
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
!
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
!
         DO K = L, N
           RV1(K) = U(I,K) / H
         END DO
!
         IF (I .EQ. M) GO TO 270
!
         DO 260 J = L, M
            S = 0.0_DP
!
            DO K = L, N
              S = S + U(J,K) * U(I,K)
            END DO
!
            DO K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
            END DO
  260    CONTINUE
!
  270    DO K = L, N
           U(I,K) = SCALE * U(I,K)
         END DO
!
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
!     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0_DP) GO TO 360
!
         DO J = L, N
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
           V(J,I) = (U(I,J) / U(I,L)) / G
         END DO
!
         DO 350 J = L, N
            S = 0.0_DP
!
            DO K = L, N
              S = S + U(I,K) * V(K,J)
            END DO
!
            DO K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
            END DO
  350    CONTINUE
!
  360    DO 380 J = L, N
            V(I,J) = 0.0_DP
            V(J,I) = 0.0_DP
  380    CONTINUE
!
  390    V(I,I) = 1.0_DP
         G = RV1(I)
         L = I
  400 CONTINUE
!     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
!     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
!
      DO 500 II = 1, MN
         I = MN + 1 - II
         L = I + 1
         G = W(I)
         IF (I .EQ. N) GO TO 430
!
         DO J = L, N
           U(I,J) = 0.0_DP
         END DO
!
  430    IF (G .EQ. 0.0_DP) GO TO 475
         IF (I .EQ. MN) GO TO 460
!
         DO 450 J = L, N
            S = 0.0_DP
!
            DO K = L, M
              S = S + U(K,I) * U(K,J)
            END DO
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
!
            DO K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
            END DO
  450    CONTINUE
!
  460    DO J = I, M
           U(J,I) = U(J,I) / G
         END DO
!
         GO TO 490
!
  475    DO J = I, M
           U(J,I) = 0.0_DP
         END DO
!
  490    U(I,I) = U(I,I) + 1.0_DP
  500 CONTINUE
!     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 CONTINUE
!     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
!     .......... TEST FOR SPLITTING.
!                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            IF (S1 + ABS(RV1(L)) .EQ. S1) GO TO 565
!     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
            IF (S1 + ABS(W(L1)) .EQ. S1) GO TO 540
  530    CONTINUE
!     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0_DP
         S = 1.0_DP
!
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (S1 + ABS(F) .EQ. S1) GO TO 565
            G = W(I)
            H = PYTHAG(F,G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (.NOT. MATU) GO TO 560
!
            DO 550 J = 1, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  550       CONTINUE
!
  560    CONTINUE
!     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
!     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5_DP * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = PYTHAG(F,1.0_DP)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
!     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0_DP
         S = 1.0_DP
!
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = PYTHAG(F,H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
            IF (.NOT. MATV) GO TO 575
!
            DO 570 J = 1, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
!
  575       Z = PYTHAG(F,H)
            W(I1) = Z
!     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0_DP) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
!
            DO 590 J = 1, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  590       CONTINUE
!
  600    CONTINUE
!
         RV1(L) = 0.0_DP
         RV1(K) = F
         W(K) = X
         GO TO 520
!     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0_DP) GO TO 700
!     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
         IF (.NOT. MATV) GO TO 700
!
         DO J = 1, N
           V(J,K) = -V(J,K)
         END DO
!
  700 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO A
!                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END SUBROUTINE


!DECK IPSORT
      SUBROUTINE IPSORT (IX, N, IPERM, KFLAG, IER)
!***BEGIN PROLOGUE  IPSORT
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1A, N6A2A
!***TYPE      INTEGER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   IPSORT returns the permutation vector IPERM generated by sorting
!   the array IX and, optionally, rearranges the values in IX.  IX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that IX(IPERM(I)) is the Ith value in the
!   rearrangement of IX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between IPSORT and its active sorting equivalent
!   ISORT is that the data are referenced indirectly rather than
!   directly.  Therefore, IPSORT should require approximately twice as
!   long to execute as ISORT.  However, IPSORT is more general.
!
!   Description of Parameters
!      IX - input/output -- integer array of values to be sorted.
!           If ABS(KFLAG) = 2, then the values in IX will be
!           rearranged on output; otherwise, they are unchanged.
!      N  - input -- number of values in array IX to be sorted.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              IX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting IX in increasing order and sort IX also.
!            =  1  means return the permutation vector resulting from
!                  sorting IX in increasing order and do not sort IX.
!            = -1  means return the permutation vector resulting from
!                  sorting IX in decreasing order and do not sort IX.
!            = -2  means return the permutation vector resulting from
!                  sorting IX in decreasing order and sort IX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   810801  Further modified by David K. Kahaner.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  IPSORT
!     .. Scalar Arguments ..
      INTEGER IER, KFLAG, N
!     .. Array Arguments ..
      INTEGER IPERM(*), IX(*)
!     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, INDX, INDX0, ISTRT, ITEMP, J, K, KK, L, LM, LMT, M, &
              NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!      EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  IPSORT
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         IER = 1
!         CALL XERMSG ('SLATEC', 'IPSORT', &
!          'The number of values to be sorted, N, is not positive.', &
!          IER, 1)
         RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         IER = 2
!         CALL XERMSG ('SLATEC', 'IPSORT', &
!          'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!          IER, 1)
         RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. 1) RETURN
!
!     Alter array IX to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
         DO 20 I=1,NN
            IX(I) = -IX(I)
   20    CONTINUE
      ENDIF
!
!     Sort IX only
!
      M = 1
      I = 1
      J = NN
      R = .375_DP
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437_DP) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875_DP
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (IX(IPERM(I)) .GT. IX(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (IX(IPERM(J)) .LT. IX(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
         IF (IX(IPERM(I)) .GT. IX(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
      IF (IX(IPERM(L)) .GT. IX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
      IF (IX(IPERM(K)) .LT. IX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
!
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (IX(IPERM(I)) .LE. IX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+1) = IPERM(K)
      K = K-1
!
      IF (IX(LM) .LT. IX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -1) THEN
         DO 130 I=1,NN
            IX(I) = -IX(I)
  130    CONTINUE
      ENDIF
!
!     Rearrange the values of IX if desired
!
      IF (KK .EQ. 2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
         DO 150 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0)  THEN
               INDX = ISTRT
               INDX0 = INDX
               ITEMP = IX(ISTRT)
  140          IF (IPERM(INDX) .GT. 0)  THEN
                  IX(INDX) = IX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               IX(INDX0) = ITEMP
            ENDIF
  150    CONTINUE
!
!        Revert the signs of the IPERM values
!
         DO 160 I=1,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE
!
      ENDIF
!
      RETURN
      END SUBROUTINE
!=====================================
      SUBROUTINE XPSORT (XX, N, IPERM, KFLAG, IER)
! ========================================================
! xpsort: edited from Slatec's IPSORT by Eric Hivon
! deal with double precision real instead of integer
! ========================================================
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   XPSORT returns the permutation vector IPERM generated by sorting
!   the array XX and, optionally, rearranges the values in XX.  XX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that XX(IPERM(I)) is the Ith value in the
!   rearrangement of XX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between SORT and its active sorting equivalent
!   ISORT is that the data are referenced indirectly rather than
!   directly.  Therefore, SORT should require approximately twice as
!   long to execute as ISORT.  However, SORT is more general.
!
!   Description of Parameters
!      XX - input/output -- integer array of values to be sorted.
!           If ABS(KFLAG) = 2, then the values in XX will be
!           rearranged on output; otherwise, they are unchanged.
!      N  - input -- number of values in array XX to be sorted.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              XX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting XX in increasing order and sort XX also.
!            =  1  means return the permutation vector resulting from
!                  sorting XX in increasing order and do not sort XX.
!            = -1  means return the permutation vector resulting from
!                  sorting XX in decreasing order and do not sort XX.
!            = -2  means return the permutation vector resulting from
!                  sorting XX in decreasing order and sort XX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   810801  Further modified by David K. Kahaner.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!     .. Scalar Arguments ..
      INTEGER IER, KFLAG, N
!     .. Array Arguments ..
      INTEGER IPERM(*)
      real(DP) :: XX(*)
!     .. Local Scalars ..
      real(DP) :: xtemp
      REAL R
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, &
              NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!      EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SORT
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         IER = 1
!         CALL XERMSG ('SLATEC', 'SORT', &
!          'The number of values to be sorted, N, is not positive.', &
!          IER, 1)
         RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         IER = 2
!         CALL XERMSG ('SLATEC', 'SORT', &
!          'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!          IER, 1)
         RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. 1) RETURN
!
!     Alter array XX to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
         DO 20 I=1,NN
            XX(I) = -XX(I)
   20    CONTINUE
      ENDIF
!
!     Sort XX only
!
      M = 1
      I = 1
      J = NN
      R = .375_DP
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437_DP) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875_DP
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (XX(IPERM(I)) .GT. XX(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (XX(IPERM(J)) .LT. XX(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
         IF (XX(IPERM(I)) .GT. XX(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
      IF (XX(IPERM(L)) .GT. XX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
      IF (XX(IPERM(K)) .LT. XX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
!
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (XX(IPERM(I)) .LE. XX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+1) = IPERM(K)
      K = K-1
!
      IF (XX(LM) .LT. XX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -1) THEN
         DO 130 I=1,NN
            XX(I) = -XX(I)
  130    CONTINUE
      ENDIF
!
!     Rearrange the values of XX if desired
!
      IF (KK .EQ. 2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
         DO 150 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0)  THEN
               INDX = ISTRT
               INDX0 = INDX
               XTEMP = XX(ISTRT)
  140          IF (IPERM(INDX) .GT. 0)  THEN
                  XX(INDX) = XX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               XX(INDX0) = XTEMP
            ENDIF
  150    CONTINUE
!
!        Revert the signs of the IPERM values
!
         DO 160 I=1,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE
!
      ENDIF
!
      RETURN
    END SUBROUTINE
!=========================================
  subroutine dsvbksb(u,w,v,m,n,mp,np,b,x)
    integer(i4b), intent(in) :: m,mp,n,np
    real(dp), intent(in) ::  b(mp),u(mp,np),v(np,np),w(np)
    real(dp), intent(out) :: x(np)
    integer(i4b) ::  i,j,jj
    real(dp) ::  s,tmp(n)
    do j=1,n
       s = 0d0
       if (w(j)/=0d0) then
          do i=1,m
             s=s+u(i,j)*b(i)
          end do
          s=s/w(j)
       endif
       tmp(j)=s
    end do
    do j=1,n
       s=0.0d0
       do jj=1,n
          s=s+v(j,jj)*tmp(jj)
       end do
       x(j)=s
    end do
  end subroutine dsvbksb

  subroutine dsvdcmp(a,m,n,mp,np,w,v)
    integer(i4b), intent(in):: m,mp,n,np
    real(dp), intent(inout) ::  a(mp,np),v(np,np),w(np)
    integer ierr
    real(dp) rv1(np),u(mp,np)

    call svd(mp,m,n,a,w,.true.,u,.true.,v,ierr,rv1)
    call assert(ierr==0,"error in svd()")
    a=u
  end subroutine dsvdcmp

  subroutine isort(n,ra)
    integer, intent(in) :: n
    integer, intent(inout), dimension(1:) :: ra
    integer indx(n),ier
    call ipsort (ra,n,indx,2,ier)
    call assert (ier==0, "error in ipsort()")
  end subroutine isort

  subroutine sort(n,ra)
    integer, intent(in) :: n
    real(dp), intent(inout), dimension(1:) :: ra
    integer indx(n),ier
    call xpsort (ra,n,indx,2,ier)
    call assert (ier==0, "error in xpsort()")
  end subroutine sort

  SUBROUTINE iindexx(n,arr,indx)
    INTEGER,                 intent(IN)  :: n
    INTEGER, dimension(1:n), intent(IN)  :: arr
    integer, dimension(1:n), intent(out) :: indx
    integer ier

    call ipsort (arr,n,indx,1,ier)
    call assert (ier==0, "error in ipsort()")
  END SUBROUTINE iindexx

  SUBROUTINE indexx(n,arr,indx)
    INTEGER,                 intent(IN)  :: n
    real(DP), dimension(1:n), intent(IN)  :: arr
    integer,  dimension(1:n), intent(out) :: indx
    integer ier

    call xpsort (arr,n,indx,1,ier)
    call assert (ier==0, "error in xpsort()")
  END SUBROUTINE indexx

  function gammp(a, x, gln)
    implicit none

    real(dp), intent(in)            :: a, x
    real(dp)                        :: gammp
    real(dp), intent(out), optional :: gln

    if (x < a+1.d0) then
       if (present(gln)) then
          gammp = gser(a,x,gln)
       else
          gammp = gser(a,x)
       end if
    else
       if (present(gln)) then
          gammp = 1.d0 - gcf(a,x,gln)
       else
          gammp = 1.d0 - gcf(a,x)
       end if
    end if

  end function gammp
  
  function gser(a, x, gln)
    implicit none

    real(dp), intent(in)            :: a, x
    real(dp), intent(out), optional :: gln
    real(dp)                        :: gser

    integer(i4b), parameter :: ITMAX = 100
    real(dp), parameter :: eps = epsilon(x)
    integer(i4b) :: n
    real(dp)     :: ap, del, summ
    
    if (x == 0.d0) then
       gser = 0.d0
       return
    end if
    ap = a
    summ = 1.d0 / a
    del  = summ
    do n = 1, ITMAX
       ap   = ap + 1.d0
       del  = del * x/ap
       summ = summ + del
       if (abs(del) < abs(summ)*eps) exit
    end do
    if (n > ITMAX) then
       write(*,*) 'a too large, ITMAX too small in gser.'
       stop
    end if
    if (present(gln)) then
       gln = gammln(a)
       gser = summ*exp(-x+a*log(x)-gln)
    else
       gser = summ*exp(-x+a*log(x)-gammln(a))
    end if

  end function gser

  function gcf(a, x, gln)
    implicit none

    real(dp), intent(in)            :: a, x
    real(dp), intent(out), optional :: gln
    real(dp)                        :: gcf

    integer(i4b), parameter :: ITMAX = 100
    real(dp), parameter :: eps = epsilon(x), FPMIN=tiny(x)/EPS
    integer(i4b) :: i
    real(dp)     :: an, b, c, d, del, h
    if (x == 0.d0) then
       gcf = 0.d0
       return
    end if
    b = x + 1.d0 - a
    c = 1.d0 / FPMIN
    d = 1.d0 / b
    h = d
    do i = 1, ITMAX
       an = -i*(i-a)
       b = b + 2.d0
       d = an*d+b
       if (abs(d) < FPMIN) d = FPMIN
       c = b + an/c
       if (abs(c) < FPMIN) c = FPMIN
       d = 1.d0 / d
       del = d*c
       h = h*del
       if (abs(del-1.d0) <= eps) exit
    end do
    if (i > ITMAX) then
       write(*,*) 'a too large, ITMAX too small in gcf'
       stop
    end if
    if (present(gln)) then
       gln = gammln(a)
       gcf = exp(-x+a*log(x)-gln)*h
    else
       gcf = exp(-x+a*log(x)-gammln(a))*h
    end if

  end function gcf

  function erf2(x)
    real(dp), intent(in) :: x
    real(dp)             :: erf2

    erf2 = gammp(0.5d0, x**2)
    if (x < 0.d0) erf2 = -erf2
  end function erf2

  real(dp) function gammln(xx)
    implicit none

    real(dp), intent(in) :: xx

    real(dp) :: x, tmp
    real(dp) :: stp = 2.5066282746310005d0
    integer(i4b) :: i

    real(dp), dimension(6) :: coef

    coef(1) = 76.18009172947146d0
    coef(2) = -86.50532032941677d0
    coef(3) = 24.01409824083091d0
    coef(4) = -1.231739572450155d0
    coef(5) = 0.001208650973866179d0
    coef(6) = -0.000005395239384953d0

    x = xx
    tmp = x + 5.5d0
    tmp = (x+0.5d0) * log(tmp) - tmp
    gammln = 0.d0
    do i = 1, 6
       gammln = gammln + coef(i)/(x+real(i,dp))
    end do
    gammln = tmp+log(stp*(gammln + 1.000000000190015d0)/x)

  end function gammln

end module num_rec
