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
! template for routine SP/DP overloading for module alm_tools

! K M A P   : map kind                 either SP or DP
! K A L M C : alm kind (complex)       either SPC or DPC
! K A L M   : alm related and cl kind  either SP or DP
!
! K L O A D : suffixe of routine name, to be replaced by either s (SP) or d (DP)
  !**************************************************************************
  !
  !             ALM2MAP
  !
  !**************************************************************************
!     computes a map form its alm for the HEALPIX pixelisation
!     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
!                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
!
!     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
!
!     * the recurrence of Ylm is the standard one (cf Num Rec)
!     * the sum over m is done by FFT
  !=======================================================================
  subroutine alm2map_sc_KLOAD(nsmax, nlmax, nmmax, alm, map)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(KMAP),   intent(OUT), dimension(    0:12*nsmax**2-1) :: map

    integer(I4B) :: l, m, ith                          ! alm related
    integer(I8B) :: istart_south, istart_north, npix   ! map related
    integer(I4B) :: nrings, nphmx

    integer(I4B)                              :: par_lm
    real(DP)                                  :: b_er, b_ei, b_or, b_oi
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_lm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(i4b)                                :: ll, l_min, l_start
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_south')

    if (.not.do_openmp()) then
       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(0:1,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(lam_lm(0:nlmax), stat = status)
       call assert_alloc(status,code,'lam_lm')
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
    endif
    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    call init_rescale()
    map = 0.0 ! set the whole map to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, &
!$OMP    cth, sth, mfac, alm, b_north, b_south) &
!$OMP private(recfac, dalm, lam_lm, b_er, b_ei, b_or, b_oi, status, ll, l, m, ithl, l_min, l_start)

       if (do_openmp()) then
          allocate(recfac(0:1,0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac')
          allocate(dalm(0:1,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm')
          allocate(lam_lm(0:nlmax), stat = status)
          call assert_alloc(status,code,'lam_lm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(alm(1,ll,m))
          enddo
          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ! compute lam_lm(theta) for all l>=m
                call do_lam_lm(nlmax, m, cth(ithl), sth(ithl), mfac(m), recfac, lam_lm)
                b_er = 0.0_dp ; b_ei = 0.0_dp ; b_or = 0.0_dp ; b_oi = 0.0_dp
                ! do the product a_lm * Y_lm for all l>=m
                l_start = max(m, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start + 1
                do l = l_start, nlmax, 2
                   b_er = b_er + lam_lm(l) * dalm(0,l)
                   b_ei = b_ei + lam_lm(l) * dalm(1,l)
                enddo

                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start + 1
                do l = l_start, nlmax, 2
                   b_or = b_or + lam_lm(l) * dalm(0,l)
                   b_oi = b_oi + lam_lm(l) * dalm(1,l)
                enddo

                b_north(m,ithl) = cmplx(b_er + b_or, b_ei + b_oi, kind=DP) ! north: Even + Odd
                b_south(m,ithl) = cmplx(b_er - b_or, b_ei - b_oi, kind=DP) ! south: Even - Odd
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (recfac,dalm,lam_lm)
       endif
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_north, b_south, nph, startpix, kphi0, map) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          map(istart_north:istart_north+nphl-1) = ring(0:nphl-1)
       enddo ! loop on ithl
!$OMP end do
!$OMP do schedule(dynamic,1)
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_south = npix-startpix(ithl)-nphl
          ith  = ithl + lchk
          if (ith < nrings) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             map(istart_south:istart_south+nphl-1) = ring(0:nphl-1)
          endif
       enddo ! loop on ithl
!$OMP end do
       if (do_openmp()) then
          deallocate(ring)
       endif
!$OMP end parallel
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate (ring,recfac,lam_lm,dalm)
    endif
    deallocate(mfac)
    deallocate(b_north, b_south)
    return
  end subroutine alm2map_sc_KLOAD

!   !=======================================================================
!   subroutine alm2map_sc_KLOAD(nsmax, nlmax, nmmax, alm, map)
!     !=======================================================================
!     !     computes a map form its alm for the HEALPIX pixelisation
!     !      for the Temperature field
!     !=======================================================================
!     integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
!     complex(KALMC), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: alm
!     real(KMAP),   intent(OUT), dimension(    0:12*nsmax**2-1) :: map

!     integer(I4B) :: l, m, ith, scalel, scalem ! alm related
!     integer(I8B) :: istart_south, istart_north, npix   ! map related
!     integer(I4B) :: nrings, nphmx

!     integer(I4B)                              :: par_lm
!     real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
!     real(DP)                                  :: ovflow, unflow
!     complex(DPC), dimension(-1:1)             :: b_ns
!     real(DP),     dimension(:,:), allocatable :: dalm
!     real(DP),     dimension(:),   allocatable :: mfac
!     real(DP),     dimension(:,:), allocatable :: recfac

!     integer(i4b)                                :: ll, l_min
!     complex(DPC), dimension(:,:), allocatable :: b_north, b_south
!     real(DP),     dimension(:),   allocatable :: ring
!     integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
!     integer(I8B), dimension(0:SMAXCHK-1) :: startpix
!     integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
!     real(DP),     dimension(0:SMAXCHK-1) :: cth, sth

!     character(LEN=*), parameter :: code = 'ALM2MAP'
!     integer(I4B) :: status
!     !=======================================================================

!     ! Healpix definitions
!     nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!     npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!     nphmx  = 4*nsmax           ! maximum number of pixels/ring

!     !     --- allocates space for arrays ---
!     nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!     chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

!     allocate(mfac(0:nmmax),stat = status)
!     call assert_alloc(status,code,'mfac')

!     allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'b_north')

!     allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'b_south')

!     if (.not.do_openmp()) then
!        allocate(recfac(0:1,0:nlmax),stat = status)
!        call assert_alloc(status,code,'recfac')
!        allocate(dalm(0:1,0:nlmax), stat = status)
!        call assert_alloc(status,code,'dalm')
!        allocate(ring(0:nphmx-1),stat = status)
!        call assert_alloc(status,code,'ring')
!     endif
!     !     ------------ initiate variables and arrays ----------------

!     call gen_mfac(nmmax,mfac)

!     call init_rescale()
!     ovflow = rescale_tab(1)
!     unflow = rescale_tab(-1)
!     map = 0.0 ! set the whole map to zero

!     ! loop on chunks
!     do ichunk = 0, nchunks-1
!        lchk = ichunk * chunksize + 1
!        uchk = min(lchk+chunksize - 1, nrings)

!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           ! get pixel location information
!           call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!        enddo
!        !        -----------------------------------------------------
!        !        for each theta, and each m, computes
!        !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
!        !        ------------------------------------------------------
!        !        lambda_mm tends to go down when m increases (risk of underflow)
!        !        lambda_lm tends to go up   when l increases (risk of overflow)


!        b_north(:,:) = 0_dpc ! pad with zeros
!        b_south(:,:) = 0_dpc

! !$OMP parallel default(none) &
! !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, ovflow, unflow, &
! !$OMP    cth, sth, mfac, alm, b_north, b_south) &
! !$OMP private(recfac, dalm, b_ns, status, m, ithl, l_min,  &
! !$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
! !$OMP   cth_ring, l)

!        if (do_openmp()) then
!           allocate(recfac(0:1,0:nlmax),stat = status)
!           call assert_alloc(status,code,'recfac')
!           allocate(dalm(0:1,0:nlmax), stat = status)
!           call assert_alloc(status,code,'dalm')
!        endif

! !$OMP do schedule(dynamic,1)
!        do m = 0, nmmax
!           ! generate recursion factors (recfac) for Ylm of degree m
!           call gen_recfac(nlmax, m, recfac)

!           ! extract needed alm under memory and CPU efficient form
!           do ll = m, nlmax
!              dalm(0,ll) =  real(alm(1,ll,m),kind=dp)
!              dalm(1,ll) = aimag(alm(1,ll,m))
!           enddo
!           do ithl = 0, uchk - lchk
!              l_min = l_min_ylm(m, sth(ithl))
!              if (nlmax >= l_min) then ! skip calculations when Ylm too small
!                 ! determine lam_mm
!                 call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
               
!                 !           ---------- l = m ----------
!                 par_lm = 1  ! = (-1)^(l+m)

!                 if (m >= l_min) then
!                    lam_lm = corfac * lam_mm
!                    b_ns( 1) = cmplx(lam_lm * dalm(0,m), lam_lm * dalm(1,m), kind=DP)
!                    b_ns(-1) = 0.0_dpc
!                 else
!                    b_ns = 0.0_dpc
!                 endif

!                 !           ---------- l > m ----------
!                 lam_0 = 0.0_dp
!                 lam_1 = 1.0_dp
!                 scalel=0
!                 cth_ring = cth(ithl)
!                 lam_2 = cth_ring * lam_1 * recfac(0,m)
!                 do l = m+1, nlmax
!                    par_lm = - par_lm  ! = (-1)^(l+m)
!                    if (l >= l_min) then
!                       lam_lm = lam_2 * corfac * lam_mm                   
!                       b_ns(par_lm)  = b_ns(par_lm) &
!                            &        + cmplx(lam_lm * dalm(0,l), lam_lm * dalm(1,l), kind=DP)  ! increment even or odd
!                    endif

!                    lam_0 = lam_1 * recfac(1,l-1)
!                    lam_1 = lam_2
!                    lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
!                    if (abs(lam_2) > OVFLOW) then
!                       lam_1 = lam_1*UNFLOW
!                       lam_2 = lam_2*UNFLOW
!                       scalel= scalel + 1
!                       corfac = rescale_tab(max(scalel+scalem,RSMIN))
!                    elseif (abs(lam_2) < UNFLOW) then
!                       lam_1 = lam_1*OVFLOW
!                       lam_2 = lam_2*OVFLOW
!                       scalel= scalel - 1
!                       corfac = rescale_tab(max(scalel+scalem,RSMIN))
!                    endif
                   
!                 enddo ! loop on l

!                 b_north(m,ithl) = b_ns(1) + b_ns(-1)
!                 b_south(m,ithl) = b_ns(1) - b_ns(-1)
!              endif ! test on nlmax
!           enddo ! loop on rings (ithl)
!        enddo ! loop on m
! !$OMP end do
!        if (do_openmp()) then
!           deallocate (recfac,dalm)
!        endif
! !$OMP end parallel

! !$OMP parallel default(none) &
! !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
! !$OMP      lchk, uchk, b_north, b_south, nph, startpix, kphi0, map) &
! !$OMP  private(ithl, nphl, istart_north, istart_south, &
! !$OMP      ith, ring, status)
!        if (do_openmp()) then
!           allocate(ring(0:nphmx-1),stat = status)
!           call assert_alloc(status,code,'ring')
!        endif
! !$OMP do schedule(dynamic,1)
!        do ithl = 0, uchk - lchk
!           nphl = nph(ithl)
!           istart_north = startpix(ithl)
!           istart_south = npix-istart_north-nphl
!           ith  = ithl + lchk
!           !        ---------------------------------------------------------------
!           !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
!           !        ---------------------------------------------------------------
!           call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
!           map(istart_north:istart_north+nphl-1) = ring(0:nphl-1)

!           if (ith < nrings) then
!              call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
!              map(istart_south:istart_south+nphl-1) = ring(0:nphl-1)
!           endif
!        enddo ! loop on ithl
! !$OMP end do
!        if (do_openmp()) then
!           deallocate(ring)
!        endif
! !$OMP end parallel
!     enddo    ! loop on chunks

!     !     --------------------
!     !     free memory and exit
!     !     --------------------
!     if (.not.do_openmp()) then
!        deallocate (ring,recfac,dalm)
!     endif
!     deallocate(mfac)
!     deallocate(b_north, b_south)
!     return
!   end subroutine alm2map_sc_KLOAD


  !=======================================================================
  subroutine alm2map_sc_pre_KLOAD(nsmax, nlmax, nmmax, alm, map, plm)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field, with precomputed Ylm
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(KMAP),   intent(OUT), dimension(    0:12*nsmax**2-1) :: map
    real(DP),     intent(IN),  dimension(0:)                  :: plm

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx

    integer(I8B) :: n_lm, n_plm, i_mm
    complex(DPC), dimension(-1:1)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(i4b)                              :: ll, l_min, l_start
    real(DP)                                  :: cth
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_south')

    if (.not. do_openmp()) then
       allocate(dalm(0:1,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
    endif
    !     ------------ initiate variables and arrays ----------------

    map = 0.0 ! set the whole map to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------

       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, sth, lchk, uchk, plm, alm, b_north, b_south, n_lm) &
!$OMP private(dalm, b_ns, status, m, ith, ithl, l_min, l_start, l, ll, i_mm)

       if (do_openmp()) then
          allocate(dalm(0:1,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(alm(1,ll,m))
          enddo
          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------

                if (m >= l_min) then
                   b_ns( 1) = cmplx(plm(i_mm) * dalm(0,m), plm(i_mm) * dalm(1,m), kind=DP)
                   b_ns(-1) = 0.0_dpc
                else
                   b_ns = 0.0_dpc
                endif

                !           ---------- l > m ----------
                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(-1) = b_ns(-1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l), kind=DP)
                enddo ! loop on l

                l_start = max(m+2, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(1)  = b_ns(1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l),kind=DP) 
                enddo ! loop on l

                b_north(m,ithl) = b_ns(1) + b_ns(-1)
                b_south(m,ithl) = b_ns(1) - b_ns(-1)
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (dalm)
       endif
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_north, b_south, nph, startpix, kphi0, map) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          map(istart_north:istart_north+nphl-1) = ring(0:nphl-1)

          if (ith < nrings) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             map(istart_south:istart_south+nphl-1) = ring(0:nphl-1)
          endif
       enddo ! loop on ithl
!$OMP end do
       if (do_openmp()) then
          deallocate(ring)
       endif
!$OMP end parallel
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate (ring,dalm)
    endif
    deallocate(b_north, b_south)
    return
  end subroutine alm2map_sc_pre_KLOAD

  !=======================================================================
  subroutine alm2map_pol_KLOAD(nsmax, nlmax, nmmax, alm_TGC, map_TQU)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), intent(IN),  dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(KMAP),   intent(OUT), dimension(0:12*nsmax**2-1,1:3) :: map_TQU

    integer(I4B) :: l, m, ith                    ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, npix, nphmx

    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: recfac, dalm, lam_lm
    real(DP),     dimension(:),   allocatable :: lam_fact, mfac

    integer(i4b)                                :: ll, l_min, l_start
    real(DP),     dimension(:),     allocatable :: normal_l
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
!    real(sp) :: time0, time1, time2, time3, time4, tt1, tt2

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    if (.not. do_openmp()) then
       allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac, dalm & lam_fact')
       allocate(lam_lm(1:3,0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_lm')
       allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
       call assert_alloc(status,code,'ring & bsub')
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    call init_rescale()
    map_TQU = 0.0 ! set the whole map to zero
    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------


       b_TQU = 0_dpc ! pad with zeros

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, &
!$OMP    rescale_tab, &
!$OMP    cth, sth, mfac, normal_l, alm_TGC, b_TQU) &
!$OMP private(recfac, dalm, lam_fact, lam_lm, b_ns, status, &
!$OMP   m, l, ll, ithl, l_min, l_start)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac, dalm & lam_fact')
          allocate(lam_lm(1:3,0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_lm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm_TGC(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(alm_TGC(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(alm_TGC(2,ll,m),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(alm_TGC(2,ll,m))        *normal_l(ll)
             dalm(4,ll) =  real(alm_TGC(3,ll,m),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(alm_TGC(3,ll,m))        *normal_l(ll)
          enddo
          do ithl = 0, uchk - lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ! compute lam_lm(p,theta) for all l>=m
                call do_lam_lm_pol(nlmax, m, cth(ithl), sth(ithl), mfac(m), recfac, lam_fact, lam_lm)
                   !      T =   alm_T * Ylm
                   !      Q = - alm_E * Wlm - i alm_B * Xlm
                   !      U = i alm_E * Xlm -   alm_B * Wlm

                b_ns = 0.0_dp
                l_start = max(m, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start + 1
                do l = l_start, nlmax, 2
                   b_ns(3: 4) = b_ns(3: 4) + lam_lm(1,l) * dalm(0:1,l) ! T even
                   b_ns(5: 8) = b_ns(5: 8) - lam_lm(2,l) * dalm(2:5,l) ! Q, U  even
                   b_ns(-1)   = b_ns(-1)   + lam_lm(3,l) * dalm(5,l) ! Q odd
                   b_ns( 0)   = b_ns( 0)   - lam_lm(3,l) * dalm(4,l)
                   b_ns( 1)   = b_ns( 1)   - lam_lm(3,l) * dalm(3,l) ! U odd
                   b_ns( 2)   = b_ns( 2)   + lam_lm(3,l) * dalm(2,l)
                enddo

                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start + 1
                do l = l_start, nlmax, 2
                   b_ns(-3:-2) = b_ns(-3:-2) + lam_lm(1,l) * dalm(0:1,l) ! T odd
                   b_ns(-1: 2) = b_ns(-1: 2) - lam_lm(2,l) * dalm(2:5,l) ! Q, U  odd
                   b_ns( 5)    = b_ns( 5)    + lam_lm(3,l) * dalm(5,l) ! Q even
                   b_ns( 6)    = b_ns( 6)    - lam_lm(3,l) * dalm(4,l)
                   b_ns( 7)    = b_ns( 7)    - lam_lm(3,l) * dalm(3,l) ! U even
                   b_ns( 8)    = b_ns( 8)    + lam_lm(3,l) * dalm(2,l)
                enddo
               
                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) deallocate (recfac,dalm, lam_fact, lam_lm)
!$OMP end parallel


!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_TQU, nph, startpix, kphi0, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, bsub, i, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
          call assert_alloc(status,code,'ring & bsub')
       endif

!$OMP do schedule(dynamic,1)
!        call wall_clock_time(time2)
!        tt2 = tt2 + time2 - time3
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             map_TQU(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo

          if (ith < nrings) then
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                map_TQU(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
!$OMP end do
       if (do_openmp()) deallocate(ring, bsub)
!$OMP end parallel

!        call wall_clock_time(time3)
!        tt1 = tt1 + time3 - time2

    enddo    ! loop on chunks
!     print*,'FFT:',tt1
!     print*,'the rest:',tt2

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate(recfac, dalm, lam_fact, lam_lm, ring, bsub)
    endif
    deallocate(mfac, normal_l)
    deallocate(b_TQU)
    return
  end subroutine alm2map_pol_KLOAD
!   !=======================================================================
!   subroutine alm2map_pol_KLOAD(nsmax, nlmax, nmmax, alm_TGC, map_TQU)
!     !=======================================================================
!     !     computes a map form its alm for the HEALPIX pixelisation
!     !      for the Temperature+Polarization field
!     !=======================================================================
!     integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
!     complex(KALMC), intent(IN),  dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
!     real(KMAP),   intent(OUT), dimension(0:12*nsmax**2-1,1:3) :: map_TQU

!     integer(I4B) :: l, m, ith, scalel, scalem ! alm related
!     integer(I4B) :: istart_south, istart_north   ! map related
!     integer(I4B) :: nrings, npix, nphmx

!     integer(I4B)                              :: par_lm
!     real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
!     real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
!     real(DP)                 :: fm2, fl, flm1
!     real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
!     real(DP)                 :: ovflow, unflow
!     real(DP),     dimension(-3:8)             :: b_ns
!     real(DP),     dimension(:,:), allocatable :: recfac, dalm
!     real(DP),     dimension(:),   allocatable :: lam_fact
!     real(DP),     dimension(:), allocatable :: mfac

!     integer(i4b)                                :: ll, l_min
!     real(DP),     dimension(:),     allocatable :: normal_l
!     complex(DPC), dimension(:,:,:), allocatable :: b_TQU
!     complex(DPC), dimension(:),     allocatable :: bsub
!     real(DP),     dimension(:),     allocatable :: ring
!     integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
!     integer(I8B), dimension(0:SMAXCHK) :: startpix
!     integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
!     real(DP),     dimension(0:SMAXCHK) :: cth, sth
!     real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2
! !    real(sp) :: time0, time1, time2, time3, time4, tt1, tt2

!     character(LEN=*), parameter :: code = 'ALM2MAP'
!     integer(I4B) :: status
!     !=======================================================================

!     ! Healpix definitions
!     nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!     npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!     nphmx  = 4*nsmax           ! maximum number of pixels/ring

!     !     --- allocates space for arrays ---
!     nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!     chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

!     allocate(mfac(0:nmmax),stat = status)
!     call assert_alloc(status,code,'mfac')

!     allocate(normal_l(0:nlmax),stat = status)
!     call assert_alloc(status,code,'normal_l')

!     allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'b_TQU')

!     if (.not. do_openmp()) then
!        allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
!        call assert_alloc(status,code,'recfac, dalm & lam_fact')
!        allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
!        call assert_alloc(status,code,'ring & bsub')
!     endif

!     !     ------------ initiate variables and arrays ----------------

!     call gen_mfac(nmmax,mfac)

!     call init_rescale()
!     ovflow = rescale_tab(1)
!     unflow = rescale_tab(-1)
!     map_TQU = 0.0 ! set the whole map to zero
!     ! generate Polarization normalisation
!     call gen_normpol(nlmax, normal_l)

! !     call wall_clock_time(time3)
! !     tt1 = 0 ; tt2 = 0
!     ! loop on chunks
!     do ichunk = 0, nchunks-1
!        lchk = ichunk * chunksize + 1
!        uchk = min(lchk+chunksize - 1, nrings)

!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           ! get pixel location information
!           call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!           one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
!             c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
!        enddo
!        !        -----------------------------------------------------
!        !        for each theta, and each m, computes
!        !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
!        !        ------------------------------------------------------
!        !        lambda_mm tends to go down when m increases (risk of underflow)
!        !        lambda_lm tends to go up   when l increases (risk of overflow)


!        b_TQU = 0_dpc ! pad with zeros

! !$OMP parallel default(none) &
! !$OMP shared(nlmax, nmmax, lchk, uchk, &
! !$OMP    rescale_tab, ovflow, unflow, &
! !$OMP    cth, sth, mfac, normal_l, alm_TGC, b_TQU, one_on_s2, c_on_s2) &
! !$OMP private(recfac, dalm, lam_fact, b_ns, status, &
! !$OMP   m, ll, fm2, normal_m, ithl, l_min, &
! !$OMP   scalem, scalel, corfac, par_lm, &
! !$OMP   lam_mm, lam_lm, lam_lm1m, lambda_w, lambda_x, lam_0, lam_1, lam_2, &
! !$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
! !$OMP   l, fl, flm1)

!        if (do_openmp()) then
!           ! allocate thread safe arrays
!           allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
!           call assert_alloc(status,code,'recfac, dalm & lam_fact')
!        endif

! !$OMP do schedule(dynamic,1)
!        do m = 0, nmmax
!           ! generate recursion factors (recfac) for Ylm of degree m
!           call gen_recfac(nlmax, m, recfac)
!           ! generate Ylm relation factor for degree m
!           call gen_lamfac(nlmax, m, lam_fact)
!           fm2 = real(m * m, kind = DP)
!           normal_m = (2.0_dp * m) * (1 - m)

!           ! extract needed alm under memory and CPU efficient form
!           do ll = m, nlmax
!              dalm(0,ll) =  real(alm_TGC(1,ll,m),kind=dp) ! T, real
!              dalm(1,ll) = aimag(alm_TGC(1,ll,m))         ! T, imag
!              dalm(2,ll) =  real(alm_TGC(2,ll,m),kind=dp)*normal_l(ll) ! G, real
!              dalm(3,ll) = aimag(alm_TGC(2,ll,m))        *normal_l(ll)
!              dalm(4,ll) =  real(alm_TGC(3,ll,m),kind=dp)*normal_l(ll) ! C, real
!              dalm(5,ll) = aimag(alm_TGC(3,ll,m))        *normal_l(ll)
!           enddo
!           do ithl = 0, uchk - lchk
!              l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
!              if (nlmax >= l_min) then
!                 ! determine lam_mm
!                 call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

!                 !           ---------- l = m ----------
!                 par_lm = 3  ! = 3 * (-1)^(l+m)
!                 lam_lm = corfac * lam_mm !Actual lam_mm 

!                 if (m >= l_min) then ! skip Ymm if too small
!                    lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
!                    lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

!                    !      T =   alm_T * Ylm
!                    !      Q = - alm_E * Wlm - i alm_B * Xlm
!                    !      U = i alm_E * Xlm -   alm_B * Wlm
!                    b_ns(-3:-2) = 0.0_dp               ! T odd
!                    b_ns(-1: 0) =   lambda_x * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
!                    b_ns( 1: 2) =   lambda_x * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
! !                    b_ns(-1) =   lambda_x * dalm(5,m) ! Q odd
! !                    b_ns( 0) =  -lambda_x * dalm(4,m) ! Q odd
! !                    b_ns( 1) =  -lambda_x * dalm(3,m) ! U odd
! !                    b_ns( 2) =   lambda_x * dalm(2,m) ! U odd

!                    b_ns( 3: 4) =   lam_lm   * dalm(0:1,m) ! T even
!                    b_ns( 5: 6) = - lambda_w * dalm(2:3,m) ! Q even
!                    b_ns( 7: 8) = - lambda_w * dalm(4:5,m) ! U even
!                 else
!                    b_ns = 0.0_dp
!                 endif

!                 !           ---------- l > m ----------
!                 lam_0 = 0.0_dp
!                 lam_1 = 1.0_dp
!                 scalel=0
!                 cth_ring = cth(ithl)
!                 lam_2 = cth_ring * lam_1 * recfac(0,m)
!                 fm_on_s2     =      m * one_on_s2(ithl)
!                 two_on_s2    = 2.0_dp * one_on_s2(ithl)
!                 two_cth_ring = 2.0_dp * cth_ring
!                 b_w          =  c_on_s2(ithl) 
!                 do l = m+1, nlmax
!                    par_lm   = - par_lm  ! = 3 * (-1)^(l+m)
!                    lam_lm1m = lam_lm * lam_fact(l) ! must be incremented even if not used
!                    lam_lm   = lam_2 * corfac * lam_mm
!                    if (l >= l_min) then

!                       fl = real(l, kind = DP)
!                       flm1 = fl - 1.0_dp
!                       a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
!                       a_x =  two_cth_ring * flm1
!                       lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
!                       lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

!                       b_ns(par_lm:  par_lm+1) = b_ns(par_lm:  par_lm+1) + lam_lm   * dalm(0:1,l) ! T even or odd
!                       b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - lambda_w * dalm(2:5,l) ! Q, U  even or odd
!                       b_ns(2-par_lm) = b_ns(2-par_lm) + lambda_x * dalm(5,l) ! Q odd (or even)
!                       b_ns(3-par_lm) = b_ns(3-par_lm) - lambda_x * dalm(4,l)
!                       b_ns(4-par_lm) = b_ns(4-par_lm) - lambda_x * dalm(3,l) ! U odd (or even)
!                       b_ns(5-par_lm) = b_ns(5-par_lm) + lambda_x * dalm(2,l)
!                    endif

!                    lam_0 = lam_1 * recfac(1,l-1)
!                    lam_1 = lam_2
!                    lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
!                    if (abs(lam_2) > OVFLOW) then
!                       lam_1 = lam_1*UNFLOW
!                       lam_2 = lam_2*UNFLOW
!                       scalel= scalel + 1
!                       corfac = rescale_tab(max(scalel+scalem,RSMIN))
!                    elseif (abs(lam_2) < UNFLOW) then
!                       lam_1 = lam_1*OVFLOW
!                       lam_2 = lam_2*OVFLOW
!                       scalel= scalel - 1
!                       corfac = rescale_tab(max(scalel+scalem,RSMIN))
!                    endif

!                 enddo ! loop on l

!                 b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
!                 b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
!                 b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
!                 b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
!                 b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
!                 b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
!              endif ! test on nlmax
!           enddo ! loop on rings (ithl)
!        enddo ! loop on m
! !$OMP end do
!        if (do_openmp()) deallocate (recfac,dalm, lam_fact)
! !$OMP end parallel


! !$OMP parallel default(none) &
! !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
! !$OMP      lchk, uchk, b_TQU, nph, startpix, kphi0, map_TQU) &
! !$OMP  private(ithl, nphl, istart_north, istart_south, &
! !$OMP      ith, ring, bsub, i, status)

!        if (do_openmp()) then
!           allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
!           call assert_alloc(status,code,'ring & bsub')
!        endif

! !$OMP do schedule(dynamic,1)
! !        call wall_clock_time(time2)
! !        tt2 = tt2 + time2 - time3
!        do ithl = 0, uchk - lchk
!           nphl = nph(ithl)
!           istart_north = startpix(ithl)
!           istart_south = npix-istart_north-nphl
!           ith  = ithl + lchk
!           !        ---------------------------------------------------------------
!           !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
!           !        ---------------------------------------------------------------
!           do i=1,3
!              bsub = b_TQU(i,:,ithl)
!              call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
!              map_TQU(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
!           enddo

!           if (ith < nrings) then
!              do i=1,3
!                 bsub = b_TQU(3+i,:,ithl)
!                 call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
!                 map_TQU(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
!              enddo
!           endif
!        enddo ! loop on ithl
! !$OMP end do
!        if (do_openmp()) deallocate(ring, bsub)
! !$OMP end parallel

! !        call wall_clock_time(time3)
! !        tt1 = tt1 + time3 - time2

!     enddo    ! loop on chunks
! !     print*,'FFT:',tt1
! !     print*,'the rest:',tt2

!     !     --------------------
!     !     free memory and exit
!     !     --------------------
!     if (.not. do_openmp()) then
!        deallocate(recfac, dalm, lam_fact, ring, bsub)
!     endif
!     deallocate(mfac, normal_l)
!     deallocate(b_TQU)
!     return
!   end subroutine alm2map_pol_KLOAD

  !=======================================================================
  subroutine alm2map_pol_pre1_KLOAD(nsmax, nlmax, nmmax, alm_TGC, map_TQU, plm)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), intent(IN),  dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(KMAP),   intent(OUT), dimension(0:12*nsmax**2-1,1:3) :: map_TQU
    real(DP),     intent(IN),  dimension(0:) :: plm

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, n_plm, i_mm
    integer(I4B) :: nrings, nphmx

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_lm, cth_ring, lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                                :: ll, l_min
    real(DP),     dimension(:),     allocatable :: normal_l
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    if (.not. do_openmp()) then
       allocate(dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm & lam_fact')
       allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
       call assert_alloc(status,code,'ring & bsub')
    endif

    !     ------------ initiate variables and arrays ----------------

    map_TQU = 0.0 ! set the whole map to zero
    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_TQU = 0_dpc ! pad with zeros

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, n_lm, ith, &
!$OMP    cth, sth, normal_l, alm_TGC, b_TQU, plm, one_on_s2, c_on_s2) &
!$OMP private(dalm, lam_fact, b_ns, status, i_mm, &
!$OMP   m, ll, fm2, normal_m, ithl, l_min, &
!$OMP   par_lm, lam_lm, lam_lm1m, lambda_w, lambda_x, &
!$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
!$OMP   l, fl, flm1)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(dalm(0:5,0:nlmax), lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'dalm & lam_fact')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm_TGC(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(alm_TGC(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(alm_TGC(2,ll,m),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(alm_TGC(2,ll,m))        *normal_l(ll)
             dalm(4,ll) =  real(alm_TGC(3,ll,m),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(alm_TGC(3,ll,m))        *normal_l(ll)
          enddo
          do ithl = 0, uchk - lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = plm(i_mm)

                if (m >= l_min) then ! skip Ymm if too small
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   !      T =   alm_T * Ylm
                   !      Q = - alm_E * Wlm - i alm_B * Xlm
                   !      U = i alm_E * Xlm -   alm_B * Wlm
                   b_ns(-3:-2) = 0.0_dp               ! T odd
                   b_ns(-1: 0) =   lambda_x * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
                   b_ns( 1: 2) =   lambda_x * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
                   b_ns( 3: 4) =   lam_lm   * dalm(0:1,m) ! T even
                   b_ns( 5: 6) = - lambda_w * dalm(2:3,m) ! Q even
                   b_ns( 7: 8) = - lambda_w * dalm(4:5,m) ! U even
                else
                   b_ns = 0.0_dp
                endif

                !           ---------- l > m ----------
                cth_ring = cth(ithl)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                   if (l >= l_min) then
                      lam_lm1m = plm(i_mm+l-m-1) * lam_fact(l)
                      lam_lm   = plm(i_mm+l-m)

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                      b_ns(par_lm:  par_lm+1) = b_ns(par_lm:  par_lm+1) + lam_lm   * dalm(0:1,l) ! T even or odd
                      b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - lambda_w * dalm(2:5,l) ! Q, U  even or odd
                      b_ns(2-par_lm) = b_ns(2-par_lm) + lambda_x * dalm(5,l) ! Q odd (or even)
                      b_ns(3-par_lm) = b_ns(3-par_lm) - lambda_x * dalm(4,l)
                      b_ns(4-par_lm) = b_ns(4-par_lm) - lambda_x * dalm(3,l) ! U odd (or even)
                      b_ns(5-par_lm) = b_ns(5-par_lm) + lambda_x * dalm(2,l)
                    endif

                enddo ! loop on l

                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) deallocate (dalm, lam_fact)
!$OMP end parallel


!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_TQU, nph, startpix, kphi0, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, bsub, i, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
          call assert_alloc(status,code,'ring & bsub')
       endif

!$OMP do schedule(dynamic,1)
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             map_TQU(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo

          if (ith < nrings) then
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                map_TQU(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
!$OMP end do
       if (do_openmp()) deallocate(ring, bsub)
!$OMP end parallel

    enddo    ! loop on chunks
    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate(dalm, lam_fact, ring, bsub)
    endif
    deallocate(normal_l)
    deallocate(b_TQU)
    return
  end subroutine alm2map_pol_pre1_KLOAD

  !=======================================================================
  subroutine alm2map_pol_pre2_KLOAD(nsmax, nlmax, nmmax, alm_TGC, map_TQU, plm)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature+Polarization field
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), intent(IN),  dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(KMAP),   intent(OUT), dimension(0:12*nsmax**2-1,1:3) :: map_TQU
    real(DP),     intent(IN),  dimension(0:,1:) :: plm

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx

    integer(I8B) :: n_lm, n_plm, i_mm, i_up
    real(DP),     dimension(-3:8)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm, plm_sub

    integer(i4b)                                :: ll, l_min, l_start
    complex(DPC), dimension(:,:,:), allocatable :: b_TQU
    complex(DPC), dimension(:),     allocatable :: bsub
    real(DP),     dimension(:),     allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i
    real(DP)                           :: cth
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    real(DP),     dimension(0:SMAXCHK) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(b_TQU(1:6, 0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_TQU')

    if (.not. do_openmp()) then
       allocate(dalm(0:5,0:nlmax),plm_sub(1:3,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm & plm_sub')
       allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
       call assert_alloc(status,code,'ring & bsub')
    endif

    !     ------------ initiate variables and arrays ----------------
    map_TQU = 0.0 ! set the whole map to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)


       b_TQU = 0_dpc ! pad with zeros

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, cth, sth, alm_TGC, b_TQU, plm, n_lm, ith) &
!$OMP private(dalm, b_ns, plm_sub, status, m, ll, ithl, l_min, &
!$OMP   l_start, l, i_mm, i_up)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(dalm(0:5,0:nlmax), plm_sub(1:3,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm & plm_sub')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm_TGC(1,ll,m),kind=dp) ! T, real
             dalm(1,ll) = aimag(alm_TGC(1,ll,m))         ! T, imag
             dalm(2,ll) =  real(alm_TGC(2,ll,m),kind=dp) ! G, real
             dalm(3,ll) = aimag(alm_TGC(2,ll,m))        
             dalm(4,ll) =  real(alm_TGC(3,ll,m),kind=dp) ! C, real
             dalm(5,ll) = aimag(alm_TGC(3,ll,m))        
          enddo
          do ithl = 0, uchk - lchk
             l_min =  l_min_ylm(m, sth(ithl)) ! lowest l of non-negligible Ylm
             if (nlmax >= l_min) then
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith
                i_up = i_mm + nlmax - m
                plm_sub(1,m:nlmax) = plm(i_mm:i_up,1)
                plm_sub(2,m:nlmax) = plm(i_mm:i_up,2)
                plm_sub(3,m:nlmax) = plm(i_mm:i_up,3)
                !           ---------- l = m ----------

                if (m >= l_min) then ! skip Ymm if too small
                   !      T =   alm_T * Ylm
                   !      Q = - alm_E * Wlm - i alm_B * Xlm
                   !      U = i alm_E * Xlm -   alm_B * Wlm
                   b_ns(-3:-2) = 0.0_dp               ! T odd
                   b_ns(-1: 0) =   plm_sub(3,m) * (/   dalm(5,m), - dalm(4,m) /) ! Q odd
                   b_ns( 1: 2) =   plm_sub(3,m) * (/ - dalm(3,m),   dalm(2,m) /) ! U odd
                   b_ns( 3: 4) =   plm_sub(1,m)   * dalm(0:1,m) ! T even
                   b_ns( 5: 6) = - plm_sub(2,m) * dalm(2:3,m) ! Q even
                   b_ns( 7: 8) = - plm_sub(2,m) * dalm(4:5,m) ! U even
                else
                   b_ns = 0.0_dp
                endif

                !           ---------- l > m ----------
!                 do l = m+1, nlmax
!                    par_lm = - par_lm  ! = 3 * (-1)^(l+m)
!                    if (l >= l_min) then
!                       b_ns(par_lm:  par_lm+1) = b_ns(par_lm:  par_lm+1) + plm_sub(1,l) * dalm(0:1,l) ! T even or odd
!                       b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - plm_sub(2,l) * dalm(2:5,l) ! Q, U  even or odd
!                       b_ns(2-par_lm) = b_ns(2-par_lm) + plm_sub(3,l) * dalm(5,l) ! Q odd (or even)
!                       b_ns(3-par_lm) = b_ns(3-par_lm) - plm_sub(3,l) * dalm(4,l)
!                       b_ns(4-par_lm) = b_ns(4-par_lm) - plm_sub(3,l) * dalm(3,l) ! U odd (or even)
!                       b_ns(5-par_lm) = b_ns(5-par_lm) + plm_sub(3,l) * dalm(2,l)
!                    endif
                   
                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(-3:-2) = b_ns(-3:-2) + plm_sub(1,l) * dalm(0:1,l) ! T odd
                   b_ns(-1: 2) = b_ns(-1: 2) - plm_sub(2,l) * dalm(2:5,l) ! Q, U  odd
                   b_ns(5)     = b_ns(5)     + plm_sub(3,l) * dalm(5,l) ! Q even
                   b_ns(6)     = b_ns(6)     - plm_sub(3,l) * dalm(4,l)
                   b_ns(7)     = b_ns(7)     - plm_sub(3,l) * dalm(3,l) ! U even
                   b_ns(8)     = b_ns(8)     + plm_sub(3,l) * dalm(2,l)
                enddo ! loop on l

                l_start = max(m+2, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(3:4) = b_ns(3:4) + plm_sub(1,l) * dalm(0:1,l) ! T even 
                   b_ns(5:8) = b_ns(5:8) - plm_sub(2,l) * dalm(2:5,l) ! Q, U  even
                   b_ns(-1)  = b_ns(-1)  + plm_sub(3,l) * dalm(5,l) ! Q odd
                   b_ns( 0)  = b_ns( 0)  - plm_sub(3,l) * dalm(4,l)
                   b_ns( 1)  = b_ns( 1)  - plm_sub(3,l) * dalm(3,l) ! U odd
                   b_ns( 2)  = b_ns( 2)  + plm_sub(3,l) * dalm(2,l)
                enddo ! loop on l

                b_TQU(1,m,ithl) = cmplx(b_ns(3) + b_ns(-3), b_ns(4) + b_ns(-2), kind = DP) ! T north
                b_TQU(4,m,ithl) = cmplx(b_ns(3) - b_ns(-3), b_ns(4) - b_ns(-2), kind = DP) ! T south
                b_TQU(2,m,ithl) = cmplx(b_ns(5) + b_ns(-1), b_ns(6) + b_ns( 0), kind = DP) ! Q n
                b_TQU(5,m,ithl) = cmplx(b_ns(5) - b_ns(-1), b_ns(6) - b_ns( 0), kind = DP) ! Q s
                b_TQU(3,m,ithl) = cmplx(b_ns(7) + b_ns( 1), b_ns(8) + b_ns( 2), kind = DP) ! U n
                b_TQU(6,m,ithl) = cmplx(b_ns(7) - b_ns( 1), b_ns(8) - b_ns( 2), kind = DP) ! U s
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) deallocate (dalm)
!$OMP end parallel


!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_TQU, nph, startpix, kphi0, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, bsub, i, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),bsub(0:nlmax),stat = status)
          call assert_alloc(status,code,'ring & bsub')
       endif

!$OMP do schedule(dynamic,1)
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk
          !        ---------------------------------------------------------------
          !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta) by FFT
          !        ---------------------------------------------------------------
          do i=1,3
             bsub = b_TQU(i,:,ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             map_TQU(istart_north:istart_north+nphl-1, i) = ring(0:nphl-1)
          enddo

          if (ith < nrings) then
             do i=1,3
                bsub = b_TQU(3+i,:,ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bsub,nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
                map_TQU(istart_south:istart_south+nphl-1, i) = ring(0:nphl-1)
             enddo
          endif
       enddo ! loop on ithl
!$OMP end do
       if (do_openmp()) deallocate(ring, bsub)
!$OMP end parallel

    enddo    ! loop on chunks
    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate(dalm, ring, bsub)
    endif
    deallocate(b_TQU)
    return
  end subroutine alm2map_pol_pre2_KLOAD


  !=======================================================================
  subroutine alm2map_sc_der_KLOAD(nsmax, nlmax, nmmax, alm, map, der1, der2)
    !=======================================================================
    !     computes a map and its 1st and 2nd derivative 
    !       from its alm for the HEALPIX pixelisation
    !      for the Temperature field
    !
    ! der1 = dT/d(theta), dT/d(phi)
    ! der2 = d^2T/d(theta^2), d^2T/d(theta)/d(phi), d^2T/d(phi^2)
    !
    ! derivatives of Ylm are obtained from Qu.Th. of Ang.Mom. Varshalovich et al
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), INTENT(IN),  dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1)     :: map
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1,1:2) :: der1
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1,1:3), optional :: der2

    integer(I4B) :: l, m, ith, scalel, scalem ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, npix, nphmx

    integer(I4B)                              :: par_lm
    real(DP)             :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                  :: ovflow, unflow
    real(DP)                  :: f2m, fm2, fm, fllp1, lam_lm1m
    real(dp)                  :: cotanth, one_on_s1, one_on_s2
    real(DP), dimension(-1:2) :: b_ns
    real(DP), dimension(-1:2) :: b_ns_p, b_ns_t
    real(DP), dimension(-1:2) :: b_ns_pp, b_ns_tt, b_ns_tp
    real(DP), dimension( 0:1) :: factor, dfdt, d2fdt2

    integer(i4b)                            :: ll, l_min
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north,    b_south
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north_t,  b_south_t
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north_p,  b_south_p
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north_tt, b_south_tt
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north_tp, b_south_tp
    complex(DPC), dimension(:,:), ALLOCATABLE :: b_north_pp, b_south_pp
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    real(DP), dimension(:,:), ALLOCATABLE :: recfac, dalm
    real(DP), dimension(:),   ALLOCATABLE :: lam_fact, mfac
    real(DP), dimension(:),   ALLOCATABLE :: ring

    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph,kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth, one_on_s

    character(len=*), parameter :: code = 'ALM2MAP_DER'
    logical(lgt) :: do_d2
    integer(I4B) :: status
    !=======================================================================
    do_d2 = (present(der2)) ! do 2nd derivative

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(b_north(0:nmmax,0:chunksize-1), b_north_t(0:nmmax,0:chunksize-1), &
         &   b_north_p(0:nmmax,0:chunksize-1),stat = status)
    allocate(b_south(0:nmmax,0:chunksize-1), b_south_t(0:nmmax,0:chunksize-1), &
         &   b_south_p(0:nmmax,0:chunksize-1),stat = status)
    if (do_d2) then
       allocate(b_north_tt(0:nmmax,0:chunksize-1),b_north_tp(0:nmmax,0:chunksize-1), &
            &   b_north_pp(0:nmmax,0:chunksize-1),                    stat = status)
       allocate(b_south_tt(0:nmmax,0:chunksize-1),b_south_tp(0:nmmax,0:chunksize-1), &
            &   b_south_pp(0:nmmax,0:chunksize-1),                    stat = status)
    endif
    call assert_alloc(status,code,'b_north & b_south')

    if (.not. do_openmp()) then
       allocate(recfac(0:1,0:nlmax), dalm(0:1,0:nlmax), lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac, dalm & lam_fact')
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
    endif

    !     ------------ initiate arrays ----------------
    call gen_mfac(nmmax,mfac)

    call init_rescale()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    map = 0.0 ! set the whole map to zero
    map = 0.0
    der1 = 0.0
    if (do_d2) der2 = 0.0

    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s(ithl)  = 1.0_dp / sth(ithl)
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       b_north    = 0_dpc ; b_south = 0_dpc    ! pad with zeros
       b_north_t  = 0_dpc ; b_south_t  = 0_dpc
       b_north_p  = 0_dpc ; b_south_p  = 0_dpc
       if (do_d2) then
          b_north_tt = 0_dpc ; b_south_tt = 0_dpc
          b_north_tp = 0_dpc ; b_south_tp = 0_dpc
          b_north_pp = 0_dpc ; b_south_pp = 0_dpc
       endif

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, &
!$OMP    rescale_tab, ovflow, unflow, &
!$OMP    cth, sth, mfac, alm, one_on_s,  &
!$OMP      b_north, b_south, b_north_t, b_south_t, &
!$OMP      b_north_p,  b_south_p,  b_north_tt, b_south_tt, &
!$OMP      b_north_tp, b_south_tp, b_north_pp, b_south_pp, do_d2 ) &
!$OMP private(recfac, dalm, lam_fact, status, &
!$OMP   m, ll, fm, f2m, fm2, ithl, l_min, &
!$OMP   scalem, scalel, corfac, par_lm, &
!$OMP   lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2, &
!$OMP   cth_ring, one_on_s1, one_on_s2, cotanth, factor, dfdt, d2fdt2, &
!$OMP   b_ns, b_ns_t, b_ns_p, b_ns_tt, b_ns_tp, b_ns_pp,   &
!$OMP   l, fllp1)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(recfac(0:1,0:nlmax), dalm(0:1,0:nlmax), lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac, dalm & lam_fact')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac_der(nlmax, m, lam_fact)
          f2m = 2.0_dp * m
          fm2 = real(m*m, kind=dp)
          fm  = real(m,   kind=dp)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(alm(1,ll,m))
          enddo

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ! determine lam_mm
                call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                cth_ring = cth(ithl)
                one_on_s2 = one_on_s(ithl)**2
                cotanth   = cth_ring * one_on_s(ithl)
                !           ---------- l = m ----------
                par_lm = 1  ! = (-1)^(l+m)

                lam_lm = corfac * lam_mm
                if (m >= l_min) then
                   !-----------------------
                   ! f = Y_lm * a_lm
                   factor(0:1) = lam_lm * dalm(0:1,m)
                   b_ns( 1:2) = factor(0:1) ! even
                   b_ns(-1:0) = 0.0_dp ! odd
                   !------------------------- 1st derivatives
                   ! df/dtheta = (l/tan(theta)*Y_lm - fact/sin(theta)*Y_l-1m)*a_lm
                   dfdt(0:1)   =   (m * cotanth) * factor(0:1)
                   b_ns_t( 1:2) = 0.0_dp
                   b_ns_t(-1:0) = dfdt ! different theta-parity
                   ! df/dphi = i * m * Y_lm * a_lm
                   b_ns_p( 1:2) = m  * (/ -factor(1), factor(0) /)
                   b_ns_p(-1:0) = 0.0_dp
                   !------------------------- 2nd derivatives
                   if (do_d2) then
                      ! d^2f/dtheta^2    = [-(l(l+1) - m^2/sin^2(theta)) Y_lm - cotan(theta) dY_lm/dtheta] * a_lm
                      b_ns_tt( 1:2) = (fm2 * one_on_s2 - fm2 - fm) * factor - cotanth * dfdt
                      b_ns_tt(-1:0) = 0.0_dp
                      ! d^2f/dtheta/dphi = i * m * df/dtheta
                      b_ns_tp( 1:2) = 0.0_dp
                      b_ns_tp(-1:0) = m * (/ -dfdt(1), dfdt(0) /)! different theta-parity
                      ! d^2f/dphi^2      = -m^2 * Y_lm * a_lm
                      b_ns_pp( 1:2) = - fm2 * factor
                      b_ns_pp(-1:0) = 0.0_dp
                   endif
                   !-------------------------
                else
                   b_ns   = 0.0_dp
                   b_ns_t = 0.0_dp ; b_ns_p = 0.0_dp
                   b_ns_tt = 0.0_dp ; b_ns_tp = 0.0_dp ; b_ns_pp = 0.0_dp
                endif

                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                lam_2 = cth_ring * lam_1 * recfac(0,m)
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = (-1)^(l+m)
                   fllp1 = real(l*l + l, kind=dp) ! l(l+1)
                   lam_lm1m = lam_lm * lam_fact(l) ! actual lambda_l-1,m
                   lam_lm = lam_2 * corfac * lam_mm
                   if (l >= l_min) then
                      !--------------------------
                      ! f = Y_lm * a_lm
                      factor(0:1) = lam_lm * dalm(0:1,l)
                      b_ns(par_lm:par_lm+1) = b_ns(par_lm:par_lm+1) + factor(0:1)
                      !-------------------------- 1st derivatives
                      ! df/dtheta = (l/tan(theta)*Y_lm - fact/sin(theta)*Y_l-1m)*a_lm
                      dfdt(0:1)   =   (l * cotanth) * factor(0:1) &
                           &        - (one_on_s(ithl) * lam_lm1m) * dalm(0:1,l)
                      b_ns_t(-par_lm:1-par_lm) = b_ns_t(-par_lm:1-par_lm) + dfdt(0:1)
                      ! df/dphi = i * m * Y_lm * a_lm
                      b_ns_p(par_lm  ) = b_ns_p(par_lm  ) - m * factor(1)
                      b_ns_p(par_lm+1) = b_ns_p(par_lm+1) + m * factor(0)
                      !-------------------------- 2nd derivatives
                      if (do_d2) then
                         ! d^2f/dtheta^2    = [-(l(l+1) - m^2/sin^2(theta)) Y_lm - cotan(theta) dY_lm/dtheta] * a_lm
                         d2fdt2(0:1) = (fm2*one_on_s2 - fllp1) * factor(0:1) - cotanth * dfdt(0:1)
                         b_ns_tt(par_lm:par_lm+1) = b_ns_tt(par_lm:par_lm+1) + d2fdt2(0:1)
                         ! d^2f/dtheta/dphi = i * m * df/dtheta
                         b_ns_tp( -par_lm) = b_ns_tp( -par_lm) - m * dfdt(1)
                         b_ns_tp(1-par_lm) = b_ns_tp(1-par_lm) + m * dfdt(0)
                         ! d^2f/dphi^2      = -m^2 * Y_lm * a_lm
                         b_ns_pp(par_lm:par_lm+1) = b_ns_pp(par_lm:par_lm+1) - fm2 * factor(0:1)
                      endif
                      !-------------------------
                   endif
                   lam_0 = lam_1 * recfac(1,l-1)
                   lam_1 = lam_2
                   lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                   if (abs(lam_2) > OVFLOW) then
                      lam_1 = lam_1*UNFLOW
                      lam_2 = lam_2*UNFLOW
                      scalel= scalel + 1
                      corfac = rescale_tab(max(scalel+scalem,RSMIN))
                   elseif (abs(lam_2) < UNFLOW) then
                      lam_1 = lam_1*OVFLOW
                      lam_2 = lam_2*OVFLOW
                      scalel= scalel - 1
                      corfac = rescale_tab(max(scalel+scalem,RSMIN))
                   endif
                enddo ! loop on l

                one_on_s1 = one_on_s(ithl)
                b_north(m,ithl)   = cmplx(b_ns(1) + b_ns(-1), b_ns(2) + b_ns(0), kind=DP)
                b_south(m,ithl)   = cmplx(b_ns(1) - b_ns(-1), b_ns(2) - b_ns(0), kind=DP)
                b_north_t(m,ithl) = cmplx(b_ns_t(1) + b_ns_t(-1), b_ns_t(2) + b_ns_t(0), kind=DP)
                b_south_t(m,ithl) = cmplx(b_ns_t(1) - b_ns_t(-1), b_ns_t(2) - b_ns_t(0), kind=DP)
                b_north_p(m,ithl) = cmplx(b_ns_p(1) + b_ns_p(-1), b_ns_p(2) + b_ns_p(0), kind=DP)*one_on_s1
                b_south_p(m,ithl) = cmplx(b_ns_p(1) - b_ns_p(-1), b_ns_p(2) - b_ns_p(0), kind=DP)*one_on_s1
                if (do_d2) then
                   b_north_tt(m,ithl) = cmplx(b_ns_tt(1) + b_ns_tt(-1), b_ns_tt(2) + b_ns_tt(0), kind=DP)
                   b_south_tt(m,ithl) = cmplx(b_ns_tt(1) - b_ns_tt(-1), b_ns_tt(2) - b_ns_tt(0), kind=DP)
                   b_north_tp(m,ithl) = cmplx(b_ns_tp(1) + b_ns_tp(-1), b_ns_tp(2) + b_ns_tp(0), kind=DP)*one_on_s1
                   b_south_tp(m,ithl) = cmplx(b_ns_tp(1) - b_ns_tp(-1), b_ns_tp(2) - b_ns_tp(0), kind=DP)*one_on_s1
                   b_north_pp(m,ithl) = cmplx(b_ns_pp(1) + b_ns_pp(-1), b_ns_pp(2) + b_ns_pp(0), kind=DP)*one_on_s2
                   b_south_pp(m,ithl) = cmplx(b_ns_pp(1) - b_ns_pp(-1), b_ns_pp(2) - b_ns_pp(0), kind=DP)*one_on_s2
                endif
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) deallocate (recfac,dalm, lam_fact)
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_north, b_south, b_north_t, b_south_t, &
!$OMP      b_north_p,  b_south_p,  b_north_tt, b_south_tt, &
!$OMP      b_north_tp, b_south_tp, b_north_pp, b_south_pp, &
!$OMP      nph, startpix, kphi0, map, der1, der2, do_d2) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1), stat = status)
          call assert_alloc(status,code,'ring')
       endif

!$OMP do schedule(dynamic,1)

       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk

       !        ---------------------------------------------------------------
       !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta)
       !        ---------------------------------------------------------------
          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),   nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          map(istart_north:istart_north+nphl-1) = ring(0:nphl-1)

          call ring_synthesis(nsmax,nlmax,nmmax,b_north_t(0,ithl), nphl,ring,kphi0(ithl))
          der1(istart_north:istart_north+nphl-1,1) = ring(0:nphl-1)
          call ring_synthesis(nsmax,nlmax,nmmax,b_north_p(0,ithl), nphl,ring,kphi0(ithl))
          der1(istart_north:istart_north+nphl-1,2) = ring(0:nphl-1)
          if (do_d2) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_north_tt(0,ithl),nphl,ring,kphi0(ithl))
             der2(istart_north:istart_north+nphl-1,1) = ring(0:nphl-1)
             call ring_synthesis(nsmax,nlmax,nmmax,b_north_tp(0,ithl),nphl,ring,kphi0(ithl))
             der2(istart_north:istart_north+nphl-1,2) = ring(0:nphl-1)
             call ring_synthesis(nsmax,nlmax,nmmax,b_north_pp(0,ithl),nphl,ring,kphi0(ithl))
             der2(istart_north:istart_north+nphl-1,3) = ring(0:nphl-1)
          endif
          
          if (ith < nrings) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),  nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             map(istart_south:istart_south+nphl-1) = ring(0:nphl-1)
             
             call ring_synthesis(nsmax,nlmax,nmmax,b_south_t(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             der1(istart_south:istart_south+nphl-1,1) = ring(0:nphl-1)
             call ring_synthesis(nsmax,nlmax,nmmax,b_south_p(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             der1(istart_south:istart_south+nphl-1,2) = ring(0:nphl-1)
             if (do_d2) then
                call ring_synthesis(nsmax,nlmax,nmmax,b_south_tt(0,ithl),nphl,ring,kphi0(ithl))
                der2(istart_south:istart_south+nphl-1,1) = ring(0:nphl-1)
                call ring_synthesis(nsmax,nlmax,nmmax,b_south_tp(0,ithl),nphl,ring,kphi0(ithl))
                der2(istart_south:istart_south+nphl-1,2) = ring(0:nphl-1)
                call ring_synthesis(nsmax,nlmax,nmmax,b_south_pp(0,ithl),nphl,ring,kphi0(ithl))
                der2(istart_south:istart_south+nphl-1,3) = ring(0:nphl-1)
             endif
          endif

       enddo    ! loop on ithl
!$OMP end do
       if (do_openmp()) deallocate(ring)
!$OMP end parallel

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate(recfac, dalm, lam_fact, ring)
    endif

    deallocate(mfac)
    deallocate(b_north, b_north_t, b_north_p)
    deallocate(b_south, b_south_t, b_south_p)
    if (do_d2) then
       deallocate(b_north_tt, b_north_tp, b_north_pp)
       deallocate(b_south_tt, b_south_tp, b_south_pp)
    endif

    return
  end subroutine alm2map_sc_der_KLOAD

  !=======================================================================
  subroutine alm2map_pol_der_KLOAD(nsmax, nlmax, nmmax, alm_TGC, map_TQU, der1, der2)
    !=======================================================================
    !     computes a map and its 1st and 2nd derivative 
    !       from its alm for the HEALPIX pixelisation
    !      for the Temperature and polarisation fields
    !
    ! der1 = dX/d(theta), dX/d(phi)
    ! der2 = d^2X/d(theta^2), d^2X/d(theta)/d(phi), d^2X/d(phi^2)
    ! where X = (T,Q,U)
    !
    ! derivatives of Ylm are obtained from Qu.Th. of Ang.Mom. Varshalovich et al
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(KALMC), INTENT(IN),  dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1,1:3)     :: map_TQU
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1,1:6) :: der1
    real(KMAP),   INTENT(OUT), dimension(0:12*nsmax**2-1,1:9), optional :: der2

    integer(I4B) :: l, m, ith ! alm related
    integer(I4B) :: istart_south, istart_north   ! map related
    integer(I4B) :: nrings, npix, nphmx

    integer(I4B)                              :: par_lm, k, k0, k1, di
    real(DP)             :: lam_mm, cth_ring
    real(DP)                  :: f2m, fm2, fm, fllp1, lam_lm1m, f2, f3, fl
    real(dp)                  :: cotanth, one_on_s1, one_on_s2
    real(dp)                  :: a0, xp, at, aq, derW, derX, derY
    real(dp)                  :: b0t, b0p, bx, der2W, der2X, der2Y
    real(DP), dimension(-3:8) :: b_ns
    real(DP), dimension(-3:8) :: b_ns_p, b_ns_t
    real(DP), dimension(-3:8) :: b_ns_pp, b_ns_tt, b_ns_tp
    real(DP), dimension( 0:1) :: factor

    integer(i4b)                            :: ll, l_min
    complex(DPC), dimension(:,:,:), ALLOCATABLE :: b_d1, b_d2
    complex(DPC), dimension(:), allocatable :: bline
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    real(DP), dimension(:,:), ALLOCATABLE :: recfac, dalm, lam_lm
    real(DP), dimension(:),   ALLOCATABLE :: lam_fact, mfac, lam_fact_der, normal_l
    real(DP), dimension(:),   ALLOCATABLE :: ring

    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph,kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth, one_on_s

    character(len=*), parameter :: code = 'ALM2MAP_DER'
    logical(lgt) :: do_d2
    integer(I4B) :: status
    !=======================================================================
    do_d2 = (present(der2)) ! do 2nd derivative

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(b_d1(0:17,0:nmmax,0:chunksize-1), stat = status)
    if (do_d2) then
       allocate(b_d2(0:17,0:nmmax,0:chunksize-1),  stat = status)
    endif
    call assert_alloc(status,code,'b_d1 & b_d2')

    if (.not. do_openmp()) then
       allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax), lam_fact_der(0:nlmax), stat = status)
       call assert_alloc(status,code,'recfac, dalm & lam_fact')
       allocate(ring(0:nphmx-1), bline(0:nmmax),stat = status)
       call assert_alloc(status,code,'ring & bline')
       allocate(lam_lm(1:3,0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_lm')
    endif

    !     ------------ initiate arrays ----------------
    call gen_mfac(nmmax,mfac)

    call init_rescale()

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    map_TQU = 0.0 ! set the whole map to zero
    der1 = 0.0
    if (do_d2) der2 = 0.0

    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s(ithl)  = 1.0_dp / sth(ithl)
       enddo
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m)
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       b_d1   = 0_dpc   ! pad with zeros
       if (do_d2) b_d2 = 0_dpc

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, &
!$OMP    rescale_tab, normal_l, &
!$OMP    cth, sth, mfac, alm_TGC, one_on_s,  &
!$OMP      b_d1, b_d2, do_d2 ) &
!$OMP private(recfac, dalm, lam_fact, lam_fact_der, status, &
!$OMP   m, ll, fm, f2m, fm2, ithl, l_min, k, k0, k1, &
!$OMP   par_lm, lam_lm, lam_lm1m, &
!$OMP   cth_ring, one_on_s1, one_on_s2, cotanth, factor, &
!$OMP   b_ns, b_ns_t, b_ns_p, b_ns_tt, b_ns_tp, b_ns_pp,   &
!$OMP   l, fllp1, fl, a0, xp, at, aq, derW, derX, derY, f2, f3, b0t, b0p, bx, der2W, der2X, der2Y)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(recfac(0:1,0:nlmax), dalm(0:5,0:nlmax), lam_fact(0:nlmax), lam_fact_der(0:nlmax), &
               & stat = status)
          call assert_alloc(status,code,'recfac, dalm & lam_fact')
          allocate(lam_lm(1:3,0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_lm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac_der(nlmax, m, lam_fact_der)
          call gen_lamfac    (nlmax, m, lam_fact)
          f2m = 2.0_dp * m
          fm2 = real(m*m, kind=dp)
          fm  = real(m,   kind=dp)

          ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             dalm(0,ll) =  real(alm_TGC(1,ll,m),kind=dp)
             dalm(1,ll) = aimag(alm_TGC(1,ll,m))
             dalm(2,ll) =  real(alm_TGC(2,ll,m),kind=dp)*normal_l(ll) ! G, real
             dalm(3,ll) = aimag(alm_TGC(2,ll,m))        *normal_l(ll)
             dalm(4,ll) =  real(alm_TGC(3,ll,m),kind=dp)*normal_l(ll) ! C, real
             dalm(5,ll) = aimag(alm_TGC(3,ll,m))        *normal_l(ll)
          enddo

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small

                ! compute lam_lm(p,theta) for all l>=m
                call do_lam_lm_pol(nlmax, m, cth(ithl), sth(ithl), mfac(m), recfac, lam_fact, lam_lm)

                cth_ring = cth(ithl)
                one_on_s1 = one_on_s(ithl)
                one_on_s2 = one_on_s(ithl)**2
                cotanth   = cth_ring * one_on_s(ithl)

                b_ns   = 0.0_dp
                b_ns_t = 0.0_dp  ; b_ns_p  = 0.0_dp
                b_ns_tt = 0.0_dp ; b_ns_tp = 0.0_dp ; b_ns_pp = 0.0_dp
                do l = l_min, nlmax
                   fl    = real(l, kind=dp) ! l
                   fllp1 = real(l*l + l, kind=dp) ! l(l+1)
                   par_lm = 3 ! = (-1)^(l+m)
                   if (mod(l+m,2) == 1) par_lm = -par_lm

                   !--------------------------
                   ! f = Y_lm * a_lm
                   factor(0:1) = lam_lm(1,l) * dalm(0:1,l)
                   b_ns(par_lm:par_lm+1)   = b_ns(par_lm:  par_lm+1) + factor(0:1) ! T even
                   b_ns(par_lm+2:par_lm+5) = b_ns(par_lm+2:par_lm+5) - lam_lm(2,l) * dalm(2:5,l) ! Q, U  even
                   b_ns(2-par_lm)          = b_ns(2-par_lm)          + lam_lm(3,l) * dalm(5,l) ! Q odd
                   b_ns(3-par_lm)          = b_ns(3-par_lm)          - lam_lm(3,l) * dalm(4,l)
                   b_ns(4-par_lm)          = b_ns(4-par_lm)          - lam_lm(3,l) * dalm(3,l) ! U odd
                   b_ns(5-par_lm)          = b_ns(5-par_lm)          + lam_lm(3,l) * dalm(2,l)
                   !-------------------------- 1st derivatives
                   if (l > 0) then
                      ! df/dphi = i * m * Y_lm * a_lm
                      f2 = m * lam_lm(2,l)
                      f3 = m * lam_lm(3,l)
                      b_ns_p(par_lm  ) = b_ns_p(par_lm  ) - m * factor(1)
                      b_ns_p(par_lm+1) = b_ns_p(par_lm+1) + m * factor(0)
                      b_ns_p(par_lm+2) = b_ns_p(par_lm+2) + f2 * dalm(3,l)
                      b_ns_p(par_lm+3) = b_ns_p(par_lm+3) - f2 * dalm(2,l)
                      b_ns_p(par_lm+4) = b_ns_p(par_lm+4) + f2 * dalm(5,l)
                      b_ns_p(par_lm+5) = b_ns_p(par_lm+5) - f2 * dalm(4,l)
                      b_ns_p(2-par_lm) = b_ns_p(2-par_lm) + f3 * dalm(4,l) ! Q odd
                      b_ns_p(3-par_lm) = b_ns_p(3-par_lm) + f3 * dalm(5,l)
                      b_ns_p(4-par_lm) = b_ns_p(4-par_lm) - f3 * dalm(2,l) ! U odd
                      b_ns_p(5-par_lm) = b_ns_p(5-par_lm) - f3 * dalm(3,l)

                      ! dY_lm/dtheta = (l/tan(theta)*Y_lm                         -fact/sin(theta)*Y_l-1m)
                      ! dW_lm/dtheta = (l/tan(theta)*W_lm - S*m/l/sin(theta)*X_lm -fact/sin(theta)*sqrt(1-S^2/l^2)*W_l-1m
                      ! dX_lm/dtheta = (l/tan(theta)*X_lm - S*m/l/sin(theta)*W_lm -fact/sin(theta)*sqrt(1-S^2/l^2)*X_l-1m
                      a0 = fl * cotanth          ! l/tan(theta)
                      at = lam_fact_der(l) * one_on_s1 ! sqrt((2l+1)/(2l-1)*(l^2-m^2))/sin(theta)
                      derY   = a0 * lam_lm(1,l) - at * lam_lm(1,l-1)                   
                      b_ns_t( -par_lm:1-par_lm) = b_ns_t( -par_lm:1-par_lm) + derY * dalm(0:1,l) ! T odd
                   endif
                   if (l > 1) then
                      xp = (2*m)*one_on_s1/fl  ! spin m / (l sin(theta))
                      aq = at * sqrt(1.0_dp - 4.0_dp/(fl*fl)) ! at * sqrt(l^2-spin^2)/l
                      derW   = a0 * lam_lm(2,l) - aq * lam_lm(2,l-1) + xp * lam_lm(3,l)
                      derX   = a0 * lam_lm(3,l) - aq * lam_lm(3,l-1) + xp * lam_lm(2,l)
                      b_ns_t(2-par_lm:5-par_lm) = b_ns_t(2-par_lm:5-par_lm) - derW * dalm(2:5,l) ! Q, U  odd
                      b_ns_t(2+par_lm)          = b_ns_t(2+par_lm)          + derX * dalm(5,l) ! Q even
                      b_ns_t(3+par_lm)          = b_ns_t(3+par_lm)          - derX * dalm(4,l)
                      b_ns_t(4+par_lm)          = b_ns_t(4+par_lm)          - derX * dalm(3,l) ! U even
                      b_ns_t(5+par_lm)          = b_ns_t(5+par_lm)          + derX * dalm(2,l)
                   endif
                   !-------------------------- 2nd derivatives
                   if (do_d2 .and. l > 0) then
                      ! d^2 Y/dtheta^2    = -[l(l+1) - (m^2)/sin^2(theta)] Y                                     - cotan(theta) dY/dtheta
                      ! d^2 W/dtheta^2    = -[l(l+1) - (m^2+s^2)/sin^2(theta)] W - 2ms cos(theta)/sin^2(theta) X - cotan(theta) dW/dtheta
                      ! d^2 X/dtheta^2    = -[l(l+1) - (m^2+s^2)/sin^2(theta)] X - 2ms cos(theta)/sin^2(theta) W - cotan(theta) dX/dtheta
                      ! s = -2
                      b0t = fm2*one_on_s2 - fllp1
                      der2Y = b0t * lam_lm(1,l)                    - cotanth * derY
                      b_ns_tt(par_lm:par_lm+1)   = b_ns_tt(par_lm:par_lm+1)   + der2Y * dalm(0:1,l)
                      b_ns_tp( -par_lm)          = b_ns_tp( -par_lm)          - (fm * derY) * dalm(1,l)
                      b_ns_tp(1-par_lm)          = b_ns_tp(1-par_lm)          + (fm * derY) * dalm(0,l)
                      b_ns_pp(par_lm:par_lm+1)   = b_ns_pp(par_lm:par_lm+1)   - fm2 * factor(0:1)

                      if (l > 1) then
                         b0p = b0t + 4._dp*one_on_s2
                         bx  = 4*m * cth_ring * one_on_s2
                         der2W = b0p * lam_lm(2,l) + bx * lam_lm(3,l) - cotanth * derW
                         der2X = b0p * lam_lm(3,l) + bx * lam_lm(2,l) - cotanth * derX
                         b_ns_tt(par_lm+2:par_lm+5) = b_ns_tt(par_lm+2:par_lm+5) - der2W * dalm(2:5,l)
                         b_ns_tt(2-par_lm)          = b_ns_tt(2-par_lm)          + der2X * dalm(5,l) ! Q odd
                         b_ns_tt(3-par_lm)          = b_ns_tt(3-par_lm)          - der2X * dalm(4,l)
                         b_ns_tt(4-par_lm)          = b_ns_tt(4-par_lm)          - der2X * dalm(3,l) ! U odd
                         b_ns_tt(5-par_lm)          = b_ns_tt(5-par_lm)          + der2X * dalm(2,l)
                         ! d^2f/dtheta/dphi = i * m * df/dtheta
                         b_ns_tp(par_lm+2:par_lm+3) = b_ns_tp(par_lm+2:par_lm+3) + (fm * derX) * dalm(4:5,l)
                         b_ns_tp(par_lm+4:par_lm+5) = b_ns_tp(par_lm+4:par_lm+5) - (fm * derX) * dalm(2:3,l)
                         b_ns_tp(2-par_lm)          = b_ns_tp(2-par_lm)          + (fm * derW) * dalm(3,l)
                         b_ns_tp(3-par_lm)          = b_ns_tp(3-par_lm)          - (fm * derW) * dalm(2,l)
                         b_ns_tp(4-par_lm)          = b_ns_tp(4-par_lm)          + (fm * derW) * dalm(5,l)
                         b_ns_tp(5-par_lm)          = b_ns_tp(5-par_lm)          - (fm * derW) * dalm(4,l)
                         ! d^2f/dphi^2      = -m^2 * Y_lm * a_lm
                         b_ns_pp(par_lm+2:par_lm+5) = b_ns_pp(par_lm+2:par_lm+5) +(fm2 * lam_lm(2,l))* dalm(2:5,l) ! Q, U  even
                         b_ns_pp(2-par_lm)          = b_ns_pp(2-par_lm)          -(fm2 * lam_lm(3,l))* dalm(5,l) ! Q odd
                         b_ns_pp(3-par_lm)          = b_ns_pp(3-par_lm)          +(fm2 * lam_lm(3,l))* dalm(4,l)
                         b_ns_pp(4-par_lm)          = b_ns_pp(4-par_lm)          +(fm2 * lam_lm(3,l))* dalm(3,l) ! U odd
                         b_ns_pp(5-par_lm)          = b_ns_pp(5-par_lm)          -(fm2 * lam_lm(3,l))* dalm(2,l)
                      endif
                   endif
                   !-------------------------
                enddo ! loop on l
                do k=0,2 ! loop on T,Q,U
                   k0 = 2*k
                   k1 = k0+1
                   ! fields
                   b_d1(0+k,m,ithl)   = cmplx(b_ns(k0+3) + b_ns(k0-3), &
                        &                     b_ns(k1+3) + b_ns(k1-3), kind=DP) ! north=Even+Odd
                   b_d1(3+k,m,ithl)   = cmplx(b_ns(k0+3) - b_ns(k0-3), &
                        &                     b_ns(k1+3) - b_ns(k1-3), kind=DP) ! south=Even-Odd
                   ! dfield/dtheta
                   b_d1(6+k,m,ithl)   = cmplx(b_ns_t(k0+3) + b_ns_t(k0-3), &
                        &                     b_ns_t(k1+3) + b_ns_t(k1-3), kind=DP) ! north=Even+Odd
                   b_d1(9+k,m,ithl)   = cmplx(b_ns_t(k0+3) - b_ns_t(k0-3), &
                        &                     b_ns_t(k1+3) - b_ns_t(k1-3), kind=DP) ! south=Even-Odd
                   ! dfield/dphi/sin(theta)
                   b_d1(12+k,m,ithl)  = cmplx(b_ns_p(k0+3) + b_ns_p(k0-3), &
                        &                     b_ns_p(k1+3) + b_ns_p(k1-3), kind=DP)*one_on_s1
                   b_d1(15+k,m,ithl)  = cmplx(b_ns_p(k0+3) - b_ns_p(k0-3), &
                        &                     b_ns_p(k1+3) - b_ns_p(k1-3), kind=DP)*one_on_s1
                enddo
                if (do_d2) then
                   do k=0,2 ! loop on T,Q,U
                      k0 = 2*k
                      k1 = k0+1
                      ! dfield/dtheta^2
                      b_d2(0+k,m,ithl)   = cmplx(b_ns_tt(k0+3) + b_ns_tt(k0-3), &
                           &                     b_ns_tt(k1+3) + b_ns_tt(k1-3), kind=DP)
                      b_d2(3+k,m,ithl)   = cmplx(b_ns_tt(k0+3) - b_ns_tt(k0-3), &
                           &                     b_ns_tt(k1+3) - b_ns_tt(k1-3), kind=DP)
                      ! dfield/dtheta/dphi/sin(theta)
                      b_d2(6+k,m,ithl)   = cmplx(b_ns_tp(k0+3) + b_ns_tp(k0-3), &
                           &                     b_ns_tp(k1+3) + b_ns_tp(k1-3), kind=DP)*one_on_s1
                      b_d2(9+k,m,ithl)   = cmplx(b_ns_tp(k0+3) - b_ns_tp(k0-3), &
                           &                     b_ns_tp(k1+3) - b_ns_tp(k1-3), kind=DP)*one_on_s1
                      ! dfield/dphi&2/sin(theta)^2
                      b_d2(12+k,m,ithl)  = cmplx(b_ns_pp(k0+3) + b_ns_pp(k0-3), &
                           &                     b_ns_pp(k1+3) + b_ns_pp(k1-3), kind=DP)*one_on_s2
                      b_d2(15+k,m,ithl)  = cmplx(b_ns_pp(k0+3) - b_ns_pp(k0-3), &
                           &                     b_ns_pp(k1+3) - b_ns_pp(k1-3), kind=DP)*one_on_s2
                   enddo
                endif
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) deallocate (recfac,dalm, lam_fact, lam_fact_der)
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, b_d1, b_d2, &
!$OMP      nph, startpix, kphi0, map_TQU, der1, der2, do_d2) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, bline, status, k0, di)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1), stat = status)
          call assert_alloc(status,code,'ring')
          allocate(bline(0:nmmax),stat = status)
          call assert_alloc(status,code,'bline')
       endif

!$OMP do schedule(dynamic,1)

       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk

       !        ---------------------------------------------------------------
       !        sum_m  b(m,theta)*exp(i*m*phi)   -> f(phi,theta)
       !        ---------------------------------------------------------------
          do k0=0,2
             bline = b_d1(k0, 0:nmmax, ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             map_TQU(istart_north:istart_north+nphl-1,k0+1) = ring(0:nphl-1)
          enddo
          do k0=0,2
             bline = b_d1(6+k0, 0:nmmax, ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             der1(istart_north:istart_north+nphl-1,k0+1) = ring(0:nphl-1)
          enddo
          do k0=0,2
             bline = b_d1(12+k0, 0:nmmax, ithl)
             call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! north hemisph. + equator
             der1(istart_north:istart_north+nphl-1,k0+4) = ring(0:nphl-1)
          enddo

          if (ith < nrings) then
             do k0=0,2
                bline = b_d1(3+k0, 0:nmmax, ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! south hemisph. w/o equat
                map_TQU(istart_south:istart_south+nphl-1,k0+1) = ring(0:nphl-1)
             enddo
             do k0=0,2
                bline = b_d1(9+k0, 0:nmmax, ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! south hemisph. w/o equat
                der1(istart_south:istart_south+nphl-1,k0+1) = ring(0:nphl-1)
             enddo
             do k0=0,2
                bline = b_d1(15+k0, 0:nmmax, ithl)
                call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! south hemisph. w/o equat
                der1(istart_south:istart_south+nphl-1,k0+4) = ring(0:nphl-1)
             enddo
          endif
          if (do_d2) then
             do di = 0,2
                do k0=0,2
                   bline = b_d2(k0+6*di, 0:nmmax, ithl)
                   call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! north hemisph. + equator
                   der2(istart_north:istart_north+nphl-1,k0+3*di+1) = ring(0:nphl-1)
                enddo
             enddo
             if (ith < nrings) then
                do di = 0,2
                   do k0=0,2
                      bline = b_d2(k0+6*di+3, 0:nmmax, ithl)
                      call ring_synthesis(nsmax,nlmax,nmmax,bline,   nphl,ring,kphi0(ithl))   ! south hemisph. w/o equat
                      der2(istart_south:istart_south+nphl-1,k0+3*di+1) = ring(0:nphl-1)
                   enddo
                enddo
             endif
          endif
       enddo    ! loop on ithl
!$OMP end do
       if (do_openmp()) then
          deallocate(ring, bline)
       endif
!$OMP end parallel

    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate(recfac, dalm, lam_fact, lam_fact_der, ring, bline)
    endif

    deallocate(mfac)
    deallocate(b_d1)
    if (do_d2) then
       deallocate(b_d2)
    endif

    return
  end subroutine alm2map_pol_der_KLOAD

  !**************************************************************************
  !
  !             MAP2ALM
  ! 
  !**************************************************************************
  !=======================================================================
  !  map2alm(nsmax, nlmax, nmmax, map, alm [, zbound, w8ring ,plm])
  !
  !     computes the a(l,m) from a map (temperature only or polarized) 
  !           for the HEALPIX pixelisation
  !
  !     For the Temperature field
  !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
  !            = int dtheta lambda_lm(theta)
  !                  * int dphi T(theta,phi) e^(-i*m*phi)
  !     For polarizaton
  !     alm_G = \int ( - Q Wlm - i U Xlm )
  !     alm_C = \int ( - U Wlm + i Q Xlm )
  !
  !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
  !
  !     * the recurrence of Ylm is the standard one (cf Num Rec)
  !     * the integral over phi is done by FFT
  !
  !     zbounds: cut apply (in cos_theta)
  !     w8ring: ring dependent weigthing scheme to improve quadrature
  !     plm: precomputed Ylm
  !=======================================================================

  ! interface with legacy code
  !=======================================================================
  subroutine map2alm_old_sc_KLOAD(nsmax, nlmax, nmmax, map, alm, cos_theta_cut, w8ring, plm)
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1) :: map
    complex(KALMC), intent(OUT), dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(DP),     intent(IN)                          :: cos_theta_cut
    real(DP),     intent(IN),  dimension(1:2*nsmax,1) :: w8ring
    real(DP),     intent(IN),  dimension(0:), optional          :: plm
    !
    real(DP),     dimension(1:2) :: zbounds
    
    call warning_oldbounds(cos_theta_cut, zbounds)
    if (present(plm)) then
       call map2alm(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring, plm)
    else
       call map2alm(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring)
    endif

    return
  end subroutine map2alm_old_sc_KLOAD
  !=======================================================================
  subroutine map2alm_old_pol_KLOAD(nsmax, nlmax, nmmax, map_TQU, alm_TGC, cos_theta_cut, w8ring, plm)
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1,1:3) :: map_TQU
    complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(DP),     intent(IN)                          :: cos_theta_cut
    real(DP),     intent(IN),  dimension(1:2*nsmax,1:3) :: w8ring
    real(DP),     intent(IN),  dimension(0:), optional          :: plm
    !
    real(DP),     dimension(1:2) :: zbounds
    
    call warning_oldbounds(cos_theta_cut, zbounds)
    if (present(plm)) then
       call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring, plm)
    else
       call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring)
    endif
    return
  end subroutine map2alm_old_pol_KLOAD

  !=======================================================================
  subroutine map2alm_old_pol2_KLOAD(nsmax, nlmax, nmmax, map_TQU, alm_TGC, cos_theta_cut, w8ring, plm)
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1,1:3) :: map_TQU
    complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(DP),     intent(IN)                          :: cos_theta_cut
    real(DP),     intent(IN),  dimension(1:2*nsmax,1:3) :: w8ring
    real(DP),     intent(IN),  dimension(0:,1:)         :: plm
    !
    real(DP),     dimension(1:2) :: zbounds
    
    call warning_oldbounds(cos_theta_cut, zbounds)
    call map2alm(nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring, plm)

    return
  end subroutine map2alm_old_pol2_KLOAD

!   !=======================================================================
!   subroutine map2alm_sc_test_KLOAD(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring)
!     !=======================================================================
!     !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
!     !        all from scratch
!     !=======================================================================
!     integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
!     real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1) :: map
!     complex(KALMC), intent(OUT), dimension(1:1,0:nlmax,0:nmmax) :: alm
!     real(DP),     intent(IN),  dimension(1:2),         optional :: zbounds
!     real(DP),     intent(IN),  dimension(1:2*nsmax,1), optional :: w8ring

!     real(DP), dimension(1:2)         :: zbounds_in
!     real(DP), dimension(1:2*nsmax,1) :: w8ring_in
!     integer(I4B) :: l, m, ith                   ! alm related
!     integer(I4B) :: istart_south, istart_north  ! map related
!     integer(I4B) :: nrings, npix, nphmx
!     real(DP)     :: omega_pix

!     real(DP),     dimension(-1:2)             :: phas_sd
!     real(DP),     dimension(:,:), allocatable :: dalm, recfac
!     real(DP),     dimension(:),   allocatable :: mfac, lam_lm

!     integer(I4B)                              :: l_min, l_start
!     complex(DPC)                              :: php, phm
!     complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
!     real(DP),     dimension(:),   allocatable :: ring
!     integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
!     integer(I8B), dimension(0:SMAXCHK-1) :: startpix
!     integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
!     real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
!     logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

!     character(LEN=*), PARAMETER :: code = 'MAP2ALM'
!     integer(I4B) :: status
!     !=======================================================================

!     zbounds_in = (/-1.d0 , 1.d0/)
!     if (present(zbounds)) zbounds_in = zbounds
!     w8ring_in  = 1.d0
!     if (present(w8ring))  w8ring_in  = w8ring

!     ! Healpix definitions
!     nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!     npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!     nphmx  = 4*nsmax           ! maximum number of pixels/ring
!     omega_pix = FOURPI / real(npix, kind=DP)  ! pixel area (identical for all pixels)

!     !     --- allocates space for arrays ---
!     nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!     chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

!     allocate(mfac(0:nmmax),stat = status)
!     call assert_alloc(status,code,'mfac')

!     allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'phas_n')

!     allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'phas_s')

!     if (.not.do_openmp()) then
!        allocate(ring(0:nphmx-1),stat = status)
!        call assert_alloc(status,code,'ring')
!        allocate(recfac(0:1,0:nlmax),stat = status)
!        call assert_alloc(status,code,'recfac')
!        allocate(dalm(1:2,0:nlmax),stat = status)
!        call assert_alloc(status,code,'dalm')
!        allocate(lam_lm(0:nlmax), stat = status)
!        call assert_alloc(status,code,'lam_lm')
!     endif
!     !     ------------ initiate variables and arrays ----------------

!     call gen_mfac(nmmax,mfac)

!     call init_rescale()
!     alm = 0.0 ! set the whole alm array to zero

!     ! loop on chunks
!     do ichunk = 0, nchunks-1
!        lchk = ichunk * chunksize + 1
!        uchk = min(lchk+chunksize - 1, nrings)

!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           ! get pixel location information
!           call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!           ! find out which rings are to be analysed
!           call select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
!        enddo

!        !-----------------------------------------------------------------------
!        !  computes the integral in phi : phas_m(theta)
!        !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
!        !-----------------------------------------------------------------------
!        phas_n = 0_dpc
!        phas_s = 0_dpc

! !$OMP parallel default(none) &
! !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
! !$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, phas_n, phas_s, &
! !$OMP      keep_north, keep_south, map) &
! !$OMP  private(ithl, nphl, istart_north, istart_south, ith, ring, status)

!        if (do_openmp())then
!           allocate(ring(0:nphmx-1),stat = status)
!           call assert_alloc(status,code,'ring')
!        endif
! !$OMP do schedule(dynamic,1)
!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           nphl = nph(ithl)
!           istart_north = startpix(ithl)
!           istart_south = npix-istart_north-nphl
!           ! do Fourier Transform on rings
!           if (keep_north(ithl)) then
!              ring(0:nphl-1) = map(istart_north:istart_north+nphl-1) * w8ring_in(ith,1)
!              call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
!           endif

!           if (ith < nrings .and. keep_south(ithl)) then
!              ring(0:nphl-1) = map(istart_south:istart_south+nphl-1) * w8ring_in(ith,1)
!              call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
!           endif
!        enddo ! loop on ring
! !$OMP end do
!        if (do_openmp()) then
!           deallocate(ring)
!        endif
! !$OMP end parallel
!        !-----------------------------------------------------------------------
!        !              computes the a_lm by integrating over theta
!        !                  lambda_lm(theta) * phas_m(theta)
!        !                         for each m and l
!        !-----------------------------------------------------------------------

! !$OMP parallel default(none) &
! !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, &
! !$OMP    cth, sth, mfac, alm, phas_n, phas_s, keep_it, omega_pix) &
! !$OMP private(recfac, dalm, lam_lm, phas_sd, status, m, l, ithl, l_min, l_start, php, phm)

!        if (do_openmp()) then
!           allocate(recfac(0:1,0:nlmax),stat = status)
!           call assert_alloc(status,code,'recfac')
!           allocate(dalm(1:2,0:nlmax),stat = status)
!           call assert_alloc(status,code,'dalm')
!           allocate(lam_lm(0:nlmax), stat = status)
!           call assert_alloc(status,code,'lam_lm')
!        endif

! !$OMP do schedule(dynamic,1)
!        do m = 0, nmmax
!           ! generate recursion factors (recfac) for Ylm of degree m
!           call gen_recfac(nlmax, m, recfac)

!           ! introduce double precision vector to perform summation over ith for each l
!           dalm(1:2, m:nlmax ) = 0.0_dp

!           do ithl = 0, uchk - lchk
!              l_min = l_min_ylm(m, sth(ithl))
!              if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
!                 ! compute lam_lm(theta) for all l>=m
!                 call do_lam_lm(nlmax, m, cth(ithl), sth(ithl), mfac(m), recfac, lam_lm)

!                 php = phas_n(m,ithl) + phas_s(m,ithl) ! sum  (if (l+m) even)
!                 phm = phas_n(m,ithl) - phas_s(m,ithl) ! diff (if (l+m) odd)
!                 phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
!                 phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
!                 l_start = max(m, l_min) ! even values of (l+m)
!                 if (mod(m+l_start,2) == 1) l_start = l_start + 1
!                 do l = l_start, nlmax, 2
!                    dalm(1:2, l) = dalm(1:2, l) + lam_lm(l) * phas_sd(1:2)
!                 enddo
!                 l_start = max(m+1, l_min) ! odd values of (l+m)
!                 if (mod(m+l_start,2) == 0) l_start = l_start + 1
!                 do l = l_start, nlmax, 2
!                    dalm(1:2, l) = dalm(1:2, l) + lam_lm(l) * phas_sd(-1:0)
!                 enddo

!              endif ! test on cut sky and nlmax
!           enddo ! loop on ithl
!           do l = m, nlmax
!              alm(1, l, m) = alm(1, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP) * omega_pix
!           enddo
!        enddo ! loop on m
! !$OMP end do
!        if (do_openmp()) then
!           deallocate (recfac,dalm,lam_lm)
!        endif
! !$OMP end parallel
!     enddo ! loop on chunks

!     !     --------------------
!     !     free memory and exit
!     !     --------------------
!     deallocate(phas_n,phas_s)
!     deallocate(mfac)
!     if (.not.do_openmp()) then
!        deallocate (recfac,dalm,lam_lm,ring)
!     endif
!     RETURN
!   END subroutine map2alm_sc_test_KLOAD

  !=======================================================================
  subroutine map2alm_sc_KLOAD(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring)
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !        all from scratch
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1) :: map
    complex(KALMC), intent(OUT), dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(DP),     intent(IN),  dimension(1:2),         optional :: zbounds
    real(DP),     intent(IN),  dimension(1:2*nsmax,1), optional :: w8ring

    real(DP), dimension(1:2)         :: zbounds_in
    real(DP), dimension(1:2*nsmax,1) :: w8ring_in
    integer(I4B) :: l, m, ith, scalem, scalel   ! alm related
    integer(I4B) :: istart_south, istart_north  ! map related
    integer(I4B) :: nrings, npix, nphmx
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                      :: ovflow, unflow
    real(DP),     dimension(-1:2)             :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    real(DP),     dimension(:),   allocatable :: ring
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    zbounds_in = (/-1.d0 , 1.d0/)
    if (present(zbounds)) zbounds_in = zbounds
    w8ring_in  = 1.d0
    if (present(w8ring))  w8ring_in  = w8ring

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, kind=DP)  ! pixel area (identical for all pixels)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_s')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(1:2,0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm')
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    call init_rescale()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    alm = 0.0 ! set the whole alm array to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_n = 0_dpc
       phas_s = 0_dpc

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, phas_n, phas_s, &
!$OMP      keep_north, keep_south, map) &
!$OMP  private(ithl, nphl, istart_north, istart_south, ith, ring, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ! do Fourier Transform on rings
          if (keep_north(ithl)) then
             ring(0:nphl-1) = map(istart_north:istart_north+nphl-1) * w8ring_in(ith,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
          endif

          if (ith < nrings .and. keep_south(ithl)) then
             ring(0:nphl-1) = map(istart_south:istart_south+nphl-1) * w8ring_in(ith,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
          endif
       enddo ! loop on ring
!$OMP end do
       if (do_openmp()) then
          deallocate(ring)
       endif
!$OMP end parallel
       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, ovflow, unflow, &
!$OMP    cth, sth, mfac, alm, phas_n, phas_s, keep_it, omega_pix) &
!$OMP private(recfac, dalm, phas_sd, status, m, ithl, l_min, &
!$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
!$OMP   cth_ring, l, php, phm)

       if (do_openmp()) then
          allocate(recfac(0:1,0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac')
          allocate(dalm(1:2,0:nlmax),stat = status)
          call assert_alloc(status,code,'dalm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                !           ---------- l = m ----------
                par_lm = 1

                php = phas_n(m,ithl) + phas_s(m,ithl) ! sum  (if (l+m) even)
                phm = phas_n(m,ithl) - phas_s(m,ithl) ! diff (if (l+m) odd)
                phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
                if (m >= l_min) then
                   lam_lm = corfac * lam_mm !Actual lam_mm 
                   dalm(1:2, m) = dalm(1:2, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                endif

                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = (-1)^(l+m)
                   if (l >= l_min) then
                      lam_lm = lam_2 * corfac * lam_mm
                      dalm(1:2, l) = dalm(1:2, l) &
                           &       + lam_lm * phas_sd(par_lm:par_lm+1)
                   endif

                   lam_0 = lam_1 * recfac(1,l-1)
                   lam_1 = lam_2
                   lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                   if (abs(lam_2) > OVFLOW) then
                      lam_1 = lam_1*UNFLOW
                      lam_2 = lam_2*UNFLOW
                      scalel= scalel + 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   elseif (abs(lam_2) < UNFLOW) then
                      lam_1 = lam_1*OVFLOW
                      lam_2 = lam_2*OVFLOW
                      scalel= scalel - 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   endif

                enddo ! loop on l
             endif ! test on cut sky and nlmax
          enddo ! loop on ithl
          do l = m, nlmax
             alm(1, l, m) = alm(1, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP) * omega_pix
          enddo
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (recfac,dalm)
       endif
!$OMP end parallel
    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(phas_n,phas_s)
    deallocate(mfac)
    if (.not.do_openmp()) then
       deallocate (ring,recfac,dalm)
    endif
    RETURN
  END subroutine map2alm_sc_KLOAD

  !=======================================================================
  subroutine map2alm_sc_pre_KLOAD(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring, plm)
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !       with precomputed Ylm(theta)
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(    0:12*nsmax**2-1) :: map
    complex(KALMC), intent(OUT), dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(DP),     intent(IN),  dimension(1:2)                 :: zbounds
    real(DP),     intent(IN),  dimension(1:2*nsmax,1) :: w8ring
    real(DP),     intent(IN),  dimension(0:)          :: plm

    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    integer(I8B)                              :: n_lm, n_plm, i_mm
    integer(I4B)                              :: l, m, ith
    real(DP),     dimension(-1:2)             :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(I4B)                              :: l_min, l_start
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    real(DP),     dimension(:),   allocatable :: ring
    real(DP)                                  :: cth
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    real(DP),     dimension(0:SMAXCHK-1) :: sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), PARAMETER :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)  ! pixel area (identical for all pixels)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_n')

    allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_s')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(dalm(1:2,0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm')
    endif
    !     ------------ initiate variables and arrays ----------------
    alm = 0.0 ! set the whole alm array to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth, zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_n = 0_dpc
       phas_s = 0_dpc

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk, nph, startpix, kphi0, w8ring, phas_n, phas_s, &
!$OMP      keep_north, keep_south, map) &
!$OMP  private(ithl, nphl, istart_north, istart_south, ith, ring, status)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ! do Fourier Transform on rings
          if (keep_north(ithl)) then
             ring(0:nphl-1) = map(istart_north:istart_north+nphl-1) * w8ring(ith,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
          endif

          if (ith < nrings .and. keep_south(ithl)) then
             ring(0:nphl-1) = map(istart_south:istart_south+nphl-1) * w8ring(ith,1)
             call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
          endif
       enddo ! loop on ring
!$OMP end do
       if (do_openmp()) then
          deallocate(ring)
       endif
!$OMP end parallel
       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, ith, &
!$OMP    sth, alm, phas_n, phas_s, keep_it, plm, n_lm, omega_pix) &
!$OMP private(dalm, phas_sd, status, m, ithl, l_min, l_start, l, php, phm, i_mm)

       if (do_openmp()) then
          allocate(dalm(1:2,0:nlmax),stat = status)
          call assert_alloc(status,code,'dalm')
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax

          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax+3-m)*m)/2 ! location of Ym,m for ring ith

                php = phas_n(m,ithl) + phas_s(m,ithl) ! sum  (if (l+m) even)
                phm = phas_n(m,ithl) - phas_s(m,ithl) ! diff (if (l+m) odd)
                phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                phas_sd( 1:2) =  (/ real(php, kind=dp), aimag(php) /)
                
                !           ---------- l >= m ----------
                l_start = max(m, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(1,l) = dalm(1,l) + plm(i_mm+l-m) * phas_sd(1)
                   dalm(2,l) = dalm(2,l) + plm(i_mm+l-m) * phas_sd(2)
                enddo ! loop on l

                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(1,l) = dalm(1,l) + plm(i_mm+l-m) * phas_sd(-1)
                   dalm(2,l) = dalm(2,l) + plm(i_mm+l-m) * phas_sd(0)
                enddo ! loop on l


             endif ! test on cut sky and nlmax
          enddo ! loop on ithl
          do l = m, nlmax
             alm(1, l, m) = alm(1, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP) * omega_pix
          enddo
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (dalm)
       endif
!$OMP end parallel
    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate(dalm, ring)
    endif
    deallocate(phas_n,phas_s)
    RETURN
  END subroutine map2alm_sc_pre_KLOAD


!   !=======================================================================
!   subroutine map2alm_pol_test_KLOAD(nsmax, nlmax, nmmax, map_TQU,&
!        &                     alm_TGC, zbounds, w8ring_TQU)
!     !=======================================================================
!     !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
!     !       all from scratch
!     !=======================================================================
!     integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
!     real(KMAP),   intent(IN),  dimension(0:12*nsmax**2-1,1:3) :: map_TQU
!     complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
!     real(DP),     intent(IN),  dimension(1:2)                 :: zbounds
!     real(DP),     intent(IN), dimension(1:2*nsmax,1:3) :: w8ring_TQU

!     integer(I4B) :: l, m, ith, nrings, nphmx
!     integer(I8B) :: istart_south, istart_north, npix
!     real(DP)     :: omega_pix

!     real(DP),     dimension(-3:8)             :: phas_sd
!     real(DP),     dimension(-3:6)             :: phas_sdx
!     real(DP),     dimension(:,:), allocatable :: recfac, dalm, lam_lm
!     real(DP),     dimension(:),   allocatable :: lam_fact, mfac

!     integer(i4b)                              :: l_min, l_start
!     real(DP),     dimension(:),   allocatable :: normal_l
!     complex(DPC)                              :: php, phm
!     complex(DPC), dimension(:,:,:), allocatable :: phas_ns
!     real(DP),     dimension(:),   allocatable :: ring
!     integer(i4b)    :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, ii
!     integer(I8B), dimension(0:SMAXCHK) :: startpix
!     integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
!     real(DP),     dimension(0:SMAXCHK) :: cth, sth
!     logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it
!     real(DP), dimension(0:3) :: xone
! !    real(sp) :: time0, time1, time2, time3, time4, tt1, tt2

!     character(LEN=*), parameter :: code = 'MAP2ALM'
!     integer(I4B) :: status
!     !=======================================================================

!     ! Healpix definitions
!     nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!     npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!     nphmx  = 4*nsmax           ! maximum number of pixels/ring
!     omega_pix = FOURPI / real(npix, kind=DP)
!     xone = (/ 1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp /)

!     !     --- allocates space for arrays ---
!     nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!     chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

!     allocate(normal_l(0:nlmax),stat = status)
!     call assert_alloc(status,code,'normal_l')

!     allocate(mfac(0:nmmax),stat = status)
!     call assert_alloc(status,code,'mfac')

!     allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
!     call assert_alloc(status,code,'phas_ns')

!     if (.not.do_openmp()) then
!        allocate(ring(0:nphmx-1),stat = status)
!        call assert_alloc(status,code,'ring')
!        allocate(recfac(0:1,0:nlmax),stat = status)
!        call assert_alloc(status,code,'recfac')
!        allocate(dalm(0:5,0:nlmax), stat = status)
!        call assert_alloc(status,code,'dalm')
!        allocate(lam_fact(0:nlmax),stat = status)
!        call assert_alloc(status,code,'lam_fact')
!        allocate(lam_lm(1:3,0:nlmax),stat = status)
!        call assert_alloc(status,code,'lam_lm')
!     endif
!     !     ------------ initiate variables and arrays ----------------

!     call gen_mfac(nmmax,mfac)

!     call init_rescale()
!     alm_TGC = 0.0 ! set the whole alm array to zero

!     ! generate Polarization normalisation
!     call gen_normpol(nlmax, normal_l)

! !    tt1 = 0
!     ! loop on chunks
!     do ichunk = 0, nchunks-1
!        lchk = ichunk * chunksize + 1
!        uchk = min(lchk+chunksize - 1, nrings)

! !       call wall_clock_time(time1)
!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           ! get pixel location information
!           call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!           ! find out which rings are to be analysed
!           call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
!        enddo

!        !-----------------------------------------------------------------------
!        !  computes the integral in phi : phas_m(theta)
!        !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
!        !-----------------------------------------------------------------------
!        phas_ns = 0_dpc

! !$OMP parallel default(none) &
! !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
! !$OMP      lchk, uchk,nph, startpix, kphi0, w8ring_TQU, phas_ns, &
! !$OMP      keep_north, keep_south, map_TQU) &
! !$OMP  private(ithl, nphl, istart_north, istart_south, &
! !$OMP      ith, ring, i, ii, status)
!        if (do_openmp()) then
!           allocate(ring(0:nphmx-1),stat = status)
!           call assert_alloc(status,code,'ring')
!        endif
! !$OMP do schedule(dynamic,1)
!        ! do Fourier Transform on rings
!        do ith = lchk, uchk
!           ithl = ith - lchk !local index
!           nphl = nph(ithl)
!           istart_north = startpix(ithl)
!           istart_south = npix-istart_north-nphl
!           if (keep_north(ithl)) then
!              do i=1,3
!                 ii = 2*(i-1) ! 0,2,4
!                 ring(0:nphl-1) = map_TQU(istart_north:istart_north+nphl-1,i) * w8ring_TQU(ith,i)
!                 call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
!              enddo
!           endif

!           if (ith < nrings .and. keep_south(ithl)) then
!              do i=1,3
!                 ii = 2*i - 1 ! 1,3,5
!                 ring(0:nphl-1) = map_TQU(istart_south:istart_south+nphl-1,i) * w8ring_TQU(ith,i)
!                 call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
!              enddo
!           endif
!        enddo
! !$OMP end do
!        if (do_openmp()) then
!           deallocate (ring)
!        endif
! !$OMP end parallel

! !       call wall_clock_time(time2)
! !       tt1 = tt1 + time2 - time1

! !$OMP parallel default(none) &
! !$OMP shared(nlmax, nmmax, lchk, uchk, &
! !$OMP    rescale_tab, omega_pix, &
! !$OMP    cth, sth, keep_it, mfac, normal_l, alm_TGC, phas_ns, xone) &
! !$OMP private(recfac, dalm, lam_fact, lam_lm, phas_sd, phas_sdx, phm, php, status, &
! !$OMP   m, l, ithl, l_min, l_start )

!        if (do_openmp()) then
!           allocate(recfac(0:1,0:nlmax),stat = status)
!           call assert_alloc(status,code,'recfac')
!           allocate(dalm(0:5,0:nlmax), stat = status)
!           call assert_alloc(status,code,'dalm')
!           allocate(lam_fact(0:nlmax),stat = status)
!           call assert_alloc(status,code,'lam_fact')
!           allocate(lam_lm(1:3,0:nlmax),stat = status)
!           call assert_alloc(status,code,'lam_lm')
!        endif
! !$OMP do schedule(dynamic,1)

!        do m = 0, nmmax
!           ! generate recursion factors (recfac) for Ylm of degree m
!           call gen_recfac(nlmax, m, recfac)
!           ! generate Ylm relation factor for degree m
!           call gen_lamfac(nlmax, m, lam_fact)

!           ! introduce double precision vector to perform summation over ith for each l
!           dalm(0:5, m:nlmax) = 0.0_dp

!           do ithl = 0, uchk - lchk
!              l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
!              if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)

!                 do i=0,2 ! loop on T, Q, U
!                    phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
!                    phas_sd(-3+2*i)   =  real(phm, kind=dp) 
!                    phas_sd(-3+2*i+1) = aimag(phm)
!                    php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
!                    phas_sd( 3+2*i)   =  real(php, kind=dp) 
!                    phas_sd( 3+2*i+1) = aimag(php)
!                 enddo
!                 phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
!                 phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)
! !                 phas_sdx(-3: 0) = xone * phas_sd(8: 5:-1)
! !                 phas_sdx( 3: 6) = xone * phas_sd(2:-1:-1)

!                 ! compute lam_lm(p,theta) for all l>=m
!                 call do_lam_lm_pol(nlmax, m, cth(ithl), sth(ithl), mfac(m), recfac, lam_fact, lam_lm)

!                 ! alm_T = \int T Ylm
!                 ! alm_G = \int ( - Q Wlm - i U Xlm )
!                 ! alm_C = \int ( - U Wlm + i Q Xlm )

!                 l_start = max(m, l_min) ! even values of (l+m)
!                 if (mod(m+l_start,2) == 1) l_start = l_start + 1
!                 do l = l_start, nlmax, 2
!                    dalm(0:1, l) = dalm(0:1, l) + lam_lm(1,l) * phas_sd(3:4)
!                    dalm(2:5, l) = dalm(2:5, l) - lam_lm(2,l) * phas_sd(5:8) &
!                         &                      + lam_lm(3,l) * phas_sdx(3:6)
!                 enddo
!                 l_start = max(m+1, l_min) ! odd values of (l+m)
!                 if (mod(m+l_start,2) == 0) l_start = l_start + 1
!                 do l = l_start, nlmax, 2
!                    dalm(0:1, l) = dalm(0:1, l) + lam_lm(1,l) * phas_sd(-3:-2)
!                    dalm(2:5, l) = dalm(2:5, l) - lam_lm(2,l) * phas_sd(-1:2) &
!                         &                      + lam_lm(3,l) * phas_sdx(-3:0)
!                 enddo

!              endif ! test on cut sky
!           enddo ! loop on ithl
!           do l = m, nlmax
!              alm_TGC(1, l, m) = alm_TGC(1, l, m) + cmplx(dalm(0, l), dalm(1, l), kind=DP) * omega_pix
!              alm_TGC(2, l, m) = alm_TGC(2, l, m) + cmplx(dalm(2, l), dalm(3, l), kind=DP) * (normal_l(l) * omega_pix)
!              alm_TGC(3, l, m) = alm_TGC(3, l, m) + cmplx(dalm(4, l), dalm(5, l), kind=DP) * (normal_l(l) * omega_pix)
!           enddo
!        enddo ! loop on m
! !$OMP end do
!        if (do_openmp()) then
!           deallocate (recfac,dalm,lam_fact,lam_lm)
!        endif
! !$OMP end parallel
!     enddo    ! loop on chunks

! !    print*,'FFT',tt1
!     !     --------------------
!     !     free memory and exit
!     !     --------------------
!     if (.not.do_openmp()) then
!        deallocate (recfac,dalm,lam_fact,lam_lm,ring)
!     endif
!     deallocate(mfac, normal_l)
!     deallocate(phas_ns)
!     return
!   end subroutine map2alm_pol_test_KLOAD

  !=======================================================================
  subroutine map2alm_pol_KLOAD(nsmax, nlmax, nmmax, map_TQU,&
       &                     alm_TGC, zbounds, w8ring_TQU)
    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       all from scratch
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(0:12*nsmax**2-1,1:3) :: map_TQU
    complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(DP),     intent(IN),  dimension(1:2)                 :: zbounds
    real(DP),     intent(IN), dimension(1:2*nsmax,1:3) :: w8ring_TQU

    integer(I4B) :: l, m, ith, scalel, scalem, nrings, nphmx
    integer(I8B) :: istart_south, istart_north, npix
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: recfac, dalm
    real(DP),     dimension(:),   allocatable :: lam_fact, mfac

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, ii
    integer(I8B), dimension(0:SMAXCHK) :: startpix
    integer(I4B), dimension(0:SMAXCHK) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it
!    real(sp) :: time0, time1, time2, time3, time4, tt1, tt2

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, kind=DP)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(0:5,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_fact')
    endif
    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    call init_rescale()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    alm_TGC = 0.0 ! set the whole alm array to zero

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

!    tt1 = 0
    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

!       call wall_clock_time(time1)
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk,nph, startpix, kphi0, w8ring_TQU, phas_ns, &
!$OMP      keep_north, keep_south, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, i, ii, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = map_TQU(istart_north:istart_north+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif

          if (ith < nrings .and. keep_south(ithl)) then
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = map_TQU(istart_south:istart_south+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo
!$OMP end do
       if (do_openmp()) then
          deallocate (ring)
       endif
!$OMP end parallel

!       call wall_clock_time(time2)
!       tt1 = tt1 + time2 - time1

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, &
!$OMP    rescale_tab, ovflow, unflow, omega_pix, &
!$OMP    cth, sth, keep_it, mfac, normal_l, alm_TGC, phas_ns, one_on_s2, c_on_s2) &
!$OMP private(recfac, dalm, lam_fact, phas_sd, phas_sdx, phm, php, status, &
!$OMP   m, fm2, normal_m, ithl, l_min, &
!$OMP   scalem, scalel, corfac, par_lm, &
!$OMP   lam_mm, lam_lm, lam_lm1m, lambda_w, lambda_x, lam_0, lam_1, lam_2, &
!$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
!$OMP   l, fl, flm1, i)

       if (do_openmp()) then
          allocate(recfac(0:1,0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac')
          allocate(dalm(0:5,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm')
          allocate(lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_fact')
       endif
!$OMP do schedule(dynamic,1)

       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
                   phas_sd( 3+2*i)   =  real(php, kind=dp) 
                   phas_sd( 3+2*i+1) = aimag(php)
                enddo
                phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = corfac * lam_mm            ! Actual lam_mm

                if (m >= l_min) then
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   ! alm_T = \int T Ylm
                   ! alm_G = \int ( - Q Wlm - i U Xlm )
                   ! alm_C = \int ( - U Wlm + i Q Xlm )
                   dalm(0:1, m) = dalm(0:1, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                   dalm(2:5, m) = dalm(2:5, m)  &
                        &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
                        &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                endif

                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
                do l = m+1, nlmax
                   par_lm   = - par_lm  ! = 3 * (-1)^(l+m)
                   lam_lm1m = lam_lm * lam_fact(l) ! must be incremented, even if not used
                   lam_lm   = lam_2 * corfac * lam_mm
                   if (l >= l_min) then

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                      dalm(0:1, l) = dalm(0:1, l) + lam_lm * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, l) = dalm(2:5, l)  &
                           &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
                           &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                   endif

                   lam_0 = lam_1 * recfac(1,l-1)
                   lam_1 = lam_2
                   lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                   if (abs(lam_2) > OVFLOW) then
                      lam_1 = lam_1*UNFLOW
                      lam_2 = lam_2*UNFLOW
                      scalel= scalel + 1
                      corfac = rescale_tab(max(scalel+scalem,RSMIN))
                   elseif (abs(lam_2) < UNFLOW) then
                      lam_1 = lam_1*OVFLOW
                      lam_2 = lam_2*OVFLOW
                      scalel= scalel - 1
                      corfac = rescale_tab(max(scalel+scalem,RSMIN))
                   endif
                   
                enddo ! loop on l
             endif ! test on cut sky
          enddo ! loop on ithl
          do l = m, nlmax
             alm_TGC(1, l, m) = alm_TGC(1, l, m) + cmplx(dalm(0, l), dalm(1, l), kind=DP) * omega_pix
             alm_TGC(2, l, m) = alm_TGC(2, l, m) + cmplx(dalm(2, l), dalm(3, l), kind=DP) * normal_l(l) * omega_pix
             alm_TGC(3, l, m) = alm_TGC(3, l, m) + cmplx(dalm(4, l), dalm(5, l), kind=DP) * normal_l(l) * omega_pix
          enddo
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (recfac,dalm,lam_fact)
       endif
!$OMP end parallel
    enddo    ! loop on chunks

!    print*,'FFT',tt1
    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate (recfac,dalm,lam_fact,ring)
    endif
    deallocate(mfac, normal_l)
    deallocate(phas_ns)
    return
  end subroutine map2alm_pol_KLOAD

  !=======================================================================
  subroutine map2alm_pol_pre1_KLOAD(nsmax, nlmax, nmmax, map_TQU,&
       &                     alm_TGC, zbounds, w8ring_TQU, plm)
    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       with precomputed Ylm_T(theta)
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(0:12*nsmax**2-1,1:3) :: map_TQU
    complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(DP),     intent(IN),  dimension(1:2)                 :: zbounds
    real(DP),     intent(IN), dimension(1:2*nsmax,1:3) :: w8ring_TQU
    real(DP),     intent(IN), dimension(0:)            :: plm

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, n_plm, i_mm
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    integer(I4B)                              :: par_lm
    real(DP)            :: lam_lm, cth_ring
    real(DP)                 :: lambda_w, lambda_x, normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: lam_fact

    integer(i4b)                              :: l_min
    real(DP),     dimension(:),   allocatable :: normal_l
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl, i, ii
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    real(DP),     dimension(0:SMAXCHK-1) :: one_on_s2, c_on_s2
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it
!    real(sp) :: time0, time1, time2, time3, time4, tt1, tt2
    integer(I4B)                       :: nphl

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, kind=DP)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(normal_l(0:nlmax),stat = status)
    call assert_alloc(status,code,'normal_l')

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(dalm(0:5,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(lam_fact(0:nlmax),stat = status)
       call assert_alloc(status,code,'lam_fact')
    endif
    !     ------------ initiate variables and arrays ----------------

    alm_TGC = 0.0 ! set the whole alm array to zero

    ! generate Polarization normalisation
    call gen_normpol(nlmax, normal_l)

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

!       call wall_clock_time(time1)
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk,nph, startpix, kphi0, w8ring_TQU, phas_ns, &
!$OMP      keep_north, keep_south, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, i, ii, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = map_TQU(istart_north:istart_north+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif

          if (ith < nrings .and. keep_south(ithl)) then
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = map_TQU(istart_south:istart_south+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo
!$OMP end do
       if (do_openmp()) then
          deallocate (ring)
       endif
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, omega_pix, n_lm, plm, ith, &
!$OMP    cth, sth, keep_it, normal_l, alm_TGC, phas_ns, one_on_s2, c_on_s2) &
!$OMP private(dalm, lam_fact, phas_sd, phas_sdx, phm, php, status, &
!$OMP   m, fm2, normal_m, ithl, l_min, par_lm, i_mm, &
!$OMP   lam_lm, lam_lm1m, lambda_w, lambda_x, &
!$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
!$OMP   l, fl, flm1, i)

       if (do_openmp()) then
          allocate(dalm(0:5,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm')
          allocate(lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_fact')
       endif
!$OMP do schedule(dynamic,1)

       do m = 0, nmmax
          ! generate Ylm relation factor for degree m
          call gen_lamfac(nlmax, m, lam_fact)
          fm2 = real(m * m, kind = DP)
          normal_m = (2.0_dp * m) * (1 - m)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
                   phas_sd( 3+2*i)   =  real(php, kind=dp) 
                   phas_sd( 3+2*i+1) = aimag(php)
                enddo
                phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                !           ---------- l = m ----------
                par_lm = 3  ! = 3 * (-1)^(l+m)
                lam_lm = plm(i_mm)

                if (m >= l_min) then
                   lambda_w =  (normal_m * lam_lm)  * ( 0.5_dp - one_on_s2(ithl) )
                   lambda_x =  (normal_m * lam_lm)  *            c_on_s2(ithl)

                   ! alm_T = \int T Ylm
                   ! alm_G = \int ( - Q Wlm - i U Xlm )
                   ! alm_C = \int ( - U Wlm + i Q Xlm )
                   dalm(0:1, m) = dalm(0:1, m) + lam_lm * phas_sd(par_lm:par_lm+1)
                   dalm(2:5, m) = dalm(2:5, m)  &
                        &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
                        &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                endif

                !           ---------- l > m ----------
                cth_ring = cth(ithl)
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
                do l = m+1, nlmax
                   par_lm = - par_lm  ! = 3 * (-1)^(l+m)
                   if (l >= l_min) then
                      lam_lm1m = plm(i_mm+l-m-1) * lam_fact(l)
                      lam_lm   = plm(i_mm+l-m)

                      fl = real(l, kind = DP)
                      flm1 = fl - 1.0_dp
                      a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                      a_x =  two_cth_ring * flm1
                      lambda_w =                b_w * lam_lm1m - a_w * lam_lm 
                      lambda_x = fm_on_s2 * (         lam_lm1m - a_x * lam_lm)

                      dalm(0:1, l) = dalm(0:1, l) + lam_lm * phas_sd(par_lm:par_lm+1)
                      dalm(2:5, l) = dalm(2:5, l)  &
                           &       - lambda_w * phas_sd(par_lm+2:par_lm+5) &
                           &       + lambda_x * phas_sdx(par_lm:par_lm+3)
                   endif

                enddo ! loop on l
             endif ! test on cut sky
          enddo ! loop on ithl
          do l = m, nlmax
             alm_TGC(1, l, m) = alm_TGC(1, l, m) + cmplx(dalm(0, l), dalm(1, l), kind=DP) * omega_pix
             alm_TGC(2, l, m) = alm_TGC(2, l, m) + cmplx(dalm(2, l), dalm(3, l), kind=DP) * normal_l(l) * omega_pix
             alm_TGC(3, l, m) = alm_TGC(3, l, m) + cmplx(dalm(4, l), dalm(5, l), kind=DP) * normal_l(l) * omega_pix
          enddo
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (dalm,lam_fact)
       endif
!$OMP end parallel
    enddo    ! loop on chunks

!    print*,'FFT',tt1
    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate (dalm,lam_fact,ring)
    endif
    deallocate(normal_l, phas_ns)
    return
  end subroutine map2alm_pol_pre1_KLOAD

  !=======================================================================
  subroutine map2alm_pol_pre2_KLOAD(nsmax, nlmax, nmmax, map_TQU,&
       &                     alm_TGC, zbounds, w8ring_TQU, plm)
    !=======================================================================
    !     computes the a(l,m) from a polarized map for the HEALPIX pixelisation
    !       with precomputed Ylm_T(theta) an dYlm_P(theta)
    !=======================================================================
    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    real(KMAP),   intent(IN),  dimension(0:12*nsmax**2-1,1:3) :: map_TQU
    complex(KALMC), intent(OUT), dimension(1:3,0:nlmax,0:nmmax) :: alm_TGC
    real(DP),     intent(IN),  dimension(1:2)                 :: zbounds
    real(DP),     intent(IN), dimension(1:2*nsmax,1:3) :: w8ring_TQU
    real(DP),     intent(IN), dimension(0:,1:)         :: plm

    integer(I4B) :: l, m, ith
    integer(I8B) :: istart_south, istart_north, npix, n_lm, n_plm, i_mm, i_up
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix

    real(DP),     dimension(-3:8)             :: phas_sd
    real(DP),     dimension(-3:6)             :: phas_sdx
    real(DP),     dimension(:,:), allocatable :: dalm, plm_sub

    integer(i4b)                              :: l_min, l_start
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:,:), allocatable :: phas_ns
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)    :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl, i, ii
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP)                           :: cth
    real(DP),     dimension(0:SMAXCHK-1) :: sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    character(LEN=*), parameter :: code = 'MAP2ALM'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, DP)
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(phas_ns(0:nlmax,0:5,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phas_ns')

    if (.not.do_openmp()) then
          allocate(dalm(0:5,0:nlmax), plm_sub(1:3,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm & plm_sub')
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
    endif
    !     ------------ initiate variables and arrays ----------------

    alm_TGC = 0.0 ! set the whole alm array to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth, zbounds, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       !-----------------------------------------------------------------------
       !  computes the integral in phi : phas_m(theta)
       !  for each parallele from Pole to Equator (use symmetry for S. Hemisphere)
       !-----------------------------------------------------------------------
       phas_ns = 0_dpc

!$OMP parallel default(none) &
!$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
!$OMP      lchk, uchk,nph, startpix, kphi0, w8ring_TQU, phas_ns, &
!$OMP      keep_north, keep_south, map_TQU) &
!$OMP  private(ithl, nphl, istart_north, istart_south, &
!$OMP      ith, ring, i, ii, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
!$OMP do schedule(dynamic,1)
       ! do Fourier Transform on rings
       do ith = lchk, uchk
          ithl = ith - lchk !local index
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          if (keep_north(ithl)) then
             do i=1,3
                ii = 2*(i-1) ! 0,2,4
                ring(0:nphl-1) = map_TQU(istart_north:istart_north+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif

          if (ith < nrings .and. keep_south(ithl)) then
             do i=1,3
                ii = 2*i - 1 ! 1,3,5
                ring(0:nphl-1) = map_TQU(istart_south:istart_south+nphl-1,i) * w8ring_TQU(ith,i)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_ns(0,ii,ithl), kphi0(ithl))
             enddo
          endif
       enddo
!$OMP end do
       if (do_openmp()) then
          deallocate (ring)
       endif
!$OMP end parallel

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, omega_pix, n_lm, ith, &
!$OMP    sth, keep_it, alm_TGC, phas_ns, plm) &
!$OMP private(dalm, phas_sd, phas_sdx, plm_sub, phm, php, status, &
!$OMP   m, ithl, l_min, l_start, l, i_mm, i_up, i)
       if (do_openmp()) then
          allocate(dalm(0:5,0:nlmax), plm_sub(1:3,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm & plm_sub')
       endif
!$OMP do schedule(dynamic,1)

       do m = 0, nmmax 

          ! introduce double precision vector to perform summation over ith for each l
          dalm(0:5, m:nlmax) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl)) ! lower l bound of non-negligible Ylm
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith
                i_up = i_mm + nlmax - m
                plm_sub(1,m:nlmax) = plm(i_mm:i_up,1)
                plm_sub(2,m:nlmax) = plm(i_mm:i_up,2)
                plm_sub(3,m:nlmax) = plm(i_mm:i_up,3)

                do i=0,2 ! loop on T, Q, U
                   phm = phas_ns(m, 2*i, ithl) - phas_ns(m, 2*i+1, ithl) ! N-S: (l+m) odd
                   phas_sd(-3+2*i)   =  real(phm, kind=dp) 
                   phas_sd(-3+2*i+1) = aimag(phm)
                   php = phas_ns(m, 2*i, ithl) + phas_ns(m, 2*i+1, ithl) ! N+S: (l+m) even
                   phas_sd( 3+2*i)   =  real(php, kind=dp) 
                   phas_sd( 3+2*i+1) = aimag(php)
                enddo
                phas_sdx(-3: 0) = (/ phas_sd(8), - phas_sd(7), - phas_sd(6),   phas_sd(5) /)
                phas_sdx( 3: 6) = (/ phas_sd(2), - phas_sd(1), - phas_sd(0),   phas_sd(-1) /)

                !           ---------- l = m ----------
                if (m >= l_min) then
                   
                   ! alm_T = \int T Ylm
                   ! alm_G = \int ( - Q Wlm - i U Xlm )
                   ! alm_C = \int ( - U Wlm + i Q Xlm )
                   dalm(0:1, m) = dalm(0:1, m) + plm_sub(1,m) * phas_sd(3:4)
                   dalm(2:5, m) = dalm(2:5, m)  &
                        &       - plm_sub(2,m) * phas_sd(5:8) &
                        &       + plm_sub(3,m) * phas_sdx(3:6)
                endif

                !           ---------- l > m ----------
                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(0:1,l) = dalm(0:1,l) + plm_sub(1,l) * phas_sd(-3:-2)
                   dalm(2:5,l) = dalm(2:5,l) - plm_sub(2,l) * phas_sd(-1:2) &
                        &                    + plm_sub(3,l) * phas_sdx(-3:0)
                enddo ! loop on l

                l_start = max(m+2, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   dalm(0:1,l) = dalm(0:1,l) + plm_sub(1,l) * phas_sd(3:4)
                   dalm(2:5,l) = dalm(2:5,l) - plm_sub(2,l) * phas_sd(5:8) &
                        &                    + plm_sub(3,l) * phas_sdx(3:6)
!                   dalm(2:5,l) = dalm(2:5,l) + plm_sub(2,l) * phas_sd(5:8) &
!                        &                    - plm_sub(3,l) * phas_sdx(3:6)
                enddo ! loop on l

             endif ! test on cut sky
          enddo ! loop on ithl
          do l = m, nlmax
             alm_TGC(1, l, m) = alm_TGC(1, l, m) + cmplx(dalm(0, l), dalm(1, l), kind=DP) * omega_pix
             alm_TGC(2, l, m) = alm_TGC(2, l, m) + cmplx(dalm(2, l), dalm(3, l), kind=DP) * omega_pix
             alm_TGC(3, l, m) = alm_TGC(3, l, m) + cmplx(dalm(4, l), dalm(5, l), kind=DP) * omega_pix
          enddo
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (dalm, plm_sub)
       endif
!$OMP end parallel
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    deallocate(phas_ns)
    if (.not.do_openmp()) then
        deallocate (ring, dalm, plm_sub)
    endif
    return
  end subroutine map2alm_pol_pre2_KLOAD

  !**************************************************************************
  !
  !             ALM treatment
  !
  !**************************************************************************

  !=======================================================================
  subroutine alter_alm_KLOAD(nsmax, nlmax, nmmax, fwhm_arcmin, alm, beam_file, window)
    !=======================================================================
    !     multiply the alm by a gaussian function of FWHM fwhm_arcmin, 
    !        or a beam_file or an arbitrary window
    !     --> equivalent to a convolution in real sphere
    !=======================================================================
    integer(I4B),     intent(in)                         :: nlmax,nsmax,nmmax
    real(KALM),       intent(in)                         :: fwhm_arcmin
    complex(KALMC),   intent(inout), dimension(1:,0:,0:) :: alm
    character(LEN=*), intent(IN),  optional              :: beam_file
    real(KALM),       intent(IN),  optional, dimension(0:,1:) :: window

    real(DP), dimension(:,:),  allocatable :: beamw
    integer(I4B)                           :: status

    integer(I4B)                           :: m, i, nl, nda, nlw, ndw, j
    character(len=*), parameter            :: code = "alter_alm"

    !-----------------------------------------------------------------------

    nda = size(alm,1)

    if (present(window)) then
       nlw = size(window,1)
       ndw = size(window,2)
       nl  = min(nlw, nlmax+1)
       do m=0, nmmax
          do i = 1, nda
             j = min(ndw, i) ! if no polarization window, replicate temperature
             alm(i, m:nl-1, m) = alm(i, m:nl-1, m) * window(m:nl-1, j)
          enddo
       enddo
       ! set to 0 alms not altered
       if (nlw <= nlmax) then
          alm(1:nda, nlw:nlmax, 0:nmmax) = 0.0_KALMC
          print*,code //' set to 0 alm with l in range ',nlw,nlmax
       endif
    else
       allocate(beamw(0:nlmax,1:nda), stat = status)
       call assert_alloc(status,code,'beamw')
       call generate_beam(real(fwhm_arcmin,kind=dp), nlmax, beamw, beam_file)
       do m = 0, nmmax
          do i = 1, nda
             alm(i, m:nlmax, m) = alm(i, m:nlmax, m) * beamw(m:nlmax, i)
          enddo
       enddo
       deallocate(beamw)
    endif

    return
  end subroutine alter_alm_KLOAD

  !=======================================================================
  subroutine create_alm_v12_KLOAD &
       &     (nsmax, nlmax, nmmax, polar, filename, iseed, fwhm_arcmin, &
       &        alm_TGC, header_PS, windowfile, units, beam_file)
    !=======================================================================
    use rngmod, only: rand_init, planck_rng
    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax

    INTEGER(I4B),       INTENT(IN) :: polar
    CHARACTER(LEN=*),   INTENT(IN)      :: filename
    INTEGER(I4B),       INTENT(INOUT)   :: iseed
    REAL(KALM),           INTENT(IN)      :: fwhm_arcmin
    COMPLEX(KALMC),       INTENT(OUT), &
         & DIMENSION(1:1+2*polar,0:nlmax,0:nmmax) :: alm_TGC
    CHARACTER(LEN=80),  INTENT(OUT),             DIMENSION(1:):: header_PS
    CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL               :: windowfile
    CHARACTER(LEN=80),  INTENT(OUT),   OPTIONAL, DIMENSION(1:):: units
    CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL               :: beam_file
    !
    type(planck_rng)       :: rng_handle


    print*,'============================================================================='
    print*,'WARNING: create_alm calling sequence has changed'
    print*,' from'
    print*,'  call create_alm (nsmax, nlmax, nmmax, polar, filename, ISEED, fwhm_arcmin,'
    print*,'                    alm_TGC, header_PS [, windowfile, units, beam_file])'
    print*,' to'
    print*,'  call create_alm (nsmax, nlmax, nmmax, polar, filename, RNG_HANDLE, fwhm_arcmin,'
    print*,'                    alm_TGC, header_PS [, windowfile, units, beam_file])'
    print*,'  '
    print*,' see documentation for details'
    print*,'============================================================================='
    if (iseed < 0) then
       call rand_init(rng_handle, iseed) ! start new sequence
       call  create_alm &
            &     (nsmax, nlmax, nmmax, polar, filename, rng_handle, fwhm_arcmin, &
            &        alm_TGC, header_PS, windowfile, units, beam_file)
       iseed = abs(iseed)
    else
       print*,'ERROR: old calling sequence can only be used with a new seed (ISEED < 0).'
       print*,' see documentation for details on new interface'
       call fatal_error('create_alm_v12')
    endif

    return
  end subroutine create_alm_v12_KLOAD
   !=======================================================================
  subroutine create_alm_KLOAD &
       &     (nsmax, nlmax, nmmax, polar, filename, rng_handle, fwhm_arcmin, &
       &        alm_TGC, header_PS, windowfile, units, beam_file)
    !=======================================================================
    !     creates the a_lm from the power spectrum,
    !     assuming they are gaussian  and complex
    !     with a variance given by C(l)
    !
    !     the input file FILENAME should contain :
    !       l, C_temp(l), [ C_grad(l), C_curl(l), C_temp_grad(l) ]
    !
    !     with *consecutive* l's (missing C(l) are set to 0.)
    !
    !     because the map is real we have : a_l-m = (-)^m conjug(a_lm)
    !     so we actually compute them only for m >= 0
    !
    !
    !     RNG_HANDLE : planck_rng structure
    !
    !     FWHM_ARCMIN
    !
    !     ALM_TGC (Complex array) : either (1:1, 0:nlmax, 0:nmmax)
    !       for temperature only (if POLAR=0)
    !                               or     (1:3, 0:nlmax, 0:nmmax)
    !       for temperature + polarisation (if POLAR =1)
    !       (respectively Temp, Grad or Electric component
    !        and Curl or Magnetic one)
    !
    !=======================================================================
    use fitstools, only: fits2cl
!    use ran_tools, only: randgauss_boxmuller
    use rngmod, only: rand_gauss, planck_rng
    USE head_fits, ONLY : add_card, merge_headers, get_card

    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B),       INTENT(IN) :: polar
    type(planck_rng),   intent(inout)                         :: rng_handle
    CHARACTER(LEN=*),   INTENT(IN)      :: filename
    REAL(KALM),           INTENT(IN)      :: fwhm_arcmin
    COMPLEX(KALMC),       INTENT(OUT), &
         & DIMENSION(1:1+2*polar,0:nlmax,0:nmmax) :: alm_TGC
    CHARACTER(LEN=80),  INTENT(OUT),             DIMENSION(1:):: header_PS
    CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL               :: windowfile
    CHARACTER(LEN=80),  INTENT(OUT),   OPTIONAL, DIMENSION(1:):: units
    CHARACTER(LEN=*),   INTENT(IN),    OPTIONAL               :: beam_file

    INTEGER(I4B) :: i, l, m, npw, l_max
    INTEGER(I4B) :: ncl, nlheader, count_tt, count_pn
    INTEGER(I4B) :: status
    REAL(DP) ::  hsqrt2, quadrupole
    REAL(DP) ::  rms_tt, rms_g1, rms_g2, rms_cc
    REAL(DP) ::  zeta1_r, zeta1_i, zeta2_r, zeta2_i, zeta3_r, zeta3_i

    LOGICAL(LGT) ::  polarisation

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: pixlw, beamw
    REAL(DP), DIMENSION(0:nlmax)          :: cls_tt, cls_tg, cls_gg, cls_cc
    REAL(DP), DIMENSION(0:nlmax,1:4)      :: cl_in
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: fact_norm

    CHARACTER(LEN=*), PARAMETER :: code = 'CREATE_ALM'
    CHARACTER(LEN=*), PARAMETER :: version = '2.0.0'
    CHARACTER(LEN=20)                            ::  string_quad
    CHARACTER(LEN=80), DIMENSION(1:180)          :: header_file
    CHARACTER(LEN=80), DIMENSION(:), allocatable :: units_power
    CHARACTER(LEN=80) :: temptype, com_tt
    CHARACTER(LEN=80) :: polnorm, com_pn

    !=======================================================================
    ! Is polarisation required ?
    if (polar .eq. 0) then
       polarisation = .false.
    elseif (polar .eq. 1) then
       polarisation = .true.
    else
       print*,'wrong choice of polar'
       call fatal_error
    endif

    ! set maximum multipole (depends on user choice and pixel window function)
    l_max = nlmax
    npw = 4*nsmax + 1
    if (present(windowfile)) then
       l_max = min(l_max, npw - 1)
    else
       npw = l_max + 1
    endif
    if (l_max < nlmax) then
       print*,'a_lm are only generated for 0 < l <= ',l_max
    endif
    !-----------------------------------------------------------------------

    !     ----- set the alm array to zero ------
    if (polarisation) then
       alm_TGC(1:3, 0:nlmax, 0:nmmax) = CMPLX(0.0,0.0, kind=KALM)
    else
       alm_TGC(1  , 0:nlmax, 0:nmmax) = CMPLX(0.0,0.0, kind=KALM)
    endif
    !     --------------------------------------
    cls_tt = 0.0_dp
    cls_gg = 0.0_dp
    cls_cc = 0.0_dp
    cls_tg = 0.0_dp
    cl_in  = 0.0_dp

    !     --- reads the C(l) file ---
    ncl = 1
    if (polarisation) ncl = 4
    nlheader = SIZE(header_file)
    allocate(units_power(1:ncl), stat = status)
    call assert_alloc(status,code,'units_power')

    call fits2cl(filename, cl_in(0:l_max,1:ncl), l_max, ncl, &
         &       header_file, units=units_power)

    !    ---- check input power spectra consistency  ----
    do i=1,ncl
       do l=0,l_max
          if (i < 4) then ! auto correlation
             if (cl_in(l,i) < 0.0) then
                print*,code,'> Negative input power spectrum at l =', l,', index = ',i
                print*,code,'> ',cl_in(l,i)
                call fatal_error
             endif
          else ! cross correlation
             if (abs(cl_in(l,4)) > sqrt(cl_in(l,1))*sqrt(cl_in(l,2)) ) then
                print *,code,'> Inconsistent T, E and ExT ower spectrum terms at l = ',l
                print*,code,'> ',cl_in(l,1),cl_in(l,2),cl_in(l,4)
             endif
          endif
       enddo
    enddo

    allocate( fact_norm(0:nlmax),stat = status)
    call assert_alloc(status,code,'fact_norm')
    do l=1, nlmax
       fact_norm(l) = 1.0_dp         !   Normalisation=1, 'pure' C_ls are expected
    enddo
    fact_norm(0) = 1.0_dp

    !     we always expect 4 columns, in the order T, G(=E), C(=B) and TG(=TE)
    cls_tt(0:l_max) = cl_in(0:l_max,1) / fact_norm(0:l_max)
    cls_gg(0:l_max) = cl_in(0:l_max,2) / fact_norm(0:l_max)
    cls_cc(0:l_max) = cl_in(0:l_max,3) / fact_norm(0:l_max)
    cls_tg(0:l_max) = cl_in(0:l_max,4) / fact_norm(0:l_max)
    deallocate(fact_norm)

    ! converts units for power spectrum into units for maps
    if (present(units)) then
       call pow2alm_units(units_power, units)
    endif

    ! finds out temperature type (THERMO or ANTENNA)
    ! THERMO is the default
    call get_card(header_file,"TEMPTYPE",temptype,com_tt,count=count_tt)
    if (count_tt < 1) then
       temptype = "THERMO"
       com_tt = " temperature : either THERMO or ANTENNA"
       print*,'   Will assume '//temptype
    endif

    if (polar > 0) then
       call get_card(header_file,"POLNORM",polnorm,com_pn,count=count_pn)
       if (count_pn < 1) then
          print*,' ********************************************************* '
          print*,'            The convention for normalization '
          print*,'        of the input polarization power spectra '
          print*,'                      is unknown. '
          print*,' The code will proceed as if the convention was that of CMBFAST'
          print*,'  See the Healpix PDF of Html documentation for details'
          print*,' ********************************************************* '
       endif
    endif

    quadrupole = cls_tt(2)
    !     --- creates the header relative to power spectrum ---
    header_PS = ''
    call add_card(header_PS)
    call add_card(header_PS,'HISTORY',' alm generated by ' &
         &        //code//'_'//version//' from following power spectrum')
    call add_card(header_PS)
    call add_card(header_PS,'COMMENT','----------------------------------------------------')
    call add_card(header_PS,'COMMENT','Planck Power Spectrum Description Specific Keywords')
    call add_card(header_PS,'COMMENT','----------------------------------------------------')
    call add_card(header_PS,'COMMENT','Input power spectrum in : ')
    call add_card(header_PS,'COMMENT',TRIM(filename))
    call add_card(header_PS)
    call add_card(header_PS,'COMMENT','Quadrupole')
    write(string_quad,'(1pe15.6)') quadrupole
    call add_card(header_PS,'COMMENT','  C(2) = '//string_quad//' '//trim(units_power(1)))
    call add_card(header_PS)
    call add_card(header_PS,"TEMPTYPE",temptype,com_tt)
    call add_card(header_PS)
    ! ------- insert header read in power spectrum file -------
    call merge_headers(header_file, header_PS)
    deallocate(units_power)

    !     --- smoothes the initial power spectrum ---
    !       beam (gaussian or external file) + pixel (external file)

    ! define beam
    allocate(beamw(0:l_max,1:3), stat = status)
    call assert_alloc(status,code,'beamw')
    call generate_beam(real(fwhm_arcmin,kind=dp), l_max, beamw, beam_file)
!     call gaussbeam(real(fwhm_arcmin,kind=dp), l_max, beamw)

    ! get the pixel window function
    allocate(pixlw(0:npw-1,1:3), stat = status)
    call assert_alloc(status,code,'pixlw')
    if(present(windowfile)) then
       call pixel_window(pixlw, windowfile=windowfile)
    else
       pixlw(:,:)  = 1.0_dp
    endif
    ! multiply beam and pixel windows
    beamw(0:l_max,1:3) = beamw(0:l_max,1:3) * pixlw(0:l_max,1:3)

    cls_tt(0:l_max) = cls_tt(0:l_max) * beamw(0:l_max,1)**2
    cls_tg(0:l_max) = cls_tg(0:l_max) * beamw(0:l_max,1)*beamw(0:l_max,2)
    cls_gg(0:l_max) = cls_gg(0:l_max) * beamw(0:l_max,2)**2
    cls_cc(0:l_max) = cls_cc(0:l_max) * beamw(0:l_max,3)**2
    deallocate(pixlw)
    deallocate(beamw)


    !     --- generates randomly the alm according to their power spectrum ---
    !     alm_T = zeta1 * rms_tt
    !     alm_G = zeta1 * rms_g1 + zeta2 * rms_g2
    !     alm_C = zeta3 * rms_cc
    hsqrt2 = SQRT2 / 2.0_dp

    do l = 0, l_max
       rms_tt = 0.0_dp
       rms_g1 = 0.0_dp
       if (cls_tt(l) .ne. 0) then
          rms_tt   = sqrt(cls_tt(l))
          rms_g1   = cls_tg(l) / rms_tt
       endif

       !        ------ m = 0 ------
!        zeta1_r = randgauss_boxmuller(iseed)  ! gaussian deviate based on ran_mwc
       zeta1_r = rand_gauss(rng_handle)
       zeta1_i = 0.0_dp
       alm_TGC(1, l, 0)   = CMPLX(zeta1_r, zeta1_i, kind=KALM) * rms_tt ! T
       if (polarisation) then
          alm_TGC(2, l, 0) = CMPLX(zeta1_r, zeta1_i, kind=KALM) * rms_g1 ! G
       endif

       !        ------ m > 0 ------
       do m = 1,l
!          zeta1_r = randgauss_boxmuller(iseed) * hsqrt2
!          zeta1_i = randgauss_boxmuller(iseed) * hsqrt2
          zeta1_r = rand_gauss(rng_handle) * hsqrt2
          zeta1_i = rand_gauss(rng_handle) * hsqrt2
          alm_TGC(1, l, m) = CMPLX(zeta1_r, zeta1_i, kind=KALM) * rms_tt
          if (polarisation) then
             alm_TGC(2, l, m) = CMPLX(zeta1_r, zeta1_i, kind=KALM) * rms_g1
          endif
       enddo
    enddo

    !     the coefficient generation is separated so that the Temperature
    !     coeff. are the same (for the same seed) weither the polarisation is set or not

    if (polarisation) then
       do l = 0, l_max
          rms_g2 = 0.0_dp
          rms_cc = 0.0_dp
          if (cls_tt(l) .ne. 0) then
             rms_g2 = cls_gg(l) - (cls_tg(l)/cls_tt(l))*cls_tg(l) ! to avoid underflow
             ! test for consistency but make sure it is not due to round off error
             if (rms_g2 <= 0.0) then
                if (abs(rms_g2) > abs(1.e-8*cls_gg(l))) then
                   print*,code,'> Inconsistent TT, GG and TG spectra at l=',l
                   call fatal_error
                else ! only round off error, keep going
                   rms_g2 = 0.0
                endif
             endif
             rms_g2 = SQRT( rms_g2 )
             rms_cc = SQRT( cls_cc(l) )
          endif

          !           ------ m = 0 ------
!          zeta2_r = randgauss_boxmuller(iseed)
          zeta2_r = rand_gauss(rng_handle)
          zeta2_i = 0.0_dp
!          zeta3_r = randgauss_boxmuller(iseed)
          zeta3_r = rand_gauss(rng_handle)
          zeta3_i = 0.0_dp
          alm_TGC(2, l, 0) = alm_TGC(2, l, 0) &
               &           + CMPLX(zeta2_r, zeta2_i, kind=KALM) * rms_g2 ! G
          alm_TGC(3, l, 0) = CMPLX(zeta3_r, zeta3_i, kind=KALM) * rms_cc ! C

          !           ------ m > 0 ------
          do m = 1,l
!              zeta2_r = randgauss_boxmuller(iseed) * hsqrt2
!              zeta2_i = randgauss_boxmuller(iseed) * hsqrt2
!              zeta3_r = randgauss_boxmuller(iseed) * hsqrt2
!              zeta3_i = randgauss_boxmuller(iseed) * hsqrt2
             zeta2_r = rand_gauss(rng_handle) * hsqrt2
             zeta2_i = rand_gauss(rng_handle) * hsqrt2
             zeta3_r = rand_gauss(rng_handle) * hsqrt2
             zeta3_i = rand_gauss(rng_handle) * hsqrt2
             alm_TGC(2, l, m) = alm_TGC(2, l, m) &
                  &                 + CMPLX(zeta2_r, zeta2_i, kind=KALM) * rms_g2 ! G
             alm_TGC(3, l, m) = CMPLX(zeta3_r, zeta3_i, kind=KALM) * rms_cc ! C
          enddo
       enddo
    endif

    return
  end subroutine create_alm_KLOAD

  !========================================================
  subroutine alm2cl2_KLOAD(nlmax, nmmax, alm1, alm2, cl)
  !========================================================
    ! computes C(l) from a_lm, in the order
    ! TT, [EE, TE, BB, [TB, EB]]
    !=======================================================
    integer(I4B),                      intent(in) :: nlmax, nmmax
    complex(KALMC), dimension(1:,0:,0:), intent(in) :: alm1, alm2
    real(KALM)    , dimension(0:, 1: ),  intent(out):: cl
    ! 
    integer(I4B) :: l, ncl, na1, na2, mm
    complex(DPC)     :: dc
    real(DP), parameter :: two = 2.000000000000000000_dp
    real(DP), parameter :: one = 1.000000000000000000_dp
    logical(LGT) :: polarisation, bcoupling
    !========================================================

    ncl = size(cl, 2)
    na1 = size(alm1, 1)
    na2 = size(alm2, 1)
    polarisation = (na1 >= 3 .and. na2 >= 3 .and. ncl >=4)
    bcoupling    = (ncl >=6) .and. polarisation
    cl = 0.0_KALM

    ! TT power spectrum
    do l = 0, nlmax
       mm = min(l, nmmax)
       dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(1,l,1:mm)))
       dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(1,l,0)
       cl(l,1) = real(dc, kind=DP) / (two*l + one)
    enddo

    if (polarisation) then
       ! GG or EE power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(2,l,1:mm)*conjg(alm2(2,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm2(2,l,0)
          cl(l,2) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! CC or BB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(3,l,1:mm)*conjg(alm2(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(3,l,0)   *      alm2(3,l,0)
          cl(l,3) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! TG or TE power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(2,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(2,l,0)
          cl(l,4) = real(dc, kind=DP) / (two*l + one)
       enddo
    endif 

    if (bcoupling) then
       ! TC or TB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(3,l,0)
          cl(l,5) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! GC or EB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(2,l,1:mm)*conjg(alm2(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm2(3,l,0)
          cl(l,6) = real(dc, kind=DP) / (two*l + one)
       enddo
    endif

    return
  end subroutine alm2cl2_KLOAD
  !========================================================
  subroutine alm2cl1_KLOAD(nlmax, nmmax, alm1, cl)
  !========================================================
    ! computes C(l) from a_lm, in the order
    ! TT, [EE, TE, BB, [TB, EB]]
    !=======================================================
    integer(I4B),                      intent(in) :: nlmax, nmmax
    complex(KALMC), dimension(1:,0:,0:), intent(in) :: alm1
    real(KALM)    , dimension(0:, 1: ),  intent(out):: cl
    ! 
    integer(I4B) :: l, ncl, na1, na2, mm
    complex(DPC)     :: dc
    real(DP), parameter :: two = 2.000000000000000000_dp
    real(DP), parameter :: one = 1.000000000000000000_dp
    logical(LGT) :: polarisation, bcoupling
    !========================================================

    ncl = size(cl, 2)
    na1 = size(alm1, 1)
    polarisation = (na1 >= 3 .and. ncl >=4)
    bcoupling    = (ncl >=6) .and. polarisation
    cl = 0.0_KALM

    ! TT power spectrum
    do l = 0, nlmax
       mm = min(l, nmmax)
       dc =          sum(      alm1(1,l,1:mm)*conjg(alm1(1,l,1:mm)))
       dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm1(1,l,0)
       cl(l,1) = real(dc, kind=DP) / (two*l + one)
    enddo

    if (polarisation) then
       ! GG or EE power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(2,l,1:mm)*conjg(alm1(2,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm1(2,l,0)
          cl(l,2) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! CC or BB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(3,l,1:mm)*conjg(alm1(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(3,l,0)   *      alm1(3,l,0)
          cl(l,3) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! TG or TE power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(1,l,1:mm)*conjg(alm1(2,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm1(2,l,0)
          cl(l,4) = real(dc, kind=DP) / (two*l + one)
       enddo
    endif 

    if (bcoupling) then
       ! TC or TB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(1,l,1:mm)*conjg(alm1(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm1(3,l,0)
          cl(l,5) = real(dc, kind=DP) / (two*l + one)
       enddo
       
       ! GC or EB power spectrum
       do l = 0, nlmax
          mm = min(l, nmmax)
          dc =          sum(      alm1(2,l,1:mm)*conjg(alm1(3,l,1:mm)))
          dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm1(3,l,0)
          cl(l,6) = real(dc, kind=DP) / (two*l + one)
       enddo
    endif

    return
  end subroutine alm2cl1_KLOAD

  !========================================================
  subroutine rotate_alm_KLOAD(lmax, alm, psi, theta, phi)
    !=========================================================
    !Input: Complex array alm(p,l,m) with (l,m) in [0,lmax]^2, and p in [1,nd]
    !Euler rotation angles psi, theta, phi in radians
    !Output: Rotated array alm(p, l,m)
    !
    ! Euler angle convention  is right handed, active rotation
    ! psi is the first rotation about the z-axis (vertical), in [-2pi,2pi]
    ! then theta about the ORIGINAL (unrotated) y-axis, in [-2pi,2pi]
    ! then phi  about the ORIGINAL (unrotated) z-axis (vertical), in [-2pi,2pi]
    !
    ! Equivalently
    ! phi is the first rotation about the z-axis (vertical)
    ! then theta  about the NEW   y-axis (line of nodes)
    ! then psi    about the FINAL z-axis (figure axis)
    ! ---
    ! the recursion on the Wigner d matrix is inspired from the very stable
    ! double sided one described in Risbo (1996, J. of Geodesy, 70, 383)
    ! based on equation (4.4.1) in Edmonds (1957).
    ! the Risbo's scatter scheme has been repladed by a gather scheme for
    ! better computing efficiency
    ! the size of the matrix is divided by 2 using Edmonds Eq.(4.2.5) 
    ! to speed up calculations
    ! the loop on j has been unrolled for further speed-up
    ! EH, March--April 2005
    !=========================================================
    integer(I4B),   intent(in) :: lmax
    complex(KALMC), intent(inout), dimension(1:,0:,0:) :: alm
    real(DP),       intent(in) :: psi, theta, phi
    ! local variables
    complex(DPC), dimension(0:lmax) :: exppsi, expphi
    complex(DPC), dimension(:,:), allocatable :: alm1, alm2
    real(DP),     dimension(:,:), allocatable :: d, dd
    real(DP),     dimension(:),   allocatable :: sqt, rsqt
    real(DP),     dimension(:),   allocatable :: tsign
    integer(I4B) :: status
    integer(I4B) :: mm, ll, na1, na2, nd
    integer(I4B) :: i, j, k, kd, hj
    real(DP)     :: p, q, pj, qj, fj, temp
    character(len=*), parameter :: code = 'ROTATE_ALM'
    !==========================================================
    
    if (abs(psi) > 2.d0*PI .or. abs(phi) > 2.d0*PI .or. abs(theta) > 2.d0*PI) then
       write(*,'(a,3(g12.4))') code,psi,theta,phi
       call fatal_error(code//': angles should be in Radians')
    endif

    nd = size(alm,1)
    na1 = size(alm,2)
    na2 = size(alm,3)
    if (na1 < (lmax+1) .or. na2 < (lmax+1)) then
       call fatal_error(code//': unconsistent alm array size and lmax')
    endif

    allocate(d (-1:2*lmax,   -1:lmax),   stat = status)
    call assert_alloc(status,code,'d')
    allocate(dd(-1:2*lmax, -1:lmax), stat = status)
    call assert_alloc(status,code,'dd')
    allocate(sqt(0:2*lmax), rsqt(0:2*lmax), stat = status)
    call assert_alloc(status,code,'sqt & rsqt')
    allocate(alm1(1:nd,0:lmax), alm2(1:nd,0:lmax), stat = status)
    call assert_alloc(status,code,'alm1 & alm2')
    allocate(tsign(0:lmax+1), stat = status)
    call assert_alloc(status,code,'tsign')
    
    do i=0, lmax,2
       tsign(i)   =  1.0_dp
       tsign(i+1) = -1.0_dp
    enddo
    !     initialization of square-root  table
    do i=0,2*lmax
       sqt(i) = SQRT(DBLE(i))
    enddo

    ! initialisation of exponential table
    exppsi(0)=cmplx(1, 0, kind=DPC)
    expphi(0)=cmplx(1, 0, kind=DPC)

    do i=1,lmax
       exppsi(i)= cmplx(cos(psi*i), -sin(psi*i), kind=DPC)
       expphi(i)= cmplx(cos(phi*i), -sin(phi*i), kind=DPC)
    enddo

    ! Note: theta has the correct sign.
    p = sin(theta/2.d0)
    q = cos(theta/2.d0)

    d  = 0.0_dp ! very important for gather scheme
    dd = 0.0_dp
    do ll=0,lmax

       ! ------ build d-matrix of order l ------
       if (ll == 0) then
          d(0,0) = 1.d0
          goto 2000
       endif
       if (ll == 1) then
          !     initialize d-matrix degree 1/2
          dd(0,0)  =  q
          dd(1,0)  = -p
          dd(0,1)  =  p
          dd(1,1)  =  q
          goto 1000
       endif

       !  l - 1 --> l - 1/2
       j = 2*ll - 1
       rsqt(0:j) = sqt(j:0:-1)
       fj = DBLE(j)
       qj = q / fj
       pj = p / fj
!$OMP parallel default(none) &
!$OMP   shared(j, fj, d, dd, rsqt, sqt, q, p, qj, pj) &
!$OMP   private(k)
!$OMP do schedule(dynamic,100)
       do k = 0, j/2 ! keep only m' <= -1/2
          dd(0:j,k) = rsqt(0:j) * ( d(0:j,k)      * (sqt(j-k)  * qj)   &
               &                  + d(0:j,k-1)    * (sqt(k)    * pj) ) &
               &    +  sqt(0:j) * ( d(-1:j-1,k-1) * (sqt(k)    * qj)   &
               &                  - d(-1:j-1,k)   * (sqt(j-k)  * pj) )
       enddo ! loop on k
!$OMP end do
!$OMP end parallel
       ! l=half-integer, reconstruct m'= 1/2 by symmetry
       hj = ll-1
       if (mod(ll,2) == 0) then
          do k = 0, j-1, 2
             dd(k,   ll) =   dd(j-k,   hj)
             dd(k+1, ll) = - dd(j-k-1, hj)
          enddo
       else
          do k = 0, j-1, 2
             dd(k,   ll) = - dd(j-k,   hj)
             dd(k+1, ll) =   dd(j-k-1, hj)
          enddo
       endif

1000   continue

       !  l - 1/2 --> l
       j = 2*ll
       rsqt(0:j) = sqt(j:0:-1)
       fj = DBLE(j)
       qj = q / fj
       pj = p / fj
!$OMP parallel default(none) &
!$OMP   shared(j, fj, d, dd, rsqt, sqt, q, p, qj, pj) &
!$OMP   private(k)
!$OMP do schedule(dynamic,100)
       do k = 0, j/2 ! keep only m' <= 0
          d (0:j,k) = rsqt(0:j) * ( dd(0:j,k)      * (sqt(j-k)  * qj)   &
               &                  + dd(0:j,k-1)    * (sqt(k)    * pj) ) &
               &    +  sqt(0:j) * ( dd(-1:j-1,k-1) * (sqt(k)    * qj)   &
               &                  - dd(-1:j-1,k)   * (sqt(j-k)  * pj) )
       enddo ! loop on k
!$OMP end do
!$OMP end parallel

2000   continue
       ! ------- apply rotation matrix -------
       do kd = 1, nd
          alm1(kd,0:ll)  = alm(kd,ll,0:ll) * exppsi(0:ll)
       enddo

       ! m = 0
       do kd = 1, nd
          alm2(kd,0:ll) = alm1(kd,0) * d(ll:2*ll,ll)
       enddo

!$OMP parallel default(none) &
!$OMP   shared(d, alm1, alm2, tsign, nd, ll) &
!$OMP   private(mm, kd)
!$OMP do schedule(dynamic,100)
       do mm = 0, ll
          do kd = 1, nd
             alm2(kd, mm) = alm2(kd,mm) + sum(alm1(kd,1:ll) *                d(ll-1:0:-1,ll-mm)) &
                  &                +conjg(sum(alm1(kd,1:ll) * (tsign(1:ll) * d(ll+1:2*ll,ll-mm))))
          enddo
       enddo
!$OMP end do
!$OMP end parallel

       ! new alm for ll
       do kd = 1,nd
          alm(kd,ll,0:ll) = alm2(kd,0:ll)*expphi(0:ll)
       enddo

    enddo ! loop on ll

    deallocate(d)
    deallocate(dd)
    deallocate(sqt, rsqt)
    deallocate(alm1, alm2)

  end subroutine rotate_alm_KLOAD

