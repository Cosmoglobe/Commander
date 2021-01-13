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
module alm_tools
  !   Scalar+Open_MP implementation
  !
  !   function do_opemp
  !   subroutine init_rescale
  !   subroutine get_pixel_layout
  !   subroutine select_rings
  !   subroutine gen_recfac
  !   subroutine gen_lamfac
  !   subroutine gen_lamfac_der
  !   subroutine gen_mfac
  !   subroutine compute_lam_mm
  !   subroutine do_lam_lm
  !   subroutine do_lam_lm_pol
  !   subroutine gen_normpol
  !   function l_min_ylm
  !   subroutine warning_oldbounds
  !
  !   subroutine ring_synthesis
  !   subroutine ring_analysis
  !
  !   -------------------- in include files (see alm_map_template.f90) ---------
  !   subroutine alm2map_sc
  !   subroutine alm2map_sc_pre
  !   subroutine alm2map_pol
  !   subroutine alm2map_pol_pre1
  !   subroutine alm2map_pol_pre2
  !   subroutine alm2map_sc_der
  !   subroutine alm2map_pol_der
  !
  !   subroutine map2alm_sc
  !   subroutine map2alm_sc_pre
  !   subroutine map2alm_pol
  !   subroutine map2alm_pol_pre1
  !   subroutine map2alm_pol_pre2
  !
  !   subroutine alter_alm
  !   subroutine create_alm
  !   subroutine alm2cl
  !
  !   subroutine rotate_alm
  !   ------------------------------------------------------------------
  !
  !   subroutine plm_gen
  !
  !   subroutine pow2alm_units
  !   subroutine generate_beam
  !   subroutine gaussbeam
  !   subroutine pixel_window
  !
  use misc_utils, only: assert_alloc, assert_present, assert, fatal_error, string, strupcase !, wall_clock_time
  use healpix_types
  use healpix_fft, only: real_fft2, planck_fft2_plan, make_fft2_plan, destroy_fft2_plan, &
       &                 fft2_forward, fft2_backward
  IMPLICIT none

  ! keep everything private unless stated otherwise
  private
  !--------------------------------------------------------------------------
  ! define large and small numbers used to renormalise the recursion on the Legendre Polynomials
  integer(I4B),      private, parameter :: LOG2LG   = 100
  real(KIND=DP),     private, parameter :: FL_LARGE = 2.0_dp **   LOG2LG
  real(KIND=DP),     private, parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
  ! declare array for dynamic rescaling of the Ylm
  integer(kind=i4b), private, parameter :: RSMAX = 20, RSMIN = -20
  real(dp),          private, dimension(RSMIN:RSMAX) :: rescale_tab
  real(DP),          private, parameter :: ALN2_INV = 1.4426950408889634073599246810_dp ! 1/log(2)
  ! misc
  integer(i4b),      private, parameter :: SMAXCHK = 50 ! maximum size of chunk (in number of ring pairs)
  ! parameters of Ylm short-cut
  integer(kind=i4b), private, parameter :: HPX_MXL0 = 40 ! minimum used, choose <=0 to do full loop
  real   (kind=dp),  private, parameter :: HPX_MXL1 = 1.35_dp
  !--------------------------------------------------------------------------

  ! make (front end) routines public
  public :: alm2map, map2alm, alm2map_der
  public :: alter_alm, create_alm, alm2cl, rotate_alm
  public :: plm_gen
  public :: ring_synthesis, ring_analysis
  public :: generate_beam, gaussbeam, pixel_window, pow2alm_units

  interface alm2cl
     module procedure alm2cl2_s, alm2cl2_d, alm2cl1_s, alm2cl1_d
  end interface
 
  interface rotate_alm
     module procedure rotate_alm_s, rotate_alm_d
  end interface
 
  interface alter_alm
     module procedure alter_alm_s, alter_alm_d
  end interface
 
  interface create_alm
     module procedure create_alm_s, create_alm_d, create_alm_v12_s, create_alm_v12_d
  end interface
 
  interface alm2map_der
     module procedure alm2map_sc_der_s, alm2map_sc_der_d, alm2map_pol_der_s, alm2map_pol_der_d
  end interface

  interface alm2map
     ! Call the correct alm2map routine, depending on whether polarisation is included
     ! (determined from the rank of the map_TQU array) or precomputed plms
     ! (scalar or tensor, determined from the rank of the plm array) are included,
     ! or whether the map and alm are single or double precision
     module procedure alm2map_sc_s, alm2map_sc_pre_s, &
          &           alm2map_pol_s, alm2map_pol_pre1_s, alm2map_pol_pre2_s, &
          &           alm2map_sc_d, alm2map_sc_pre_d, &
          &           alm2map_pol_d, alm2map_pol_pre1_d, alm2map_pol_pre2_d
  end interface

  interface map2alm
     ! Call the correct map2alm routine, depending on whether polarisation is included
     ! (determined from the rank of the map_TQU array) or precomputed plms
     ! (scalar or tensor, determined from the rank of the plm array) are included,
     ! or whether the map and alm are single or double precision
     module procedure map2alm_old_sc_s,  map2alm_old_pol_s,  map2alm_old_pol2_s, &
          &           map2alm_sc_s, map2alm_sc_pre_s, &
          &           map2alm_pol_s, map2alm_pol_pre1_s, map2alm_pol_pre2_s, &
          &           map2alm_old_sc_d,  map2alm_old_pol_d,  map2alm_old_pol2_d, &
          &           map2alm_sc_d, map2alm_sc_pre_d, &
          &           map2alm_pol_d, map2alm_pol_pre1_d, map2alm_pol_pre2_d
  end interface

  ! make routines public as most of them are called by mpi_alm* routines
  public :: do_openmp
  public :: init_rescale
  public :: do_lam_lm, do_lam_lm_pol
  ! needed by mpi_alm*
  public :: l_min_ylm, select_rings
  public :: get_pixel_layout
  public :: gen_recfac, gen_lamfac, gen_lamfac_der, gen_mfac, compute_lam_mm, gen_normpol
  !=========================================================

  !
  ! Aug 14, 2000, EH, Caltech, switch to CMBFAST normalisation of polarisation power spectra
  !  (see Zaldarriaga, astro-ph/9709271; ApJ, 503, 1 (1998))
  !   the factor normal_l (used for synthesis and analysis) is reduced by sqrt(2)
  !
  ! Oct 17, 2001
  !   polarised alm expression is larger by a factor 2
  !   introduced FL_LARGE and FL_SMALL parameters
  !
  ! Sept 2001, EH, Caltech: map2alm* : skip calculations for cut out pixels
  !                                  : use explicit loops for multiplication of alm by omega_pix
  !
  ! Nov 2001, EH, Caltech : added KIND=DP in all CMPLX commands that missed it
  !                         turned 1.d0 in 1.0_dp and so on
  !                         moved 'use healpix_types' and 'implicit none' to module top level
  !
  ! Dec 2001, EH, Caltech : corrected error on inner loop upper bound of lam_fact initialisation
  !   in alm2map_pol_pre1, alm2map_pol, map2alm_pol_pre1, map2alm_pol
  !   pointed out by M. Ashdown to Level S
  !
  !   added pow2alm_units, gaussbeam, generate_beam and alter_alm (moved for smoothing/alter_alm.f90)
  !
  !
  ! Nov 2002, EH, reordered routines to simplify maintenance of scalar vs parallel routines
  ! ------post 1.21
  ! Dec 2003-Jan 2004, EH, remove 'trig' array from ring_analysis and ring_synthesis
  ! 2004-05-28, EH, edition of pixel_window
  ! Aug 2004 : put ring in DP in ring_synthesis
  !            add alm2map_sc_der for calculation of derivatives
  !            moved plm_gen here from plmgen main code
  ! Oct 2004: reduced size of recfac and lam_fact
  !         : introduced m_max_syn and m_max_ana
  ! Nov 2004: streamlining of *alm* routines to gain CPU time
  ! Dec 2004: editions to alter_alm (new argument), pixel_window
  ! April 2005: separate Lambda_lm calculation (with do_lam_lm[_pol]) from scalar
  !  product in alm2map_sc and alm2map_pol for speed up.
  !  does not improve the speed of map2alm though, so stick to former code for those
  ! May 2005: pixel window returns 1. if nside = 0 (interpreted as infinitely small pixel)
  ! Aug 2005: added alm2map_pol_der
  ! ------post 2.00
  ! Sep-Oct 2005: made internal subroutines public so they can be called by mpi_alm, 
  !              corrected OpenMP+SGI bug in alm2map_pol_der
  ! =====================================================
  ! about the F90 compilers 'features'
  ! - Intel compiler (ifc) (and maybe the other compilers as well)
  !  is apparently very slow at allocating arrays
  ! therefore as much as possible the 'allocate' shoud be done outside the loops.
  ! This what is done below, expect when OpenMP require thread safe arrays that much
  ! be created each time we enter the loop.
  !
  ! - NagF95 (f95) is very slow at creating vectors with the operator (/ ... /)
  ! therefore, in a loop, the operation
  ! x(0:4) = x(0:4) + (/ a, b, c, d /) should be replaced by the equivalent
  ! but MUCH faster, explicit form
  ! x(0) = x(0) + a
  ! x(1) = x(1) + b
  ! ...
  ! unless (/ a,b,c,d /) can be constructed before entering the loop
  !
contains

  !**************************************************************************
  !
  !             ALM2MAP/MAP2ALM    SIDEKICKS
  !
  !**************************************************************************
  !================================================
  function do_openmp()
    !================================================
    ! returns .true. if code was compiled with OpenMP
    !================================================
    logical(LGT) :: do_openmp
    !------------------------------------------------

    do_openmp = .false.
! DO NOT REMOVE FOLLOWING LINES
!$    do_openmp = .true.  ! Intel f90
!IBMP do_openmp = .true.  ! IBM xlf90
! -----------------------------

    return
  end function do_openmp
  
  !================================================
  subroutine init_rescale()
    !================================================
    ! local variables
    integer(i4b) :: s, smax
    real(dp) :: logOVFLOW
    character(len=*), parameter :: code = 'gen_rescale'
    !------------------------------------------------
    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    if (smax > (RSMAX-1)) then
       print*,'Array rescale_tab too small in '//code
       print*,smax ,'>', RSMAX
       stop
    endif

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    return
  end subroutine init_rescale

  !=======================================================================
  subroutine get_pixel_layout(nside, ith, cth, sth, nphi, startpix, kphi0)
    !=======================================================================
    ! output Healpix pixel layout for the ring ith in [0,2*nside]
    !=======================================================================
    integer(I4B), intent(IN)  :: nside, ith
    real(DP)    , intent(OUT) :: cth, sth
    integer(I4B), intent(OUT) :: nphi, kphi0
    integer(I8B), intent(OUT) :: startpix
    !
    integer(I4B) :: nrings
    real(DP)     :: dth1, dth2, dst1
    !=======================================================================

    nrings = 2*nside
    if (ith < 1 .or. ith> nrings) then
       print*,'ith out of bounds ',ith,1,nrings
       call fatal_error
    endif

    dth1 = 1.0_dp / (3.0_dp*DBLE(nside)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nside))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nside) )

    if (ith < nside) then  ! polar cap (north)
       cth = 1.0_dp  - DBLE(ith)**2 * dth1
       nphi = 4*ith
       kphi0 = 1
       sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       startpix = 2*ith*(ith-1_I8B)
    else                   ! tropical band (north) + equator
       cth = DBLE(2*nside-ith) * dth2
       nphi = 4*nside
       kphi0 = MOD(ith+1-nside,2)
       sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       startpix = 2*nside*(nside-1_I8B) + (ith-nside)*int(nphi,kind=I8B)
    endif

    return
  end subroutine get_pixel_layout
  !=======================================================================
  subroutine select_rings(z, zbounds, keep_north, keep_south, keep_either)
    !=======================================================================
    ! select rings laying within zbounds
    ! if zbounds(1) < zbounds(2) : keep  zbounds(1) < z < zbounds(2)
    ! if zbounds(2) < zbounds(1) : keep z < zbounds(2) Union  zbounds(1) < z
    ! if zbounds(1)=zbounds(2) : keep everything
    ! input z should be >= 0
    !=======================================================================
    real(DP)    , intent(in)  :: z
    real(DP)    , intent(in), dimension(1:2)  :: zbounds
    logical(LGT), intent(OUT) :: keep_north, keep_south, keep_either
    !
    real(DP) :: zn, zs
    !=======================================================================

    ! if (zbounds(1) = zbounds(2)) keep everything
    if (abs(zbounds(1)-zbounds(2)) < 1.e-6) then
       keep_north    = .true.
       keep_south    = .true.
       keep_either   = .true.
       return
    endif

    zn = abs(z)
    zs = -zn

    if (zbounds(1) < zbounds(2)) then
       ! inner ring
       keep_north = (zn >= zbounds(1) .and. zn <= zbounds(2))
       keep_south = (zs >= zbounds(1) .and. zs <= zbounds(2))

    else
       ! outter ring
       keep_north = (zn > zbounds(1) .or. zn < zbounds(2))
       keep_south = (zs > zbounds(1) .or. zs < zbounds(2))
    endif
    keep_either   = keep_north .or. keep_south


    return
  end subroutine select_rings

  !=======================================================================
  subroutine gen_recfac( l_max, m, recfac)
  !=======================================================================
    ! generates recursion factors used to computes the Ylm of degree m 
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                            :: l_max, m
    real(DP),     intent(OUT), dimension(0:1, 0:l_max)  :: recfac
    !
    real(DP) :: fm2, fl2
    integer(I4B) :: l

    recfac(0:1,0:m) = 0.0_dp
    fm2 = DBLE(m) **2
    do l = m, l_max
       fl2 = DBLE(l+1) **2
       recfac(0,l) = DSQRT( (4.0_dp * fl2 - 1.0_dp) / (fl2-fm2) )
    enddo
    ! put outside the loop because of problem on some compilers
    recfac(1,m:l_max) = 1.0_dp / recfac(0,m:l_max)


    return
  end subroutine gen_recfac

  !=======================================================================
  subroutine gen_lamfac( l_max, m, lam_fact)
  !=======================================================================
    ! generates factor relating scalar Ylm to polar Ylm
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max, m
    real(DP),     intent(OUT), dimension(0:l_max)  :: lam_fact
    !
    real(DP) :: fm2, fl, fl2
    integer(I4B) :: l

    lam_fact(0:m) = 0.
    fm2 = real(m * m, kind=DP)
    do l = max(2,m+1), l_max
       fl  = real(l, kind=dp)
       fl2 = fl * fl
       lam_fact(l) = 2.0_dp * SQRT( (2.0_dp * fl + 1.0_dp) / (2.0_dp * fl - 1.0_dp) * (fl2-fm2) )
    enddo
    
    return
  end subroutine gen_lamfac

  !=======================================================================
  subroutine gen_lamfac_der(l_max, m, lam_fact)
    !=======================================================================
    ! generates factor relating scalar Ylm to its derivatives
    ! for all l in m<=l<=l_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max, m
    real(DP),     intent(OUT), dimension(0:l_max)  :: lam_fact
    !
    real(DP) :: fm2, fl, fl2
    integer(I4B) :: l

    lam_fact(0:m) = 0.
    fm2 = real(m * m, kind=DP)
    do l = max(1,m+1), l_max ! different lower bound than pol. factor
       fl  = real(l, kind=dp)
       fl2 = fl * fl
       lam_fact(l) = SQRT( (2.0_dp * fl + 1.0_dp) / (2.0_dp * fl - 1.0_dp) * (fl2-fm2) )
       ! different normalization than polarization factor
    enddo
    
    return
  end subroutine gen_lamfac_der

  !=======================================================================
  subroutine gen_mfac( m_max, m_fact)
  !=======================================================================
    ! generates factor used in lam_mm calculation
    ! for all m in 0<=m<=m_max
    !=======================================================================
    integer(I4B), intent(IN)                       :: m_max
    real(DP),     intent(OUT), dimension(0:m_max)  :: m_fact
    !
    integer(I4B) :: m

    ! fact(m) = fact(m-1) * sqrt( (2m+1) / (2m) )
    m_fact(0) = 1.0_dp
    do m=1,m_max
      m_fact(m) = m_fact(m-1)*sqrt(dble(2*m+1)/dble(2*m))
    end do

    ! Log_2 ( fact(m) / sqrt(4 Pi) )
    do m=0,m_max
       m_fact(m) = log(SQ4PI_INV * m_fact(m)) * ALN2_INV 
    enddo

    return
  end subroutine gen_mfac
  !=======================================================================
  subroutine compute_lam_mm(mfac, m, sth, lam_mm, corfac, scalem)
    !=======================================================================
    ! computes lam_mm
    ! the true lam_mm is     lam_mm * corfac
    !=======================================================================
    integer(I4B),            intent(in)  :: m
    real(DP),                intent(in)  :: sth, mfac
    real(DP),                intent(out) :: lam_mm, corfac
    integer(I4B),            intent(out) :: scalem
    !
    real(DP) :: log2val, dlog2lg

    dlog2lg = real(LOG2LG, kind=DP)

    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalem = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalem,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalem * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m

    return
  end subroutine compute_lam_mm
  
  !=======================================================================
  subroutine do_lam_lm(lmax, m, cth, sth, mfac, recfac, lam_lm)
    !=======================================================================
    ! computes scalar lambda_lm(theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    ! output: lam_lm
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(out) :: lam_lm
    !
    real(DP) :: log2val, dlog2lg
    real(DP) :: ovflow, unflow, corfac
    real(DP) :: lam_mm, lam_0, lam_1, lam_2
    integer(I4B) :: scalel, l, l_min
    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    ! computes lamba_mm
    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalel = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalel,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalel * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m
    
    lam_lm(m:lmax) = 0.0_dp
    ! --- l = m ---
    lam_lm(m) = lam_mm * corfac

    ! --- l > m ---
    lam_0 = 0.0_dp
    lam_1 = 1.0_dp
    lam_2 = cth * lam_1 * recfac(0,m)
    do l = m+1, lmax
       ! do recursion
       if (l >= l_min) then
          lam_lm(l) = lam_2 * corfac * lam_mm
       endif
       lam_0 = lam_1 * recfac(1,l-1)
       lam_1 = lam_2
       lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

       ! do dynamic rescaling
       if (abs(lam_2) > ovflow) then
          lam_1 = lam_1*unflow
          lam_2 = lam_2*unflow
          scalel= scalel + 1
          corfac = rescale_tab(max(scalel,RSMIN))
       elseif (abs(lam_2) < unflow) then
          lam_1 = lam_1*ovflow
          lam_2 = lam_2*ovflow
          scalel= scalel - 1
          corfac = rescale_tab(max(scalel,RSMIN))
       endif
                   
    enddo ! loop on l
  end subroutine do_lam_lm
  !=======================================================================
  subroutine do_lam_lm_pol(lmax, m, cth, sth, mfac, recfac, lam_fact, lam_lm)
    !=======================================================================
    ! computes temperature&polar lambda_lm(p,theta) for all l in [m,lmax] for a given m, and given theta
    ! input: lmax, m, cos(theta), sin(theta)
    !        mfac: precomputed (by gen_mfac) quantity useful for lambda_mm calculation
    !        recfac: precomputed (by gen_recfac) quantities useful for 
    !            lambda_lm recursion for a given m
    !        lam_fact: precomputed (by gen_lamfac) factor useful for polarised lambda recursion
    ! output: lam_lm for T and P
    ! the routine also needs the array rescale_tac initialized by init_rescale
    !=======================================================================
    integer(I4B),                    intent(in)  :: lmax,  m
    real(DP),                        intent(in)  :: cth, sth, mfac
    real(DP), dimension(0:1,0:lmax), intent(in)  :: recfac
    real(DP), dimension(    0:lmax), intent(in)  :: lam_fact
    real(DP), dimension(1:3,0:lmax), intent(out) :: lam_lm
    !
    real(DP) :: log2val, dlog2lg
    real(DP) :: ovflow, unflow, corfac
    real(DP) :: lam_mm, lam_0, lam_1, lam_2, lam_lm1m
    integer(I4B) :: scalel, l, l_min
    real(DP) :: normal_m, fm2, fl, flm1
    real(DP) :: two_cth, one_on_s2, fm_on_s2, two_on_s2, c_on_s2
    real(DP) :: a_w, a_x, b_w
    !---------------------------------------------------------------

    ! define constants
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    l_min = l_min_ylm(m, sth)
    dlog2lg = real(LOG2LG, kind=DP)
    
    fm2       = real(m * m, kind = DP)
    normal_m  = (2.0_dp * m) * (1 - m)
    two_cth   = 2.0_dp * cth
    one_on_s2 = 1.0_dp / (sth * sth)
    fm_on_s2  =      m * one_on_s2
    two_on_s2 = 2.0_dp * one_on_s2
    c_on_s2   = cth    * one_on_s2
    b_w       =  c_on_s2 
    

    ! computes lamba_mm
    log2val = mfac + m*log(sth) * ALN2_INV ! log_2(lam_mm)
    scalel = int (log2val / dlog2lg)
    corfac = rescale_tab(max(scalel,RSMIN))
    lam_mm = 2.0_dp **(log2val - scalel * dlog2lg) ! rescaled lam_mm
    if (IAND(m,1)>0) lam_mm = -lam_mm ! negative for odd m
    
    lam_lm(1:3,m:lmax) = 0.0_dp
    ! --- l = m ---
    lam_lm(1,m) = corfac * lam_mm !Actual lam_mm 

    if (m >= l_min) then ! skip Ymm if too small
       lam_lm(2,m) =  (normal_m * lam_lm(1,m))  * ( 0.5_dp - one_on_s2 )
       lam_lm(3,m) =  (normal_m * lam_lm(1,m))  *            c_on_s2
    endif

    ! --- l > m ---
    lam_0 = 0.0_dp
    lam_1 = 1.0_dp
    lam_2 = cth * lam_1 * recfac(0,m)

    do l = m+1, lmax
       ! do recursion
       lam_lm1m = lam_lm(1,l-1) * lam_fact(l) ! must be incremented even if not used
       lam_lm(1,l) = lam_2 * corfac * lam_mm
       if (l >= l_min) then
          fl = real(l, kind = DP)
          flm1 = fl - 1.0_dp
          a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
          a_x =  two_cth * flm1
          lam_lm(2,l) =                b_w * lam_lm1m - a_w * lam_lm(1,l)
          lam_lm(3,l) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm(1,l))
       endif
       lam_0 = lam_1 * recfac(1,l-1)
       lam_1 = lam_2
       lam_2 = (cth * lam_1 - lam_0) * recfac(0,l)

       ! do dynamic rescaling
       if (abs(lam_2) > ovflow) then
          lam_1 = lam_1*unflow
          lam_2 = lam_2*unflow
          scalel= scalel + 1
          corfac = rescale_tab(max(scalel,RSMIN))
       elseif (abs(lam_2) < unflow) then
          lam_1 = lam_1*ovflow
          lam_2 = lam_2*ovflow
          scalel= scalel - 1
          corfac = rescale_tab(max(scalel,RSMIN))
       endif
                   
    enddo ! loop on l
  end subroutine do_lam_lm_pol
  !=======================================================================
  subroutine gen_normpol(l_max, normal_l)
    !=======================================================================
    ! generates normalisaton factors for polarisation basis functions
    ! for all l 
    !=======================================================================
    integer(I4B), intent(IN)                       :: l_max
    real(DP),     intent(OUT), dimension(0:l_max)  :: normal_l
    !
    integer(I4B) :: l
    real(DP)     :: fl, xx

    normal_l(0:1) = 0.0_dp
    do l = 2, l_max
       fl = DBLE(l)
       xx = (fl+2.0_dp) * (fl+1.0_dp) * fl * (fl-1.0_dp)
       normal_l(l) = SQRT ( KvS / xx)
       ! either CMBFAST (KvS=1) or Kamionkowski et al. (KvS=2) definition
    enddo

    return
  end subroutine gen_normpol

  !================================================================
  function l_min_ylm(m, sth) result(lmin)
  !================================================================
    ! returns minimal order l at which to keep Ylm
    ! |Ylm| < eps * Y00 ==>
    ! m_cut(theta, l) = theta * l * e / 2 + | ln(eps)| + ln(l)/2
    ! if eps = 1.e-15 and l < 1.e4
    ! m_cut(theta, l) = theta * l * 1.35 + 40
    ! the choice of 1.35 (or larger) 
    ! also insures that the equatorial rings will have all their Ylm's computed
    ! default parameters are HPX_MXL0 = 40 and HPX_MXL1 = 1.35_DP
    !======================================================
    ! parameters of short-cut: defined in module header
    ! dummy variables
    integer(I4B)             :: lmin
    integer(I4B), intent(IN) :: m
    real(DP),     intent(IN) :: sth

    lmin = m ! default
    if (HPX_MXL0 > 0) lmin = max(lmin, int((m - HPX_MXL0)/(HPX_MXL1 * sth)))

    return
  end function l_min_ylm
  !=========================================================
  subroutine warning_oldbounds(cos_theta_cut, zbounds)
    !=========================================================
    real(DP), intent(in)                  :: cos_theta_cut
    real(DP), intent(out), dimension(1:2) :: zbounds

    if (cos_theta_cut <= 0.0_dp) then ! no cut
       zbounds(1:2) = 0.0_dp
    else
       zbounds(1) =   cos_theta_cut
       zbounds(2) = - cos_theta_cut
    endif
    print*,' -------------------------------------'
    print*,'WARNING: obsolete interface to MAP2ALM: '
    print*,'    cos_theta_cut (6th argument) currently a DP scalar with value'
    write(*,9000) '    ',cos_theta_cut
    print*,'    shoud now be replaced with a 2-element vector with values:'
    write(*,9001) '    ',zbounds(1),zbounds(2)
    print*,'    See documentation for details.'
    print*,' -------------------------------------'
9000 format (a,g12.6)
9001 format (a,g12.6,g12.6)

    return
  end subroutine warning_oldbounds

  !**************************************************************************
  !
  !             FOURIER TRANSFORM ON RINGS
  !
  !**************************************************************************
  !=======================================================================
  subroutine ring_synthesis(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
    !=======================================================================
    !     RING_SYNTHESIS
    !       called by alm2map
    !       calls     real_fft
    !
    !     dataout(j) = sum_m datain(m) * exp(i*m*phi(j))
    !     with phi(j) = j*2pi/nph + kphi0*pi/nph and kphi0 =0 or 1
    !
    !     as the set of frequencies {m} is larger than nph,
    !     we wrap frequencies within {0..nph-1}
    !     ie  m = k*nph + m' with m' in {0..nph-1}
    !     then
    !     noting bw(m') = exp(i*m'*phi0)
    !                   * sum_k (datain(k*nph+m') exp(i*k*pi*kphi0))
    !        with bw(nph-m') = CONJ(bw(m')) (if datain(-m) = CONJ(datain(m)))
    !     dataout(j) = sum_m' [ bw(m') exp (i*j*m'*2pi/nph) ]
    !                = Fourier Transform of bw
    !        is real
    !
    !         NB nph is not necessarily a power of 2
    !
    !=======================================================================


    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(IN)  :: datain
    REAL(DP),     DIMENSION(0:nph-1), INTENT(OUT) :: dataout

    INTEGER(I4B) :: iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    type(planck_fft2_plan) :: plan
    COMPLEX(DPC) :: dat
    real(DP)     :: arg
    !=======================================================================

    ksign = + 1
    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, nmmax                        ! in -nmmax, nmmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo

    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    dataout(0)=REAL(bw(0), kind=DP)
    do iw = 1, nph/2-1
       m = ksign*(iw)
       if(kphi0==1) then
          arg = m * PI / dble(nph)
          dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
       else
          dat =bw(iw)
       endif
       dataout(iw*2-1) = REAL(dat, kind=DP)
       dataout(iw*2  ) = AIMAG(dat)
    enddo
    iw=nph/2
    m = ksign*(iw)
    if(kphi0==1) then
       arg = m * PI / dble(nph)
       dat =bw(iw) * CMPLX( DCOS(arg), DSIN(arg), kind=DP)
    else
       dat =bw(iw)
    endif
    dataout(iw*2-1) = REAL(dat, kind=DP)

    call make_fft2_plan(plan,length=nph,direction=fft2_backward)
    call real_fft2 (plan, dataout)
    !     ^^^^^^^^^^^^
    call destroy_fft2_plan (plan)
    RETURN
  END subroutine ring_synthesis

  !=======================================================================
  subroutine ring_analysis(nsmax,nlmax,nmmax,datain,nph,dataout,kphi0)
    !=======================================================================
    !     ring_analysis
    !       called by map2alm
    !       calls     real_fft
    !
    !     integrates (data * phi-dependence-of-Ylm) over phi
    !     --> function of m can be computed by FFT
    !     with  0<= m <= npoints/2 (: Nyquist)
    !     because the data is real the negative m are the conjugate of the
    !     positive ones
    !=======================================================================

    INTEGER(I4B), INTENT(IN) :: nsmax, nlmax, nmmax
    INTEGER(I4B), INTENT(IN) :: nph, kphi0
    REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
    COMPLEX(DPC), DIMENSION(0:nmmax), INTENT(OUT) :: dataout

    INTEGER(I4B) :: i,m,im_max,ksign
!    REAL(DP) :: phi0
    REAL(DP), DIMENSION(0:nph-1) :: data
    type(planck_fft2_plan) :: plan
    real(DP)  :: arg

    !=======================================================================

    ksign = - 1
    data=0.
    data(0:nph-1)=datain(0:nph-1)

    call make_fft2_plan(plan,length=nph,direction=fft2_forward)
    call real_fft2(plan,data)
    call destroy_fft2_plan(plan)

    im_max = MIN(nph/2,nmmax)

    ! m = 0,  i=0
    dataout(0)=CMPLX(data(0),0.0_dp,kind=DP)

    ! 1 <= m <= im_max, --> i=1,2,3,..., im_max
    do i = 1, im_max*2-3, 2
       dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP)
    enddo

    if(im_max==nph/2) then
       dataout(im_max)= CMPLX( data(nph-1),0.0_dp,kind=DP) ! m = Nyquist
    else
       dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
    endif

    if(im_max==nmmax) goto 1000 ! m_max <= Nyquist

    ! if (m > Nyquist)
    do i =  im_max+1,min(nph,nmmax)
       dataout(i) = conjg(dataout(2*im_max-i) )
    end do

    if(min(nph,nmmax)==nmmax) goto 1000 ! nph > nmmax

    do i =  2*im_max+1,nmmax
       dataout(i) = dataout(mod(i,2*im_max))
    end do

1000 continue

    if(kphi0==1)then
       do i =  0,nmmax
          m = ksign*i
!           dataout(i)=dataout(i)* CONJG(trig(-m,nph/4))
          arg = m * PI / dble(nph)
          dataout(i)=dataout(i)* CMPLX( DCOS(arg), DSIN(arg), kind=DP)
       enddo
    end if

    RETURN
  END subroutine ring_analysis

  !**************************************************************************
  !   ALM2MAP
  !   MAP2ALM
  !   ALM creation and alteration
  !   import overloaded routines
  !**************************************************************************

  ! single precision routines
  include 'alm_map_ss_inc.f90'

  ! double precision routines
  include 'alm_map_dd_inc.f90'
  

  !**************************************************************************
  !
  !             PLM GENERATION
  !
  !**************************************************************************
  !========================================================
  subroutine plm_gen(nsmax, nlmax, nmmax, plm)
    !========================================================
    integer(i4b),             intent(IN) :: nsmax, nlmax, nmmax
    real(dp), dimension(0:,1:), intent(OUT):: plm

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nd2, nrings
    integer(I8B) :: nd1, n_lm, n_plm, i_mm, i_up
    real(DP)            :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
!!    real(DP)                 :: lambda_w, lambda_x
    real(DP)                 :: normal_m, lam_lm1m
    real(DP)                 :: fm2, fl, flm1, fpol
    real(DP)                 :: a_w, b_w, a_x, fm_on_s2, two_on_s2, two_cth_ring
    real(DP)                 :: ovflow, unflow
    real(DP),     dimension(:,:,:), allocatable :: plm_sub
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP),     dimension(:),   allocatable :: lam_fact
    real(DP),     dimension(:), allocatable :: mfac

    real(DP),     dimension(:),     allocatable :: normal_l
    integer(i4b)        :: nchunks, chunksize, ichunk, lchk, uchk, ithl
    integer(i4b)        :: nph, kphi0, i
    integer(i8b)        :: startpix
    real(DP),     dimension(0:SMAXCHK) :: cth, sth
    real(DP),     dimension(0:SMAXCHK) :: one_on_s2, c_on_s2

    INTEGER(I4B) :: status
    LOGICAL(LGT) :: polarisation
    character(len=*), parameter :: code = 'PLM_GEN'

    !=================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat

    n_lm  = ((nmmax+1)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm = n_lm * nrings
    nd1 = size(plm, 1)
    nd2 = size(plm, 2)

    if (nd1 < n_plm) then
       print*,code//' > Plm array too small:', nd1, n_plm
       stop
    endif
    if (nd2 /= 1 .and. nd2 /= 3) then
       print*,code//' > Plm array should have dimension 1 or 3:',nd2
       stop
    endif
    polarisation = (nd2 == 3)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')
    if (polarisation) then
       allocate(normal_l(0:nlmax),stat = status)
       call assert_alloc(status,code,'normal_l')
    endif

    if (.not. do_openmp()) then
       allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
       call assert_alloc(status,code,'recfac & plm_sub')
       if (polarisation) then
          allocate(lam_fact(0:nlmax),stat = status)
          call assert_alloc(status,code,'lam_fact')
       endif
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax, mfac)
    ! generate Polarization normalisation
    if (polarisation) call gen_normpol(nlmax, normal_l)
    call init_rescale()
    ovflow = rescale_tab(1)
    unflow = rescale_tab(-1)
    plm = 0.0_dp

    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph, startpix, kphi0)
          one_on_s2(ithl) =    1.0_dp / sth(ithl)**2
            c_on_s2(ithl) = cth(ithl) / sth(ithl)**2
       enddo

!$OMP parallel default(none) &
!$OMP shared(nlmax, nmmax, lchk, uchk, nd2, chunksize, n_lm, &
!$OMP    rescale_tab, ovflow, unflow, polarisation, &
!$OMP    plm, cth, sth, mfac, normal_l, one_on_s2, c_on_s2) &
!$OMP private(plm_sub, recfac, lam_fact, status, &
!$OMP   m, fm2, normal_m, ithl, ith, i_mm, i_up, &
!$OMP   scalem, scalel, corfac, fpol, &
!$OMP   lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2, &
!$OMP   cth_ring, fm_on_s2, two_on_s2, two_cth_ring, a_w, a_x, b_w,  &
!$OMP   l, fl, flm1, i)

       if (do_openmp()) then
          ! allocate thread safe arrays
          allocate(recfac(0:1,0:nlmax), plm_sub(1:nd2,0:nlmax,0:chunksize-1), stat = status)
          call assert_alloc(status,code,'recfac & plm_sub')
          if (polarisation) then
             allocate(lam_fact(0:nlmax),stat = status)
             call assert_alloc(status,code,'lam_fact')
          endif
       endif

!$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)
          if (polarisation) then
             ! generate Ylm relation factor for degree m
             call gen_lamfac(nlmax, m, lam_fact)
             fm2 = real(m * m, kind = DP)
             normal_m = (2.0_dp * m) * (1 - m)
          endif

          do ithl = 0, uchk - lchk
             ! determine lam_mm
             call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)
             ! ---------- l = m ----------
             !           temperature 
             lam_lm = corfac * lam_mm !Actual lam_mm 
             plm_sub(1, m, ithl) = lam_lm

             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             cth_ring = cth(ithl)
             lam_2 = cth_ring * lam_1 * recfac(0,m)

             if (polarisation) then
                fpol = normal_m * normal_l(m) * lam_lm
                plm_sub(2, m, ithl) =  fpol * ( 0.5_dp - one_on_s2(ithl) )
                plm_sub(3, m, ithl) =  fpol *              c_on_s2(ithl)
                !
                fm_on_s2     =      m * one_on_s2(ithl)
                two_on_s2    = 2.0_dp * one_on_s2(ithl)
                two_cth_ring = 2.0_dp * cth_ring
                b_w          =  c_on_s2(ithl) 
             endif
             ! ---------- l > m ----------
             do l = m+1, nlmax
                if (polarisation) lam_lm1m = lam_lm * lam_fact(l)
                lam_lm = lam_2 * corfac * lam_mm
                plm_sub(1, l, ithl) = lam_lm
                
                if (polarisation) then
                   fl = real(l, kind = DP)
                   flm1 = fl - 1.0_dp
                   a_w =  two_on_s2 * (fl - fm2)  + flm1 * fl
                   a_x =  two_cth_ring * flm1
                   plm_sub(2, l, ithl) =            (   b_w * lam_lm1m - a_w * lam_lm) * normal_l(l)
                   plm_sub(3, l, ithl) = fm_on_s2 * (         lam_lm1m - a_x * lam_lm) * normal_l(l)
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
          enddo ! loop on rings (ithl)

          ! do memory skipping operations outside inner loops
          do ith = lchk, uchk
             i_mm = n_lm * (ith-1) + ((2*nlmax+3-m)*m)/2 ! location of Ym,m for ring ith
             i_up = i_mm + nlmax - m ! location of Ynlmax,m for ring ith
             ithl = ith - lchk
             do i=1, nd2
                plm(i_mm:i_up, i) = plm_sub(i, m:nlmax, ithl)
             enddo
          enddo
!!!!!!!!
!           i_mm = n_lm*chunksize*ichunk + ((2*nlmax+3-m)*m)/2
!           do ithl = 0, uchk-lchk
!              i_up = i_mm+(nlmax-m+1)*ithl
!           do i=1,nd2
!              plm(i_up:i_up+nlmax-m,i) = plm_sub(i, m:nlmax, ithl)
!           enddo
!           enddo
          !------------------------------------------------
       enddo ! loop on m
!$OMP end do
       if (do_openmp()) then
          deallocate (recfac, plm_sub)
          if (polarisation) deallocate(lam_fact)
       endif
!$OMP end parallel


    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not. do_openmp()) then
       deallocate (recfac, plm_sub)
       if (polarisation) deallocate(lam_fact)
    endif
    deallocate(mfac)
    if (polarisation) deallocate(normal_l)

    return
  end subroutine plm_gen
  !**************************************************************************
  !
  !             M I S C
  !
  !**************************************************************************
  !=======================================================================
  subroutine pow2alm_units(units_pow, units_alm)
    !=======================================================================
    ! does the unit conversion between power and alm (and map)
    !=======================================================================
    character(len=*), dimension(1:), intent(in)  :: units_pow
    character(len=*), dimension(1:), intent(out) :: units_alm

    integer(kind=i4b) :: i, nu, j, ntemplate, idx
    character(len=80) :: uu, ui, uo
    character(len=5),dimension(1:7) :: template = &
         & (/ "_SQUA","-SQUA","SQUA ","^2   ","^ 2  ","**2  ","** 2 " /)
    !=======================================================================
    nu = min( size(units_pow), size(units_alm))
    ntemplate = size(template)

    units_alm = ""
    do i = 1, nu
       ui = adjustl(units_pow(i))
       uu = trim(strupcase(ui))
       uo = "unknown" ! default
       do j = 1, ntemplate
          idx = index(uu, template(j))
          if (idx > 0) then
             uo = ui(1:idx-1)
             exit
          endif
       enddo
       units_alm(i) = uo
       !         print*,i,trim(units_pow(i))," -> ",trim(units_alm(i))
    enddo

    return
  end subroutine pow2alm_units
  !****************************************************************************
  subroutine generate_beam(fwhm_arcmin, lmax, gb, beam_file)
    !==========================================================================
    !
    ! generate_beam(fwhm_arcmin, lmax, gb [, beam_file])
    ! generate the beam window function gb up to lmax
    !  either a gaussian of FHWM : fwhm_arcmin (in arcmin)
    !  or the function read from beam_file
    !
    ! July 2003, EH : replicate temperature beam for polarization when reading
    !  standard ASCII file
    !==========================================================================
!     use fitstools, only : getsize_fits, fits2cl
    use fitstools, only : fits2cl
    real(kind=DP),                   intent(in)  :: fwhm_arcmin
    real(kind=DP), dimension(0:,1:), intent(out) :: gb
    integer(kind=I4B),               intent(in)  :: lmax
    character(len=*),      optional, intent(in)  :: beam_file

    character(len=256) :: new_beam_file
    logical(kind=lgt) :: extfile
    integer(kind=i4b) :: type, nl, nd, lunit, il, i
    character(len=80), dimension(1:180) :: header
    character(len=1600) :: str
    character(len=80) :: card
    !==========================================================================
    ! test if name of external is given and valid
    extfile = .false.
    if (present(beam_file)) then
       extfile = (trim(beam_file) /= "")
    endif

    lunit = 15
    nl = size(gb, 1)
    nd = size(gb, 2)
    gb = 0.0_dp

    if (nl <= lmax) then
       print*,'Generate_beam: beam array only available up to ',nl
    endif
    nl = min(nl, lmax+1)

    if (extfile) then
       call assert_present(beam_file)
       new_beam_file = beam_file
       print*,'Reading beam information from '//trim(new_beam_file)

       ! find out file nature
       type = 1
       open(unit=lunit,file=new_beam_file,status='old', &
               &          form='formatted',action='read')
       read(lunit,'(a)') card
       close(lunit)
       card = adjustl(card)
       if (card(1:8) /= 'SIMPLE  ' .AND. card(1:8) /= 'XTENSION') type = -1

       ! read file according to its type
       if (type < 0) then 
          ! ordinary ascii file ?
          lunit = 32
          open(unit=lunit,file=new_beam_file,status='old', &
               &          form='formatted',action='read')
          do
             read(lunit,'(a)', end=100, err=100) str
             if (str(1:1) /= '#') then
                read(str,*) il, gb(il,1)
                if (il == nl-1) exit
             endif
          enddo
100       continue
          close(lunit)
          if (il < (nl-1)) then
             print*,'WARNING: Beam transfer function only available up to l= ',il !
             print*,'         The larger multipoles will be set to 0'
          endif

       else if (type == 1) then
          ! FITS file with ascii table
          call fits2cl(new_beam_file, gb, nl-1, nd, header)
       else
          print*,'the file '//trim(new_beam_file) &
               &            //' is of unknown type.'
          call fatal_error
       endif

       ! if Grad absent, replicate Temperature; if Curl absent, replicate Grad
       do i=2,nd
          if ( sum(abs(gb(:,i))) < 1.e-7 ) then
             print*,'column #',i,' empty, fill in with column #',i-1
             gb(:,i) = gb(:,i-1)
          endif
       enddo

    else
       ! gaussian beam
       print*,'Generating gaussian beam of FHWM [arcmin] = ',fwhm_arcmin
       call gaussbeam(fwhm_arcmin, nl-1, gb)
    endif

    return
  end subroutine generate_beam
  !*************************************************************
  subroutine gaussbeam(fwhm_arcmin, lmax, gb)
    !===========================================================
    ! gaussbeam(fwhm_arcmin, gb)
    !   returns in gb the window function on [0,lmax] corresponding
    !   to the gaussian beam of FWHM = fwhm_arcmin
    ! The returned beam function has up to 3 components,
    ! gb(*,1) = bt                  : temperature
    ! gb(*,2) = bt * exp(2*sigma^2) : grad
    ! gb(*,3) = bt * exp(2*sigma^2) : curl
    ! with sigma = gaussian rms in radian
    ! and bt = exp( l(l+1) * sigma^2 / 2)
    !===========================================================
    real(kind=DP),                   intent(in)  :: fwhm_arcmin
    real(kind=DP), dimension(0:,1:), intent(out) :: gb
    integer(kind=I4B),               intent(in)  :: lmax

    integer(kind=I4B) :: l, nd
    real(kind=DP)     :: sigma, arcmin2rad, sigma2fwhm, fact_pol
    !===========================================================

    call assert (fwhm_arcmin>=0.0_dp,'fwhm of gaussian beam should be positive')

    nd   = size(gb,2)

    arcmin2rad = PI / (180.0_dp * 60.0_dp)
    sigma2fwhm = sqrt(8.0_dp * log(2.0_dp))

    sigma    = fwhm_arcmin * arcmin2rad / sigma2fwhm ! in radians

    fact_pol = exp(2.0_dp*sigma**2) ! correction for polarised fields

    ! temperature
    do l=0,lmax
       gb(l,1) = exp(-.5_dp * l*(l+1.0_dp) * sigma**2)
    enddo
    ! electric or gradient
    if (nd > 1) gb(0:lmax,2) = gb(0:lmax,1) * fact_pol
    ! magnetic or curl
    if (nd > 2) gb(0:lmax,3) = gb(0:lmax,1) * fact_pol

    return
  end subroutine gaussbeam

  !=======================================================================
  subroutine pixel_window(pixlw, nside, windowfile)
    !=======================================================================
    !=======================================================================

    use fitstools, only : getsize_fits, read_dbintab
    use extension, only : getEnvironment
    use pix_tools, only : nside2npix

    real(DP),     dimension(0:,1:), intent(OUT)          :: pixlw
    character(len=*),             intent(IN), optional :: windowfile
    integer(i4b),                 intent(IN), optional :: nside

    character(len=filenamelen) :: wfile, wfile_def, healpixdir
    character(len=4) :: sstr
    integer(I4B) :: npw_file, ncolw_file, npix
    integer(I4B) :: npw, ncolw, i
    logical(LGT) :: good, bad
    REAL(DP) ::  nullval
    character(len=*), parameter :: code = 'Pixel_Window > '
    real(DP),     dimension(:,:), allocatable :: pixtmp
    !=======================================================================

    call assert (present(nside) .or. present(windowfile), &
      code//'Nside or windowfile have to be present')

    npw   = size(pixlw, dim=1)
    ncolw = size(pixlw, dim=2)

    if (present(windowfile)) then
       wfile = trim(windowfile)
       inquire(file=wfile, exist = good)
    else
       if (nside == 0) then
          pixlw = 1.0_dp
          return
       endif
       npix = nside2npix(nside)
       if (npix < 0) then
          print*,code//'invalid Nside = ',nside
          call fatal_error
       endif
       sstr = trim(string(nside,'(i4.4)'))
       wfile_def = "data/pixel_window_n"//trim(sstr)//".fits"
       call getEnvironment("HEALPIX",healpixdir)

       wfile = trim(healpixdir)//wfile_def
       inquire(file=wfile, exist = good)
       if (good) goto 10
       wfile = trim(healpixdir)//'/'//wfile_def
       inquire(file=wfile, exist = good)
       if (good) goto 10
    endif
    if (.not.good) then
       print*,code//trim(wfile)//' not found'
       call fatal_error
    endif
10  continue

    npw_file = getsize_fits(wfile, nmaps = ncolw_file)
    if (ncolw_file > 3) then
       print*,code//' Too many columns in '//trim(wfile)
       call fatal_error
    endif
    ! allocate tmp array large enough to read the whole file
    allocate(pixtmp(0:npw_file-1, 1:ncolw_file))

    if (npw_file < npw) then
       print*,'Only ',npw_file,' multipoles available in '//trim(wfile)//', expected ',npw
    endif

!     call read_dbintab(windowfile, pixlw, npw, ncolw_file, nullval, bad)
!     edited 2004-05-28
!     call read_dbintab(wfile, pixlw, npw, ncolw_file, nullval, bad)

    call read_dbintab(wfile, pixtmp, npw_file, ncolw_file, nullval, bad)

    npw = min(npw_file, npw)
    ncolw_file = min(ncolw_file, ncolw)

    do i=1, ncolw_file
       pixlw(0:npw-1, i) = pixtmp(0:npw-1, i)
    enddo
    if (ncolw_file == 1) then       ! T only in file
       if (ncolw > 1) pixlw(:,2)  = pixlw(:,1) ! w_G = w_T
       if (ncolw > 2) pixlw(:,3)  = pixlw(:,1) ! w_C = w_T
    else if (ncolw_file == 2) then  ! T and G in file
       if (ncolw > 2) pixlw(:,3)  = pixlw(:,2) ! w_C = w_G
    endif

    deallocate(pixtmp)

    return
  end subroutine pixel_window


!***********************************************************************************

end module alm_tools
