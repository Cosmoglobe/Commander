!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_tod_WMAP_mod
  !   Module which contains all the WMAP time ordered data processing and routines
  !   for a given frequency band
  !
  !   Main Methods
  !   ------------
  !   constructor(cpar, id_abs, info, tod_type)
  !       Initialization routine that reads in, allocates and associates 
  !       all data needed for TOD processing
  !   process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out)
  !       Routine which processes the time ordered data
   use comm_tod_driver_mod
   use comm_conviqt_mod
   implicit none

   private
   public comm_WMAP_tod

   type, extends(comm_tod) :: comm_WMAP_tod
      integer(i4b) :: nside_M_lowres, nmaps_M_lowres
      logical(lgt) :: comp_S
      character(len=20), allocatable, dimension(:) :: labels ! names of fields
      real(dp), allocatable, dimension(:,:)        :: M_lowres, M_diag
      character(len=512) :: noise_format
   contains
      procedure     :: process_tod             => process_WMAP_tod
      procedure     :: precompute_M_lowres
      procedure     :: apply_map_precond       => apply_wmap_precond
      procedure     :: construct_corrtemp_inst => construct_corrtemp_wmap
   end type comm_WMAP_tod

   interface comm_WMAP_tod
      procedure constructor_wmap
   end interface comm_WMAP_tod

contains



   !**************************************************
   !             Constructor
   !**************************************************
   function constructor_wmap(cpar, id, id_abs, info, tod_type) result(c)
      !
      ! Constructor function that gathers all the instrument parameters in a pointer
      ! and constructs the objects
      !
      ! Arguments:
      ! ----------
      ! cpar:     derived type
      !           Object containing parameters from the parameterfile.
      ! id_abs:   integer
      !           The index of the current band within the parameters, related to cpar
      ! info:     map_info structure
      !           Information about the maps for this band, e.g.,
      !           how the maps are distributed in memory
      ! tod_type: string
      !           Instrument specific tod type
      !
      ! Returns
      ! ----------
      ! constructor: pointer
      !              Pointer that contains all instrument data

      implicit none
      type(comm_params),      intent(in) :: cpar
      integer(i4b),           intent(in) :: id, id_abs
      class(comm_mapinfo),    target     :: info
      character(len=128),     intent(in) :: tod_type
      class(comm_WMAP_tod),   pointer    :: c

      integer(i4b) :: i, nside_beam, lmax_beam, nmaps_beam
      logical(lgt) :: pol_beam

      real(dp),  allocatable, dimension(:, :)      :: m_buf




      ! Initialize common parameters
      allocate (c)


    end function constructor_wmap

   !**************************************************
   !             Driver routine
   !**************************************************
   subroutine process_WMAP_tod(self, chaindir, chain, iter, handle, map_in, delta, map_out, rms_out, map_gain)
      !
      ! Routine that processes the WMAP time ordered data.
      ! Samples absolute and relative bandpass, gain and correlated noise in time domain,
      ! perform data selection, correct for sidelobes, compute chisquare  and outputs maps and rms.
      ! Writes maps to disc in fits format
      !
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_WMAP_tod class
      !           Points to output of the constructor
      ! chaindir: string
      !           Directory for output files
      ! chain:    integer
      !           Index number of the chain being run
      ! iter:     integer
      !           Gibbs iteration number
      ! handle:   planck_rng derived type
      !           Healpix definition for random number generation
      !           so that the same sequence can be resumed later on from that same point
      ! map_in:   array
      !           Array of dimension (ndet,ndelta) with pointer to maps,
      !           with both access to maps and changing them.
      !           ndet is the number of detectors and
      !           ndelta is the number of bandpass deltas being considered
      ! delta:    array
      !           Array of bandpass corrections with dimensions (0:ndet,npar,ndelta)
      !           where ndet is number of detectors, npar is number of parameters
      !           and ndelta is the number of bandpass deltas being considered
      !
      ! Returns:
      ! ----------
      ! map_out: comm_map class
      !          Final output map after TOD processing combined for all detectors
      ! rms_out: comm_map class
      !          Final output rms map after TOD processing combined for all detectors
      implicit none
      class(comm_WMAP_tod),             intent(inout) :: self
      character(len=*),                    intent(in) :: chaindir
      integer(i4b),                        intent(in) :: chain, iter
      type(planck_rng),                 intent(inout) :: handle
      type(map_ptr), dimension(1:, 1:), intent(inout) :: map_in    ! (ndet,ndelta)
      real(dp),  dimension(0:, 1:, 1:), intent(inout) :: delta     ! (0:ndet,npar,ndelta) BP corrections
      class(comm_map), intent(inout) :: map_out      ! Combined output map
      class(comm_map), intent(inout) :: rms_out      ! Combined output rms
      type(map_ptr), dimension(1:, 1:), intent(inout), optional :: map_gain    ! (ndet,1)

      real(dp)     :: t1, t2, monopole, sigma_mono
      integer(i4b) :: i, j, k, l, n
      integer(i4b) :: nside, npix, nmaps 
      integer(i4b) :: ierr, ndelta, t_mid=53765
      real(sp), allocatable, dimension(:, :)          :: s_buf
      real(sp), allocatable, dimension(:, :, :)       :: d_calib
      real(dp), allocatable, dimension(:, :)          :: chisq_S, m_buf
      real(dp), allocatable, dimension(:, :)          :: M_diag, buffer1
      real(dp), allocatable, dimension(:, :, :)       :: M_diag_1
      real(dp), allocatable, dimension(:)             :: II_inv, QQ_inv, UU_inv, QU_inv, det
      real(dp), allocatable, dimension(:, :, :)       :: b_map, b_mono, sys_mono, buffer2
      real(dp), allocatable, dimension(:,:,:,:)       :: b_map_1, b_map_2
      character(len=512) :: prefix, postfix
      character(len=2048) :: Sfilename

      logical(lgt)        :: select_data, sample_abs_bandpass, sample_rel_bandpass, bp_corr, output_scanlist, split
      type(comm_scandata) :: sd

      character(len=4)   :: ctext, myid_text
      character(len=6)   :: samptext, scantext
      character(len=512), allocatable, dimension(:) :: slist
      real(sp),       allocatable, dimension(:)     :: procmask, procmask2, sigma0
      real(sp),  allocatable, dimension(:, :, :, :) :: map_sky, m_gain
      class(map_ptr),     allocatable, dimension(:) :: outmaps

      ! biconjugate gradient-stab parameters
      integer(i4b) :: num_cg_iters
      real(dp) ::  epsil
      real(dp) ::  nullval
      real(dp), allocatable, dimension(:, :)    :: bicg_sol, bicg_sol_1, bicg_sol_2
      real(dp), allocatable, dimension(:, :)    :: map_full
      class(comm_map), pointer :: wmap_guess

      ! Counting data kept, lost
      integer(i8b) :: n_tot, n_flag, n_discard

      character(len=80), dimension(180) :: header



   end subroutine process_WMAP_tod


   subroutine precompute_M_lowres(self)
      !
      ! Routine that precomputes the low-resolution preconditioner, M_lowres = (P^t invN_w P)^{-1}
      !
      ! Arguments:
      ! ----------
      ! self:     pointer of comm_WMAP_tod class
      !           Points to output of the constructor
      !
      implicit none
      class(comm_WMAP_tod),             intent(inout) :: self

      integer(i4b) :: i, j, k, t, p1, p2, p1_l,p1_r,p2_l,p2_r,k1, k2, ntot, npix, npix_hi, ierr, ntod, lpix, rpix, q, nhorn, lpsi, rpsi
      real(dp)     :: var, inv_sigma, lcos2psi, lsin2psi, rcos2psi, rsin2psi
      real(dp)     :: dx, xbar, f_l, f_r, mA, mB
      real(dp), allocatable, dimension(:)   :: dl, dr, pl, pr
      real(dp), allocatable, dimension(:,:) :: M
      integer(i4b), allocatable, dimension(:)         :: flag, dgrade
      real(sp),  allocatable, dimension(:)         :: procmask
      real(dp),  allocatable, dimension(:, :)      :: m_buf
      integer(i4b), allocatable, dimension(:, :)      :: pix, psi
      type(hdf_file) :: precond_file




    end subroutine precompute_M_lowres


  subroutine apply_wmap_precond(self, map, map_out)
    implicit none
    class(comm_WMAP_tod),              intent(in)    :: self
    real(dp),        dimension(0:,1:), intent(in)    :: map
    real(dp),        dimension(0:,1:), intent(out)   :: map_out

    integer(i4b) :: i, npix_lowres, n_lowres, nmaps, n, p, q, j
    real(dp) :: determ
    real(dp), allocatable, dimension(:)   :: m_lin, m
    real(dp), allocatable, dimension(:,:) :: m_low
    !
    !   Routine follows Section 3.4.7 of Jarosik et al. 2007
    !
    !   Arguments: 
    !   ----------
    !   self:     pointer of comm_WMAP_tod class
    !             Points to output of the constructor
    !   map:      Map to be preconditioned
    !
    !   map_out:  Map after being preconditioned
    !
    !   Intermediate steps:
    !   -------------------
    !   m_lin - a linearized map, length nmaps * npix

    if (self%comp_S) then
       map_out =  map/self%M_diag
    else

       map_out = 0d0

       npix_lowres = 12*self%nside_M_lowres**2
       nmaps       = self%nmaps_M_lowres

       ! Apply lowres preconditioner
       allocate(m_lin(0:npix_lowres*nmaps-1), m(0:size(map,1)-1))
       allocate(m_low(0:size(map,1)-1, nmaps))
       do i = 1, nmaps
          m = map(:,i)
          call udgrade_ring(m, self%info%nside, m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres)
       end do
       ! m_lin is now the low resolution linearized version of the map
       m_lin = matmul(self%M_lowres, m_lin)
       ! m_lin has now been preconditioned.

       do i = 1, nmaps
          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, m_low(:,i), self%info%nside)
       end do
       
       ! Apply highres preconditioner to residual
       map_out = map - m_low
       do i = 0, size(map,1)-1
          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
          map_out(i,1) =  map_out(i,1)/self%M_diag(i,1)
          map_out(i,2) = (map_out(i,2)*self%M_diag(i,3) - map_out(i,2)*self%M_diag(i,4))/determ
          map_out(i,3) = (map_out(i,3)*self%M_diag(i,2) - map_out(i,3)*self%M_diag(i,4))/determ
       end do

       do i = 1, nmaps
          call udgrade_ring(m_lin((i-1)*npix_lowres:i*npix_lowres-1), self%nside_M_lowres, m, self%info%nside)
          map_out(:,i) = map_out(:,i) + m_low(:,i)
       end do

       deallocate(m, m_lin, m_low)
       

!       do i = 0, size(map,1)-1
!          determ       = self%M_diag(i,2)*self%M_diag(i,3) - self%M_diag(i,4)**2
!          map_out(i,1) =  map(i,1)/self%M_diag(i,1)
!          map_out(i,2) = (map(i,2)*self%M_diag(i,3) - map(i,2)*self%M_diag(i,4))/determ
!          map_out(i,3) = (map(i,3)*self%M_diag(i,2) - map(i,3)*self%M_diag(i,4))/determ
!       end do


    end if

  end subroutine apply_wmap_precond

  subroutine sample_baseline_WMAP(tod, scan, raw, s_tot, mask, handle)
    !   Sample LFI specific 1Hz spikes shapes and amplitudes
    !
    !   Arguments:
    !   ----------
    !   tod:      comm_tod derived type
    !             contains TOD-specific information
    !   scan:     local scan ID
    !   raw:      raw tod in du
    !   s_tot:    total signal model in mK
    !   mask:     list of accepted samples
    !   handle:   planck_rng derived type
    !             Healpix definition for random number generation
    implicit none
    class(comm_wmap_tod),                   intent(inout) :: tod
    integer(i4b),                           intent(in)    :: scan
    real(sp),            dimension(1:,1:),  intent(in)    :: raw, s_tot, mask
    type(planck_rng),                       intent(inout) :: handle

    integer(i4b) :: i, j, k, n
    real(dp)     :: dt, t_tot, t, A, b, mval, eta
    real(dp), allocatable, dimension(:) :: x, y 

    allocate(x(tod%scans(scan)%ntod), y(tod%scans(scan)%ntod))
    dt = 1.d0 / tod%scans(scan)%ntod


    do j = 1, tod%ndet
       if (.not. tod%scans(scan)%d(j)%accept) cycle

       t = 0.d0
       n = 0
       do k = 1, tod%scans(scan)%ntod
          t      = t + dt
          if (mask(k,j) > 0.5) then
             n    = n + 1
             x(n) = t
             y(n) = raw(k,j) - tod%scans(scan)%d(j)%gain * s_tot(k,j)
          end if
       end do

       call fit_polynomial(x(1:n), y(1:n), tod%scans(scan)%d(j)%baseline)
    end do

    deallocate(x, y)

  end subroutine sample_baseline_WMAP

  subroutine construct_corrtemp_wmap(self, scan, pix, psi, s)
    !  Construct an WMAP instrument-specific correction template; for now contains baseline
    !
    !  Arguments:
    !  ----------
    !  self: comm_tod object
    !
    !  scan: int
    !       scan number
    !  pix: int
    !       index for pixel
    !  psi: int
    !       integer label for polarization angle
    !
    !  Returns:
    !  --------
    !  s:   real (sp)
    !       output template timestream
    implicit none
    class(comm_wmap_tod),                  intent(in)    :: self
    integer(i4b),                          intent(in)    :: scan
    integer(i4b),        dimension(:,:),   intent(in)    :: pix, psi
    real(sp),            dimension(:,:),   intent(out)   :: s

    integer(i4b) :: i, j, k, nbin, b
    real(dp)     :: dt, t


    dt = 1.d0 / self%scans(scan)%ntod
    s = 0.
    do j = 1, self%ndet
       t = 0.d0
       if (.not. self%scans(scan)%d(j)%accept) cycle
       do k = 1, self%scans(scan)%ntod
          t      = t + dt
          do i = 0, self%baseline_order
             s(k,j) = s(k,j) + self%scans(scan)%d(j)%baseline(i) * t**i
          end do
       end do
    end do


  end subroutine construct_corrtemp_wmap


end module comm_tod_WMAP_mod
