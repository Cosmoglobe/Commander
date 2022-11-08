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
module comm_N_mod
  use comm_param_mod
  use comm_map_mod
  use comm_status_mod
  implicit none

  private
  public comm_N, compute_invN_lm, uniformize_rms
  
  type, abstract :: comm_N
     ! Data variables
     character(len=512)       :: type
     integer(i4b)             :: nside, nside_chisq_lowres, nmaps, np, npix, myid, comm, nprocs
     logical(lgt)             :: pol, pol_only
     logical(lgt)             :: set_noise_to_mean
     character(len=512)       :: cg_precond
     real(dp)                 :: uni_fsky
     real(dp), allocatable, dimension(:) :: alpha_nu ! (T,Q,U)
     class(comm_map),     pointer :: invN_diag => null()
     class(comm_map),     pointer :: rms_reg   => null()
     class(comm_mapinfo), pointer :: info      => null()
     class(map_ptr), allocatable, dimension(:) :: samp_group_mask
   contains
     ! Data procedures
     procedure(matmulInvN),       deferred :: invN
     procedure(matmulInvNlowres), deferred :: invN_lowres
     procedure(matmulN),          deferred :: N
     procedure(matmulSqrtInvN),   deferred :: sqrtInvN
     procedure(returnRMS),        deferred :: rms
     procedure(returnRMSpix),     deferred :: rms_pix
     procedure(update_N),         deferred :: update_N
     procedure                             :: P => apply_projection
  end type comm_N
  
  abstract interface
     ! Return map_out = invN * map
     subroutine matmulInvN(self, map, samp_group)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulInvN

     ! Return map_out = invN * map
     subroutine matmulInvNlowres(self, map, samp_group)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulInvNlowres

     ! Return map_out = N * map
     subroutine matmulN(self, map, samp_group)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulN

     ! Return map_out = sqrtInvN * map
     subroutine matmulSqrtInvN(self, map, samp_group)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulSqrtInvN

     ! Return rms map
     subroutine returnRMS(self, res, samp_group)
       import comm_map, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: res
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine returnRMS

     ! Return rms map
     function returnRMSpix(self, pix, pol, samp_group)
       import i4b, comm_N, dp, i4b
       implicit none
       class(comm_N),   intent(in)             :: self
       integer(i4b),    intent(in)             :: pix, pol
       real(dp)                                :: returnRMSpix
       integer(i4b),    intent(in),   optional :: samp_group
     end function returnRMSpix

     ! Update noise model
     subroutine update_N(self, info, handle, mask, regnoise, procmask, noisefile, map)
       import comm_N, comm_mapinfo, comm_map, dp, planck_rng
       implicit none
       class(comm_N),                      intent(inout)          :: self
       type(planck_rng),                   intent(inout)          :: handle
       class(comm_mapinfo),                intent(in)             :: info
       class(comm_map),                    intent(in),   optional :: mask
       real(dp),         dimension(0:,1:), intent(out),  optional :: regnoise
       class(comm_map),                    intent(in),   optional :: procmask
       character(len=*),                   intent(in),   optional :: noisefile
       class(comm_map),                    intent(in),   optional :: map
     end subroutine update_N

  end interface

contains

  subroutine compute_invN_lm(invN_diag)
    implicit none

    class(comm_map),  intent(inout) :: invN_diag

    real(dp)     :: l0min, l0max, l1min, l1max, npix, t1, t2
    integer(i4b) :: j, l, m, lp, l2, ier, twolmaxp2, pos, lmax
    complex(dpc) :: val(invN_diag%info%nmaps)
    real(dp), allocatable, dimension(:)   :: threej_symbols, threej_symbols_m0
    real(dp), allocatable, dimension(:,:) :: N_lm, a_l0

    lmax      = invN_diag%info%lmax
    twolmaxp2 = 2*lmax+2
    npix      = real(invN_diag%info%npix,dp)

    allocate(N_lm(0:invN_diag%info%nalm-1,invN_diag%info%nmaps))
    allocate(a_l0(0:lmax,invN_diag%info%nmaps))
    call invN_diag%YtW_scalar

    if (invN_diag%info%myid == 0) then
       ! Seriously ugly hack. Should figure out how to compute N_{lm,l'm'} properly
       ! for polarization
       a_l0(:,:) = invN_diag%alm(0:lmax,:)
    end if
    call mpi_bcast(a_l0, size(a_l0), MPI_DOUBLE_PRECISION, 0, invN_diag%info%comm, ier)

    call update_status(status, "compute_invN_lm 3j")

    !call wall_time(t1)
    !!$OMP PARALLEL PRIVATE(pos,j,m,l,threej_symbols_m0,threej_symbols,ier,val,l2,lp,l0min,l0max,l1min,l1max)
    allocate(threej_symbols(twolmaxp2))
    allocate(threej_symbols_m0(twolmaxp2))
    !!$OMP DO SCHEDULE(guided)
    do j = 1, invN_diag%info%nm
       m = invN_diag%info%ms(j)
       do l = m, lmax
       
          call DRC3JJ(real(l,dp), real(l,dp), 0.d0, 0.d0, l0min, l0max, &
               & threej_symbols_m0, twolmaxp2, ier)
          call DRC3JJ(real(l,dp), real(l,dp), real(-m,dp), real(m,dp), l1min, l1max, &
               & threej_symbols, twolmaxp2, ier)
             
          val = dcmplx(0.d0,0.d0)
          do l2 = 1, nint(l1max-l1min)+1
             lp = nint(l1min) + l2 - 1
             if (lp > lmax) exit
             ! Only m3 = 0 contributes, because m2 = -m1 => m3 = 0
             val = val + a_l0(lp,:) * sqrt(2.d0*lp+1.d0) * &
                  & threej_symbols(l2) * threej_symbols_m0(lp-nint(l0min)+1)
          end do
          val = val * (2*l+1) / sqrt(4.d0*pi) * npix / (4.d0*pi)
          if (mod(m,2)/=0) val = -val

          call invN_diag%info%lm2i(l,m,pos)
          if (m == 0) then
             N_lm(pos,:)   = real(val,dp)
          else
             N_lm(pos,:)   = real(val,dp)
             N_lm(pos+1,:) = real(val,dp)
          end if
       end do
    end do
    !!$OMP END DO
    deallocate(threej_symbols, threej_symbols_m0)
    !!$OMP END PARALLEL
    call wall_time(t2)

    call update_status(status, "compute_invN_lm done")

    invN_diag%alm = N_lm
    !write(*,*) sum(abs(invN_diag%alm))

    deallocate(N_lm, a_l0)
    
  end subroutine compute_invN_lm

  subroutine uniformize_rms(handle, rms, fsky, mask, regnoise)
    implicit none
    type(planck_rng),                   intent(inout) :: handle
    class(comm_map),                    intent(inout) :: rms
    real(dp),                           intent(in)    :: fsky
    class(comm_map),                    intent(in)    :: mask
    real(dp),         dimension(0:,1:), intent(out), optional   :: regnoise

    integer(i4b) :: i, j, nbin=1000, ierr, b
    real(dp)     :: limits(2), dx, threshold, sigma
    real(dp), allocatable, dimension(:) :: F

    if (fsky <= 0.d0) then
       if (present(regnoise)) regnoise = 0.d0
       return
    end if

    allocate(F(nbin))
    do j = 1, rms%info%nmaps
       if (all(mask%map(:,j) < 0.5d0)) cycle
       ! Find pixel histogram across cores
       limits(1) = minval(rms%map(:,j), mask%map(:,j) > 0.5d0)
       limits(2) = maxval(rms%map(:,j), mask%map(:,j) > 0.5d0)
       call mpi_allreduce(MPI_IN_PLACE, limits(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, rms%info%comm, ierr)       
       call mpi_allreduce(MPI_IN_PLACE, limits(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, rms%info%comm, ierr)       
       dx = (limits(2)-limits(1))/nbin
       if (dx == 0.d0) then
          if (present(regnoise)) regnoise(:,j) = 0.d0
          cycle
       end if
       F = 0.d0
       do i = 0, rms%info%np-1
          if (mask%map(i,j) <= 0.5d0) cycle
          b    = max(min(int((rms%map(i,j)-limits(1))/dx),nbin),1)
          F(b) = F(b) + 1.d0
       end do
       call mpi_allreduce(MPI_IN_PLACE, F, nbin, MPI_DOUBLE_PRECISION, MPI_SUM, rms%info%comm, ierr)

       ! Compute cumulative distribution
       do i = 2, nbin
          F(i) = F(i-1) + F(i)
       end do
       F = F / maxval(F)

       ! Find threshold
       i = 1
       do while (F(i) < fsky)
          i = i+1
       end do
       threshold = limits(1) + dx*(i-1)

       ! Update RMS map, and draw corresponding noise realization
       do i = 0, rms%info%np-1
          if (rms%map(i,j) < threshold .and. mask%map(i,j) > 0.5d0) then
             sigma         = sqrt(threshold**2 - rms%map(i,j)**2)
             rms%map(i,j)  = threshold                  ! Update RMS map to requested limit
             if (present(regnoise)) regnoise(i,j) = sigma * rand_gauss(handle) ! Draw corresponding noise realization
          else
             if (present(regnoise)) regnoise(i,j) = 0.d0
          end if
       end do
    end do
    deallocate(F)

  end subroutine uniformize_rms


  subroutine apply_projection(self, map)
    implicit none
    class(comm_N), intent(in)              :: self
    class(comm_map),    intent(inout)           :: map

  end subroutine apply_projection
  
end module comm_N_mod
