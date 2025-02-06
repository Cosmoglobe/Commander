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
  use comm_map_mod
  implicit none

  type :: comm_N
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
     procedure :: invN            => matmulInvN
     procedure :: invN_lowres     => matmulInvNlowres
     procedure :: N               => matmulN
     procedure :: sqrtInvN        => matmulSqrtInvN
     procedure :: rms             => returnRMS
     procedure :: rms_pix         => returnRMSpix
     procedure :: update_N        => update_N
     procedure :: P               => apply_projection
  end type comm_N

  type comm_N_ptr
     class(comm_N), pointer :: p => null()
  end type comm_N_ptr
  


contains
     ! Return map_out = invN * map
     subroutine matmulInvN(self, map, samp_group)
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulInvN

     ! Return map_out = invN * map
     subroutine matmulInvNlowres(self, map, samp_group)
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulInvNlowres

     ! Return map_out = N * map
     subroutine matmulN(self, map, samp_group)
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulN

     ! Return map_out = sqrtInvN * map
     subroutine matmulSqrtInvN(self, map, samp_group)
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: map
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine matmulSqrtInvN

     ! Return rms map
     subroutine returnRMS(self, res, samp_group)
       implicit none
       class(comm_N),   intent(in)             :: self
       class(comm_map), intent(inout)          :: res
       integer(i4b),    intent(in),   optional :: samp_group
     end subroutine returnRMS

     ! Return rms map
     function returnRMSpix(self, pix, pol, samp_group, ret_invN)
       implicit none
       class(comm_N),   intent(in)             :: self
       integer(i4b),    intent(in)             :: pix, pol
       real(dp)                                :: returnRMSpix
       integer(i4b),    intent(in),   optional :: samp_group
       logical(lgt),    intent(in),   optional :: ret_invN
       returnRMSpix = infinity
     end function returnRMSpix

     ! Update noise model
     subroutine update_N(self, info, handle, mask, regnoise, procmask, noisefile, map)
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

  subroutine compute_invN_lm(invN_diag)
    implicit none

    class(comm_map),  intent(inout) :: invN_diag

    real(dp)     :: l0min, l0max, l1min, l1max, npix, t1, t2
    integer(i4b) :: j, l, m, lp, l2, ier, twolmaxp2, pos, lmax
    complex(dpc) :: val(invN_diag%info%nmaps)
    real(dp), allocatable, dimension(:,:) :: N_lm, a_l0

    lmax      = invN_diag%info%lmax
    twolmaxp2 = 2*lmax+2
    npix      = real(invN_diag%info%npix,dp)

    allocate(N_lm(0:invN_diag%info%nalm-1,invN_diag%info%nmaps))
    allocate(a_l0(0:lmax,invN_diag%info%nmaps))
    call invN_diag%YtW_scalar


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


  end subroutine uniformize_rms


  subroutine apply_projection(self, map)
    implicit none
    class(comm_N), intent(in)              :: self
    class(comm_map),    intent(inout)           :: map

  end subroutine apply_projection
  
end module comm_N_mod
