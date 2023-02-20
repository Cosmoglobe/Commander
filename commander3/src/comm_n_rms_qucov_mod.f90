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
module comm_N_rms_QUcov_mod
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  use comm_status_mod
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This module handles the case of a white noise matrix where there is
  ! correlation between the Q and U Stokes parameters, but each pixel is
  ! independent. Explicitly, each matrix is represented mathematically as
  ! N^{-1}_{pp} = [\sum 1/\sigma^2, 0,                             0                               ]
  !               [0,              \sum \cos^2 2\psi/\sigma^2,     \sum \cos2\psi\sin2\psi/\sigma^2]
  !               [0,              \sum\cos2\psi\sin2\psi/\sigma^2 \sum \sin^2 2\psi/\sigma^2      ]
  ! If there were no off-diagonal term, the diagonal terms would be 1/\sigma_{T/Q/U}^2.
  !
  ! The N_map object contains four maps, namely the TT, QQ, UU, and QU elements of the per-pixel covariance
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private
  public comm_N_rms_QUcov, comm_N_rms_QUcov_ptr
  
  type, extends (comm_N) :: comm_N_rms_QUcov
     real(dp),        allocatable, dimension(:,:) :: siN, iN, iN_low
     class(comm_map), pointer                     :: N_map => null()
     class(comm_map), pointer                     :: N_low => null()
   contains
     ! Data procedures
     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map_lowres
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS
     procedure :: rms_pix     => returnRMSpix
     procedure :: update_N    => update_N_rms_QUcov
  end type comm_N_rms_QUcov

  interface comm_N_rms_QUcov
     procedure constructor
  end interface comm_N_rms_QUcov

  type comm_N_rms_QUcov_ptr
     type(comm_N_rms_QUcov), pointer :: p => null()
  end type comm_N_rms_QUcov_ptr

  
contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, id_smooth, mask, handle, regnoise, procmask, map)
    implicit none
    class(comm_N_rms_QUcov),                  pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask
    class(comm_map),                    pointer, intent(in), optional :: map

    integer(i4b)       :: i, ierr, nside_smooth
    real(dp)           :: sum_noise, npix
    type(comm_mapinfo), pointer :: info_smooth => null()
    type(comm_mapinfo), pointer :: info_cg => null()
    type(comm_mapinfo), pointer :: info_lowres => null()

    call update_status(status, "comm_N_rms_QUcov constructor")
    
    ! General parameters
    allocate(constructor)

    ! Component specific parameters
    constructor%type              = cpar%ds_noise_format(id_abs)
    constructor%nmaps             = info%nmaps
    constructor%pol               = info%nmaps == 3
    constructor%uni_fsky          = cpar%ds_noise_uni_fsky(id_abs)
    constructor%set_noise_to_mean = cpar%set_noise_to_mean
    constructor%cg_precond        = cpar%cg_precond
    constructor%info              => info

    if (id_smooth == 0) then
       constructor%nside        = info%nside
       constructor%nside_chisq_lowres = min(info%nside, cpar%almsamp_nside_chisq_lowres) ! Used to be n128
       constructor%np           = info%np
       if (cpar%ds_regnoise(id_abs) /= 'none') then
          constructor%rms_reg => comm_map(constructor%info, trim(cpar%ds_regnoise(id_abs)))
       end if
       if (present(procmask)) then
          call constructor%update_N(info, handle, mask, regnoise, procmask=procmask, &
               & noisefile=trim(cpar%ds_noisefile(id_abs)))
       else
          call constructor%update_N(info, handle, mask, regnoise, &
               & noisefile=trim(cpar%ds_noisefile(id_abs)))
       end if
    else
       if (present(map)) then
          constructor%nside        = info%nside
          constructor%nside_chisq_lowres = min(info%nside, cpar%almsamp_nside_chisq_lowres) ! Used to be n128
          constructor%np           = info%np
          call constructor%update_N(info, handle, mask, regnoise, map=map)
       else
          info_smooth => comm_mapinfo(info%comm, nside_smooth, cpar%lmax_smooth(id_smooth), &
               & constructor%nmaps, constructor%pol)
          constructor%nside   = info_smooth%nside
          constructor%np      = info_smooth%np
          constructor%N_map   => comm_map(info_smooth, trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)))
          info_lowres => comm_mapinfo(info%comm, constructor%nside_chisq_lowres, 0, constructor%nmaps, constructor%pol)
          constructor%N_low => comm_map(info_lowres)
          call constructor%N_map%udgrade(constructor%N_low)
          constructor%N_low%map = constructor%N_low%map / (constructor%nside/constructor%nside_chisq_lowres)**2
          allocate(constructor%iN(4,0:info_smooth%np-1))
          allocate(constructor%iN_low(4,0:info_smooth%np-1))
          allocate(constructor%siN(4,0:info_smooth%np-1))
          call initialize_iN_siN(constructor%N_map, constructor%N_low, constructor%iN, constructor%siN, constructor%iN_low)
       end if
    end if

    constructor%pol_only = all(constructor%siN(1,:) == 0.d0)
    call mpi_allreduce(mpi_in_place, constructor%pol_only, 1, MPI_LOGICAL, MPI_LAND, info%comm, ierr)

    call update_status(status, "cg sample group masks")

    ! Initialize CG sample group masks
    info_cg => comm_mapinfo(info%comm, info%nside, 0, 3, info%pol)
    allocate(constructor%samp_group_mask(cpar%cg_num_user_samp_groups+cpar%cs_ncomp)) 
    do i = 1, cpar%cg_num_user_samp_groups
       if (trim(cpar%cg_samp_group_mask(i)) == 'fullsky') cycle
       constructor%samp_group_mask(i)%p => comm_map(info_cg, trim(cpar%cg_samp_group_mask(i)), udgrade=.true.)
       where (constructor%samp_group_mask(i)%p%map > 0.d0)
          constructor%samp_group_mask(i)%p%map = 1.d0
       elsewhere
          constructor%samp_group_mask(i)%p%map = 0.d0
       end where
    end do

  end function constructor

  subroutine update_N_rms_QUcov(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_rms_QUcov),                   intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, j, ierr
    real(dp)     :: sum_tau, sum_tau2, sum_noise, npix
    class(comm_map),     pointer :: invW_tau => null(), N_tmp => null()
    class(comm_mapinfo), pointer :: info_lowres => null()
    class(comm_mapinfo), pointer :: info_pre => null()

    call update_status(status, "update_N_rms_QUcov")

    if (associated(self%rms_reg) .or. self%uni_fsky > 0) then
       call report_error( "Regularization noise not yet supported for rms_QUcov")
    end if

    ! Initialize N
    if (present(noisefile)) then
       self%N_map     => comm_map(info, noisefile)
    else if (present(map)) then
       self%N_map => comm_map(info)
       if (map%info%nmaps == 3) then
           self%N_map%map(:,1:3) = map%map
           self%N_map%map(:,4)   = 0
       else
           self%N_map%map = map%map
       end if
    else
       call report_error('Error in update_N_rms - no noisefile or map declared')
    end if


    ! Add regularization noise; rms_reg is N for this type, not RMS
!!$    if (associated(self%rms_reg)) self%N_map = self%N_map + self%rms_reg 
!!$    call uniformize_rms(handle, self%siN, self%uni_fsky, mask, regnoise)

    ! Apply mask
    if (present(mask)) then
       if (self%pol) then
           self%N_map%map(:,1:3) = self%N_map%map(:,1:3) * mask%map
           self%N_map%map(:,4)   = self%N_map%map(:,4)   * mask%map(:,2)*mask%map(:,3)
       else
           self%N_map%map = self%N_map%map * mask%map
       end if
    end if

    ! Set N to its mean; useful for debugging purposes
    if (self%set_noise_to_mean) then
       do i = 1, self%nmaps
          sum_noise = sum(self%N_map%map(:,i))
          npix      = size(self%N_map%map(:,i))
          call mpi_allreduce(MPI_IN_PLACE, sum_noise,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, npix,       1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          self%N_map%map(:,i) = sum_noise/npix
       end do
    end if


    ! Add white noise corresponding to the user-specified regularization noise map
!!$    if (associated(self%rms_reg) .and. present(regnoise)) then
!!$       write(*,*) 'Warning -- QUcov not accounted for in regnoise'
!!$       do j = 1, size(regnoise,2)
!!$          do i = 0, self%rms_reg%info%np-1
!!$             regnoise(i,j) = regnoise(i,j) + self%rms_reg%map(i,j) * rand_gauss(handle) 
!!$          end do
!!$       end do
!!$    end if
    regnoise = 0d0

    ! Boost noise rms by 20 in processing mask; only for T
    if (present(procmask)) then
       where (procmask%map(:,1) < 0.5d0)
          self%N_map%map(:,1) = self%N_map%map(:,1) * 400.d0 
       end where
    end if


    ! Set up lowres map
    call update_status(status, "N_rms_QUcov lowres")
    if (.not.associated(self%N_low)) then
       info_lowres => comm_mapinfo(self%info%comm, self%nside_chisq_lowres, 0, self%nmaps, self%pol)
       self%N_low => comm_map(info_lowres)
    end if
    call self%N_map%udgrade(self%N_low)
    self%N_low%map = self%N_low%map / (self%nside/self%nside_chisq_lowres)**2

    ! Compute invN and sqrt(invN); both are symmetric 
    if (.not. allocated(self%iN))  allocate(self%iN(4,0:info%np-1))
    if (.not. allocated(self%siN)) allocate(self%siN(4,0:info%np-1))
    if (.not. allocated(self%iN_low)) allocate(self%iN_low(4,0:info%np-1))
    call initialize_iN_siN(self%N_map, self%N_low, self%iN, self%siN, self%iN_low)

    ! Initialize preconditioner noise
    call update_status(status, "N_rms_QUcov cg_precond")
    info_pre => comm_mapinfo(info%comm, info%nside, info%lmax, &
         & min(info%nmaps,3), info%pol)
    if (trim(self%cg_precond) == 'diagonal') then
       ! Set up diagonal covariance matrix
       if (.not. associated(self%invN_diag)) self%invN_diag => comm_map(info_pre)
       do i = 1, info_pre%nmaps
          self%invN_diag%map(:,i) = self%iN(i,:)
       end do
       call compute_invN_lm(self%invN_diag)
    else if (trim(self%cg_precond) == 'pseudoinv') then
       ! Compute alpha_nu for pseudo-inverse preconditioner
       if (.not. allocated(self%alpha_nu)) allocate(self%alpha_nu(self%nmaps))
       invW_tau     => comm_map(info_pre)
       do i = 1, info_pre%nmaps
          invW_tau%map(:,i) = self%iN(i,:)
       end do
       call invW_tau%Yt()
       call invW_tau%Y()
       ! Temperature
       sum_tau  = sum(invW_tau%map(:,1))
       sum_tau2 = sum(invW_tau%map(:,1)**2)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       if (sum_tau > 0.d0) then
          self%alpha_nu(1) = sqrt(sum_tau2/sum_tau)
       else
          self%alpha_nu(1) = 0.d0
       end if

       if (info_pre%nmaps == 3) then
          sum_tau  = sum(invW_tau%map(:,2:3))
          sum_tau2 = sum(invW_tau%map(:,2:3)**2)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau,  1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, sum_tau2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          if (sum_tau > 0.d0) then
             self%alpha_nu(2:3) = sqrt(sum_tau2/sum_tau)
          else
             self%alpha_nu(2:3) = 0.d0
          end if
          call invW_tau%dealloc(); deallocate(invW_tau)
       end if
    end if


  end subroutine update_N_rms_QUcov

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i

    if (self%pol) then
       do i = 0, self%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%iN(1,i) * map%map(i,1)
          map%map(i,2) = self%iN(2,i) * buff_Q + self%iN(4,i) * buff_U
          map%map(i,3) = self%iN(4,i) * buff_Q + self%iN(3,i) * buff_U
       end do
    else
       map%map = map%map * self%iN
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulInvN_1map

  ! Return map_out = invN * map
  subroutine matmulInvN_1map_lowres(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i

    if (self%pol) then
       do i = 0, self%N_low%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%iN_low(1,i) * map%map(i,1)
          map%map(i,2) = self%iN_low(2,i) * buff_Q + self%iN_low(4,i) * buff_U
          map%map(i,3) = self%iN_low(4,i) * buff_Q + self%iN_low(3,i) * buff_U
       end do
   else
       map%map = self%iN_low
   end if
!!$    if (present(samp_group)) then
!!$       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
!!$    end if
  end subroutine matmulInvN_1map_lowres

  ! Return map_out = N * map
  subroutine matmulN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i
    write(*,*) 'N map'
    do i = 0, self%info%np-1
       buff_Q = map%map(i,2)
       buff_U = map%map(i,3)       
       map%map(i,1) = self%N_map%map(i,1) * map%map(i,1)
       map%map(i,2) = self%N_map%map(i,2) * buff_Q + self%N_map%map(i,4) * buff_U
       map%map(i,3) = self%N_map%map(i,4) * buff_Q + self%N_map%map(i,3) * buff_U
    end do
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulN_1map
  
  ! Return map_out = sqrtInvN * map
  subroutine matmulSqrtInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp)     :: buff_Q, buff_U
    integer(i4b) :: i

    !write(*,*) 'matmul sqrtInvN * map'

    if (self%pol) then
       do i = 0, self%info%np-1
          buff_Q = map%map(i,2)
          buff_U = map%map(i,3)       
          map%map(i,1) = self%siN(1,i) * map%map(i,1)
          map%map(i,2) = self%siN(2,i) * buff_Q + self%siN(4,i) * buff_U
          map%map(i,3) = self%siN(4,i) * buff_Q + self%siN(3,i) * buff_U
       end do
    else
       map%map = self%siN * map%map
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulSqrtInvN_1map

!!$  ! Return map_out = invN * map
!!$  subroutine matmulInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_rms_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$    res%map = (self%siN%map)**2 * map%map
!!$  end subroutine matmulInvN_2map
!!$  
!!$  ! Return map_out = sqrtInvN * map
!!$  subroutine matmulSqrtInvN_2map(self, map, res)
!!$    implicit none
!!$    class(comm_N_rms_QUcov), intent(in)              :: self
!!$    class(comm_map),   intent(in)              :: map
!!$    class(comm_map),   intent(inout)           :: res
!!$    res%map = self%siN%map * map%map
!!$  end subroutine matmulSqrtInvN_2map

  ! Return RMS map
  subroutine returnRMS(self, res, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: res
    integer(i4b),      intent(in),   optional  :: samp_group
    where (self%N_map%map(:,1:3) > 0.d0)
       res%map = sqrt(self%N_map%map(:,1:3))
    elsewhere
       res%map = infinity
    end where
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) then
          where (self%samp_group_mask(samp_group)%p%map == 0.d0)
             res%map = infinity
          end where
       end if
    end if
  end subroutine returnRMS
  
  ! Return rms for single pixel
  function returnRMSpix(self, pix, pol, samp_group, ret_invN)
    implicit none
    class(comm_N_rms_QUcov),   intent(in)              :: self
    integer(i4b),        intent(in)              :: pix, pol
    real(dp)                                     :: returnRMSpix
    integer(i4b),        intent(in),   optional  :: samp_group
    logical(lgt),        intent(in),   optional  :: ret_invN
    if (self%N_map%map(pix,pol) > 0.d0) then
       returnRMSpix = sqrt(self%N_map%map(pix,pol))
    else
       returnRMSpix = infinity
    end if
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) then
          if (self%samp_group_mask(samp_group)%p%map(pix,pol) == 0.d0) then
             returnRMSpix = infinity
          end if
       end if
    end if
    if (present(ret_invN)) then
       if (ret_invN) returnRMSpix = 1d0/returnRMSpix**2
    end if
  end function returnRMSpix

  subroutine initialize_iN_siN(N, N_low, iN, siN, iN_low)
    implicit none
    class(comm_map),                   intent(in)    :: N
    class(comm_map),                   intent(in)    :: N_low
    real(dp),        dimension(1:,0:), intent(inout) :: iN
    real(dp),        dimension(1:,0:), intent(inout) :: siN
    real(dp),        dimension(1:,0:), intent(inout) :: iN_low

    integer(i4b) :: i
    real(dp) :: A(2,2)
    logical(lgt)        :: pol

    if (size(N%map, dim=2) .eq. 3) then
      pol = .true.
    else
      pol = .false.
    end if

    do i = 0, N%info%np-1
       ! T component
       if (N%map(i,1) > 0.) then
          iN(1,i)  = 1.d0 / N%map(i,1)
          siN(1,i) = sqrt(iN(1,i))
       else
          iN(1,i)  = 0.d0
          siN(1,i) = 0.d0
       end if

       if (pol) then
          ! QU block; check for positive definite matrix
          if (N%map(i,2)*N%map(i,3)-N%map(i,4)**2 > 0.) then
             A(1,1)  = N%map(i,2) ! QQ
             A(1,2)  = N%map(i,4) ! QU
             A(2,1)  = N%map(i,4) ! UQ
             A(2,2)  = N%map(i,3) ! UU

             call compute_hermitian_root(A, -1.d0)
             iN(2,i)  = A(1,1)    ! QQ
             iN(3,i)  = A(2,2)    ! UU
             iN(4,i)  = A(1,2)    ! QU = UQ

             call compute_hermitian_root(A, 0.5d0)
             siN(2,i) = A(1,1)    ! QQ
             siN(3,i) = A(2,2)    ! UU
             siN(4,i) = A(1,2)    ! QU = UQ
          else
             iN(2:4,i)  = 0.d0
             siN(2:4,i) = 0.d0
          end if
       end if
    end do

    do i = 0, N_low%info%np-1
       ! T component
       if (N_low%map(i,1) > 0.) then
          iN_low(1,i)  = 1.d0 / N_low%map(i,1)
       else
          iN_low(1,i)  = 0.d0
       end if

       if (pol) then
          ! QU block; check for positive definite matrix
          if (N_low%map(i,2)*N_low%map(i,3)-N_low%map(i,4)**2 > 0.) then
             A(1,1)  = N_low%map(i,2) ! QQ
             A(1,2)  = N_low%map(i,4) ! QU
             A(2,1)  = N_low%map(i,4) ! UQ
             A(2,2)  = N_low%map(i,3) ! UU

             call compute_hermitian_root(A, -1.d0)
             iN_low(2,i)  = A(1,1)    ! QQ
             iN_low(3,i)  = A(2,2)    ! UU
             iN_low(4,i)  = A(1,2)    ! QU = UQ
          else
             iN_low(2:4,i)  = 0.d0
          end if
       end if
    end do

  end subroutine initialize_iN_siN

end module comm_N_rms_QUcov_mod
