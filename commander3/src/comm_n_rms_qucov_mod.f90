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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This module handles the case of a white noise matrix where there is
  ! correlation between the Q and U Stokes parameters, but each pixel is
  ! independent. Explicitly, each matrix is represented mathematically as
  ! N^{-1}_{pp} = [\sum 1/\sigma^2, 0,                             0                               ]
  !               [0,              \sum \cos^2 2\psi/\sigma^2,     \sum \cos2\psi\sin2\psi/\sigma^2]
  !               [0,              \sum\cos2\psi\sin2\psi/\sigma^2 \sum \sin^2 2\psi/\sigma^2      ]
  ! If there were no off-diagonal term, the diagonal terms would be 1/\sigma_{T/Q/U}^2.
  !
  ! The rms0 object contains four maps, the first three being the usual \sigma_{T/Q/U} from 
  ! inverting the diagonal of the weighted hits matrix above, and the fourth one
  ! being the off-diagonal as recorded above. In the event that no off-diagonal
  ! terms are recorded, this should reduce to the usual diagonal case.
  !
  ! The siN matrix is the square root of the inverse of N (N^{-1/2}) and is often used
  ! during the Gaussian sampling routines (see, e.g., Appendix A2 of arXiv:2011.05609).
  ! 
  ! Using the notation of the wikipedia page
  ! https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
  ! where M = ((A,  B), (C, D)), we can write
  !       R = ((A+s, B), (C, D+s))/t,
  ! where s=\sqrt{AD-BC} and t=\sqrt{(A+D)+2s}.
  !
  ! So we need to implement the N^{-1}m and N^{-1/2}m operations, specifically
  ! for QU.
  !  
  ! N^{-1} = ((\sigma_Q^{-2}, \rho), (\rho, \sigma_U^{-2}))
  !
  ! N^{-1}m = (\sigma_Q^{-2} Q + \rho U, \sigma_U^{-2} U + \rho Q)
  !
  ! s = \sqrt{\sigma_Q^{-2}\sigma_U^{-2} - \rho^2}
  ! t = \sqrt{\sigma_Q^{-2} + \sigma_U^{-2} + 2s}
  !
  ! N^{-1/2} = ((\sigma_Q^{-2} + s, \rho), (\rho, \sigma_U^{-2} + s))/t
  !
  ! N^{-1/2}m = ((\sigma_Q^{-2}+s)Q + \rho U, (\sigma_U^{-2}+s)U + \rho Q))/t
  !
  ! siN is the N^{-1/2} term.
  !
  ! A trickier one is computing Nm
  ! N = ((\sigma_U^{-2}, -\rho), (-\rho, \sigma_Q^{-2}))/\sqrt{(\sigma_Q\sigma_U)^{-2}-\rho^2}
  !
  ! Nm = (Q\sigma_U^{-2} - U \rho, U\sigma_U^{-2} - Q\rho)/\sqrt{(\sigma_Q\sigma_U)^{-2}-\rho^2}
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use comm_N_mod
  use comm_param_mod
  use comm_map_mod
  implicit none

  private
  public comm_N_rms_QUcov, comm_N_rms_QUcov_ptr
  
  type, extends (comm_N) :: comm_N_rms_QUcov
     class(comm_map), pointer :: siN        => null()
     class(comm_map), pointer :: siN_lowres => null()
     class(comm_map), pointer :: rms0       => null()
   contains
     ! Data procedures

     procedure :: invN        => matmulInvN_1map
     procedure :: invN_lowres => matmulInvN_1map_lowres
     procedure :: N           => matmulN_1map
     procedure :: sqrtInvN    => matmulSqrtInvN_1map
     procedure :: rms         => returnRMS
     procedure :: rms_pix     => returnRMSpix
     procedure :: update_N    => update_N_rms
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
    class(comm_N_rms_QUcov),               pointer       :: constructor
    type(comm_params),                  intent(in)    :: cpar
    type(comm_mapinfo), target,         intent(in)    :: info
    integer(i4b),                       intent(in)    :: id, id_abs, id_smooth
    class(comm_map),                    intent(in)    :: mask
    type(planck_rng),                   intent(inout) :: handle
    real(dp), dimension(0:,1:),         intent(out),         optional :: regnoise
    class(comm_map),                    pointer, intent(in), optional :: procmask
    class(comm_map),                    pointer, intent(in), optional :: map

    class(comm_N_rms_QUcov),               pointer       :: constructor_temp
    integer(i4b)       :: i, j, ierr, tmp, nside_smooth
    real(dp)           :: sum_noise, npix, buffer
    type(comm_mapinfo), pointer :: info_smooth => null()


    ! General parameters
    allocate(constructor)
    allocate(constructor_temp)

    ! Component specific parameters
    constructor%type              = cpar%ds_noise_format(id_abs)
    constructor%nmaps             = info%nmaps + 1
    constructor%pol               = info%nmaps == 4 ! sigmaT, sigmaQ, sigmaU, Cov(Q,U)
    constructor%uni_fsky          = cpar%ds_noise_uni_fsky(id_abs)
    constructor%set_noise_to_mean = cpar%set_noise_to_mean
    constructor%cg_precond        = cpar%cg_precond
    constructor%info              => info
    constructor%info%nmaps        = constructor%nmaps

    call update_status(status, 'reading in masks and stuff')

    if (id_smooth == 0) then
       constructor%nside        = info%nside
       constructor%nside_chisq_lowres = min(info%nside, cpar%almsamp_nside_chisq_lowres) ! Used to be n128
       constructor%np           = info%np
       if (cpar%ds_regnoise(id_abs) /= 'none') then
          constructor%rms_reg => comm_map(constructor%info, trim(cpar%ds_regnoise(id_abs)))
       end if
       call update_status(status, 'updating N')
       if (present(procmask)) then
          call constructor%update_N(constructor%info, handle, mask, regnoise, procmask=procmask, &
               & noisefile=trim(cpar%ds_noisefile(id_abs)))
       else
          call constructor%update_N(constructor%info, handle, mask, regnoise, &
               & noisefile=trim(cpar%ds_noisefile(id_abs)))
       end if
    else
       if (present(map)) then
          constructor%nside        = info%nside
          constructor%nside_chisq_lowres = min(info%nside, cpar%almsamp_nside_chisq_lowres) ! Used to be n128
          constructor%np           = info%np
          call constructor%update_N(info, handle, mask, regnoise, map=map)
       else
          tmp         =  int(getsize_fits(trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)), nside=nside_smooth), i4b)
          info_smooth => comm_mapinfo(info%comm, nside_smooth, cpar%lmax_smooth(id_smooth), &
               & constructor%nmaps + 1, constructor%pol)
          constructor%nside   = info_smooth%nside
          constructor%np      = info_smooth%np
          constructor%siN     => comm_map(info_smooth, trim(cpar%ds_noise_rms_smooth(id_abs,id_smooth)))
          
       end if
    end if

    call update_status(status, 'checking if pol only')

    constructor%pol_only = all(constructor%siN%map(:,1) == 0.d0)
    call mpi_allreduce(mpi_in_place, constructor%pol_only, 1, MPI_LOGICAL, MPI_LAND, info%comm, ierr)

    call update_status(status, 'checking CG stuff')
    ! Initialize CG sample group masks
    allocate(constructor%samp_group_mask(cpar%cg_num_user_samp_groups+cpar%cs_ncomp)) !had to add number og active components so that the array is long enough for the unique sample groups
    constructor_temp => constructor
    constructor_temp%info%nmaps = constructor_temp%info%nmaps - 1
    do i = 1, cpar%cg_num_user_samp_groups
       if (trim(cpar%cg_samp_group_mask(i)) == 'fullsky') cycle
       constructor%samp_group_mask(i)%p => comm_map(constructor_temp%info, trim(cpar%cg_samp_group_mask(i)), udgrade=.true.)
       where (constructor%samp_group_mask(i)%p%map > 0.d0)
          constructor%samp_group_mask(i)%p%map = 1.d0
       elsewhere
          constructor%samp_group_mask(i)%p%map = 0.d0
       end where
    end do

    call update_status(status, 'should have finished constructing')

  end function constructor

  subroutine update_N_rms(self, info, handle, mask, regnoise, procmask, noisefile, map)
    implicit none
    class(comm_N_rms_QUcov),                intent(inout)          :: self
    class(comm_mapinfo),                 intent(in)             :: info
    type(planck_rng),                    intent(inout)          :: handle
    class(comm_map),                     intent(in),   optional :: mask
    real(dp),          dimension(0:,1:), intent(out),  optional :: regnoise
    class(comm_map),                     intent(in),   optional :: procmask
    character(len=*),                    intent(in),   optional :: noisefile
    class(comm_map),                     intent(in),   optional :: map

    integer(i4b) :: i, j, ierr
    real(dp)     :: sum_tau, sum_tau2, sum_noise, npix, buffer
    class(comm_map),     pointer :: invW_tau => null(), iN => null()
    class(comm_mapinfo), pointer :: info_lowres => null()

    if (present(noisefile)) then
       self%rms0     => comm_map(info, noisefile)
    else if (present(map)) then
       !info%nmaps    = info%nmaps + 1
       !self%rms0     => comm_map(info)
       if (size(map%map, dim=2) == 3) then
           self%rms0%map(:,1:3) = map%map
           self%rms0%map(:,4) = 0
       else
           self%rms0%map = map%map
       end if
    else
       call report_error('Error in update_N_rms - no noisefile or map declared')
    end if
    if (associated(self%siN)) then
       self%siN%map = self%rms0%map
    else
       self%siN     => comm_map(self%rms0)
    end if
    if (associated(self%rms_reg)) then
       self%siN%map = sqrt(self%siN%map**2 + self%rms_reg%map**2) 
    end if
    call uniformize_rms(handle, self%siN, self%uni_fsky, mask, regnoise)
    ! self%siN%map = self%siN%map * mask%map ! Apply mask
    self%siN%map(:,1:3) = self%siN%map(:,1:3) * mask%map ! Apply mask
    self%siN%map(:,4)   = self%siN%map(:,4) * mask%map(:,3) ! Apply mask
    ! Store the diagonal version, and the inverted 2x2 matrix
    if (present(procmask)) then
       where (procmask%map < 0.5d0)
          self%siN%map = self%siN%map * 20.d0 ! Boost noise by 20 in processing mask
       end where
    end if


    ! Add white noise corresponding to the user-specified regularization noise map
    if (associated(self%rms_reg) .and. present(regnoise)) then
       do j = 1, self%rms_reg%info%nmaps
          do i = 0, self%rms_reg%info%np-1
             regnoise(i,j) = regnoise(i,j) + self%rms_reg%map(i,j) * rand_gauss(handle) 
          end do
       end do
    end if




    if (trim(self%cg_precond) == 'diagonal') then
       ! Set up diagonal covariance matrix
       if (.not. associated(self%invN_diag)) self%invN_diag => comm_map(info)
       !write(*,*) info%nmaps, shape(self%invN_diag%map), 'basdfa'
       self%invN_diag%map = self%siN%map(:,1:3)**2
       !call compute_invN_lm(self%invN_diag)
    else if (trim(self%cg_precond) == 'pseudoinv') then
       ! Compute alpha_nu for pseudo-inverse preconditioner
       if (.not. allocated(self%alpha_nu)) allocate(self%alpha_nu(self%nmaps))
       invW_tau     => comm_map(self%siN)
       invW_tau%map =  invW_tau%map**2
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

       if (self%nmaps == 3) then
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


    ! Set up lowres map
    if (.not.associated(self%siN_lowres)) then
       info_lowres => comm_mapinfo(self%info%comm, self%nside_chisq_lowres, 0, self%nmaps, self%pol)
       self%siN_lowres => comm_map(info_lowres)
    end if

    iN => comm_map(self%siN)
    iN%map = iN%map**2
    call iN%udgrade(self%siN_lowres)
    call iN%dealloc(); deallocate(iN)
    self%siN_lowres%map = sqrt(self%siN_lowres%map) * (self%nside/self%nside_chisq_lowres)

  end subroutine update_N_rms

  ! Return map_out = invN * map
  subroutine matmulInvN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp), dimension(:), allocatable :: buff_Q, buff_U
    integer(i4b) :: npix
    npix = size(map%map(:,1))
    allocate(buff_Q(npix), buff_U(npix))
    buff_Q = map%map(:,2)
    buff_U = map%map(:,3)
    map%map(:,1) = (self%rms0%map(:,1))**(-2) * map%map(:,1)
    map%map(:,2) = self%rms0%map(:,2)**(-2) * buff_Q + &
                 & self%rms0%map(:,4)       * buff_U
    map%map(:,3) = self%rms0%map(:,3)**(-2) * buff_U + &
                 & self%rms0%map(:,4)       * buff_Q
    deallocate(buff_Q, buff_U)
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
    ! Implement same thing as above, just downgrade the matrices appropriately.
    map%map = (self%siN_lowres%map)**2 * map%map
    if (present(samp_group)) then
       if (associated(self%samp_group_mask(samp_group)%p)) map%map = map%map * self%samp_group_mask(samp_group)%p%map
    end if
  end subroutine matmulInvN_1map_lowres

  ! Return map_out = N * map
  subroutine matmulN_1map(self, map, samp_group)
    implicit none
    class(comm_N_rms_QUcov), intent(in)              :: self
    class(comm_map),   intent(inout)           :: map
    integer(i4b),      intent(in),   optional  :: samp_group
    real(dp), dimension(:), allocatable :: buff_Q, buff_U
    integer(i4b) :: npix
    npix = size(map%map(:,1))
    allocate(buff_Q(npix), buff_U(npix))
    buff_Q = map%map(:,2)
    buff_U = map%map(:,3)
    ! Nm = (Q\sigma_U^{-2} - U \rho, U\sigma_U^{-2} - Q\rho)/\sqrt{(\sigma_Q\sigma_U)^{-2}-\rho^2}
    map%map(:,1) = map%map(:,1) / self%rms0%map(:,1)**2
    map%map(:,2) = (buff_Q/self%rms0%map(:,3)**2 - buff_U*self%rms0%map(:,4)) / &
                 & (1/(self%rms0%map(:,2)*self%rms0%map(:,3))**2 - self%rms0%map(:,4)**2)**0.5
    map%map(:,3) = (buff_U/self%rms0%map(:,2)**2 - buff_Q*self%rms0%map(:,4)) / &
                 & (1/(self%rms0%map(:,2)*self%rms0%map(:,3))**2 - self%rms0%map(:,4)**2)**0.5

    deallocate(buff_Q, buff_U)


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
    integer(i4b) :: npix, i, j
    real(dp), dimension(:), allocatable :: buff_Q, buff_U, s, t
    npix = size(map%map(:,1))
    allocate(buff_Q(npix), buff_U(npix), s(npix), t(npix))
    buff_Q = map%map(:,2)
    buff_U = map%map(:,3)
    ! s = \sqrt{\sigma_Q^{-2}\sigma_U^{-2} - \rho^2}
    ! t = \sqrt{\sigma_Q^{-2} + \sigma_U^{-2} + 2s}
    !
    ! N^{-1/2} = ((\sigma_Q^{-2} + s, \rho), (\rho, \sigma_U^{-2} + s))/t
    !
    ! N^{-1/2}m = ((\sigma_Q^{-2}+s)Q + \rho U, (\sigma_U^{-2}+s)U + \rho Q))/t

    map%map(:,1) = self%rms0%map(:,1) * map%map(:,1)
    where (self%rms0%map(:,2) > 0)
      s = ((self%rms0%map(:,2)*self%rms0%map(:,3))**(-2) - self%rms0%map(:,4)**2)**0.5
      t = (self%rms0%map(:,2)**(-2) + self%rms0%map(:,3)**(-2) + 2*s)
      map%map(:,2) = ((self%rms0%map(:,2)**(-2) + s)*buff_Q + self%rms0%map(:,4)*buff_U)/t
      map%map(:,3) = ((self%rms0%map(:,3)**(-2) + s)*buff_U + self%rms0%map(:,4)*buff_Q)/t
    end where
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
    where (self%siN%map > 0.d0)
       res%map = 1.d0/self%siN%map
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
  function returnRMSpix(self, pix, pol, samp_group)
    implicit none
    class(comm_N_rms_QUcov),   intent(in)              :: self
    integer(i4b),        intent(in)              :: pix, pol
    real(dp)                                     :: returnRMSpix
    integer(i4b),        intent(in),   optional  :: samp_group
    if (self%siN%map(pix,pol) > 0.d0) then
       returnRMSpix = 1.d0/self%siN%map(pix,pol)
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
  end function returnRMSpix

end module comm_N_rms_QUcov_mod
