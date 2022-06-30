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
module comm_camb_mod
  use comm_utils
  use comm_comp_mod
  use comm_param_mod
  use comm_cmb_comp_mod
  use comm_signal_mod
  use math_tools
  use comm_map_mod
  use comm_cl_mod
  use comm_data_mod

  use comm_camb_eval_mod
  implicit none

  private
  public comm_camb, comm_camb_sample, initialize_camb_mod

  type comm_camb_sample
     real(dp), dimension(:),   allocatable :: theta    ! Cosmological parameters
     real(dp), dimension(:,:), allocatable :: c_l      ! Angular power spectrum
     real(dp), dimension(:,:), allocatable :: s_lm     ! Mean-field alms
     real(dp), dimension(:,:), allocatable :: f_lm     ! Fluctuation alms
     real(dp), dimension(:,:), allocatable :: scaled_f_lm     ! Scaled fluctuation alms
   contains
     procedure :: dealloc => dealloc_camb_sample
     procedure :: equal   => camb_sample_set_equal
  end type comm_camb_sample

  type comm_camb
     type(comm_camb_sample), pointer :: sample_new   ! Current sample 
     type(comm_camb_sample), pointer :: sample_old   ! Previous sample 

     integer(i4b)     :: lmin               ! Minimum multipole
     integer(i4b)     :: lmax               ! Maximum multipole
     integer(i4b)     :: nalm               ! Number of alms = (l_lmax+1)**2
     integer(i4b)     :: nmaps              ! Number of polarization maps (typically {T,E,B})
     integer(i4b)     :: nspec              ! Number of power spectra (typically six, {TT,TE,TB,EE,EB,BB})
     integer(i4b)     :: total_samples      
     integer(i4b)     :: accepted_samples
     real(dp)         :: proposal_multiplier
     real(dp),         dimension(6)   :: correct_cosmo_param ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(6)   :: sigma_cosmo_param   ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(:, :), allocatable :: L_mat               ! Cholesky factor
     
   contains
     procedure sample_joint_Cl_theta_sampler
     procedure init_covariance_matrix
     procedure cosmo_param_proposal
     procedure acceptance
     procedure get_s_lm_f_lm
     procedure compute_fluctuation_acceptance
  end type comm_camb

  interface comm_camb
     procedure constructor
  end interface comm_camb

  interface comm_camb_sample
     procedure constructor_camb_sample
  end interface comm_camb_sample

contains

  function initialize_camb_mod(cpar, handle, handle_noise)
    implicit none
    type(comm_params),         intent(in) :: cpar
    class(comm_camb), pointer             :: initialize_camb_mod
    type(planck_rng), intent(inout)       :: handle, handle_noise

    class(comm_comp), pointer             :: c => null()
    integer(i4b)                          :: lmax

    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
          lmax = c%x%info%lmax
      end select
      c => c%next()
    end do

    initialize_camb_mod => comm_camb(cpar, handle, handle_noise, 2, lmax)
  end function initialize_camb_mod

  function constructor(cpar, handle, handle_noise, lmin, lmax)
    ! 
    ! Constructor for CAMB object
    ! 
    ! Arguments
    ! --------- 
    ! lmin:    int (scalar)
    !          Minimum multipole to be used for parameter fitting
    ! lmax:    int (scalar)
    !          Maximum multipole to be used for parameter fitting
    !
    ! Returns
    ! -------
    ! constructor: 
    !          pointer to CAMB object
    ! 
    implicit none
    class(comm_comp), pointer             :: c => null()
    integer(i4b),              intent(in) :: lmin, lmax
    class(comm_camb), pointer             :: constructor
    type(comm_params)                     :: cpar
    class(comm_camb_sample), pointer      :: old_sample
    type(planck_rng), intent(inout)       :: handle, handle_noise
    real(dp), dimension(4, 0: lmax)       :: init_sample_c_l, first_sample_c_l
    logical(i4b)                          :: finished
    integer(i4b)                          :: ierr, comm, l
    
    real(dp), dimension(6, 6) :: L_mat

    integer(i4b) :: samp_group

    allocate(constructor)
    allocate(constructor%L_mat(6, 6))

    finished = .false.
    samp_group = 1
    constructor%lmin = lmin
    constructor%lmax = lmax
    constructor%nmaps = 3
    constructor%nspec = constructor%nmaps*(constructor%nmaps+1)/2
    constructor%nalm  = (lmax+1)**2
    !constructor%sample_new => comm_camb_sample(lmin, lmax, constructor%nmaps)
    constructor%sample_new => constructor_camb_sample(cpar, lmin, lmax, constructor%nmaps)
    constructor%sample_old => constructor_camb_sample(cpar, lmin, lmax, constructor%nmaps)
    ! Static variables
    constructor%proposal_multiplier = 0.3d0
    constructor%correct_cosmo_param = [0.02237d0, 0.12d0, 67.36d0, 0.0544d0, 3.035d0, 0.9649d0] ! Cosmo Params used to simulate CMB power spectra
    constructor%sigma_cosmo_param   = [0.0005d0,   0.0005d0,  0.3d0,   0.0005d0, 0.01d0, 0.005d0] ! Hard coded uncertainty in cosmo param proposal
    constructor%accepted_samples    = 0
    constructor%total_samples       = 0

    ! Initialize
    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
        comm = c%x%info%comm
        if (c%x%info%myid == 0) then
          ! Root core is running this
          ! First sample
          call get_c_l_from_camb(constructor%correct_cosmo_param, init_sample_c_l)
          
          call constructor%init_covariance_matrix(L_mat)
          
          call mpi_bcast(init_sample_c_l, size(init_sample_c_l),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
          call mpi_bcast(L_mat, size(L_mat),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
          finished = .true.
          call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
        else
          loop: do while (.true.)
            call mpi_bcast(init_sample_c_l, size(init_sample_c_l),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
            call mpi_bcast(L_mat, size(L_mat),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
            call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
            if (finished) exit loop
          end do loop
        end if
      end select
      c => c%next()
    end do

    constructor%L_mat = L_mat
    constructor%sample_old%theta = constructor%correct_cosmo_param
    constructor%sample_old%c_l = init_sample_c_l

    ! First s_lm and f_lm sample
    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_diffuse_comp)
        !{TT,TE,TB,EE,EB,BB}
        ! TT, EE, BB and TE
        
        do l = 2, lmax
          c%Cl%Dl(l, 1) = init_sample_c_l(1, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 2) = init_sample_c_l(4, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 3) = 0d0 !TB=0
          c%Cl%Dl(l, 4) = init_sample_c_l(2, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 5) = 0d0 !EB=0 (bad approximation given that cosmic birefringence is real ;))
          c%Cl%Dl(l, 6) = init_sample_c_l(3, l) * l * (l+1) / (2*pi)
        end do
        call c%Cl%updateS
      end select
      c => c%next()
    end do

    call constructor%get_s_lm_f_lm(cpar, samp_group, handle, handle_noise, constructor%sample_old, constructor%sample_old, .false.)

  end function constructor

  function constructor_camb_sample(cpar, lmin, lmax, nmaps, s_init)
    ! 
    ! Constructor for CAMB sample object
    ! 
    ! Arguments
    ! --------- 
    ! lmin:    int (scalar)
    !          Minimum multipole to be used for parameter fitting
    ! lmax:    int (scalar)
    !          Maximum multipole to be used for parameter fitting
    ! nmaps:   int (scalar)
    !          Number of polarization maps (typically 3; {T,E,B})
    ! s_init:  derived type (comm_camb_sample; optional)
    !          Existing object; if passed, the returned object will be set equal to this
    !
    ! Returns
    ! -------
    ! constructor: 
    !          pointer to CAMB sample object
    ! 
    implicit none
    integer(i4b),                      intent(in)           :: lmin, lmax, nmaps
    class(comm_camb_sample),           intent(in), optional :: s_init
    class(comm_camb_sample),   pointer                      :: constructor_camb_sample
    type(comm_params)                                       :: cpar
    
    integer(i4b) :: nalm, nspec

    allocate (constructor_camb_sample)

    nspec = nmaps*(nmaps+1)/2
    nalm  = (lmax+1)**2

    allocate(constructor_camb_sample%c_l(0:lmax,nspec))
    allocate(constructor_camb_sample%theta(6))

    if (present(s_init)) then
       constructor_camb_sample%c_l   = s_init%c_l
       constructor_camb_sample%theta = cpar%cmb_theta
    else
       constructor_camb_sample%c_l   = 0d0
       constructor_camb_sample%theta = [0.02202d0, 0.1153d0, 68.2d0, 0.066d0, 3.035d0, 0.960d0] 
    end if

  end function constructor_camb_sample

  subroutine dealloc_camb_sample(self)
    ! 
    ! Cleanup routine for comm_camb_sample object
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived type (comm_camb_sample)
    !          Object to be deallocated
    ! 
    implicit none
    class(comm_camb_sample), intent(inout) :: self
  
    if (allocated(self%c_l))  deallocate(self%c_l)
    if (allocated(self%s_lm)) deallocate(self%s_lm)
    if (allocated(self%f_lm)) deallocate(self%f_lm)
    if (allocated(self%scaled_f_lm)) deallocate(self%scaled_f_lm)

  end subroutine dealloc_camb_sample

  subroutine camb_sample_set_equal(self, s_in)
    ! 
    ! Routine for current object equal to s_in
    ! 
    ! Arguments
    ! --------- 
    ! self:    derived type (comm_camb_sample)
    !          Object to overwritten
    ! s_in:    derived type (comm_camb_sample)
    !          Object to copy
    ! 
    implicit none
    class(comm_camb_sample), intent(inout) :: self
    class(comm_camb_sample), intent(in)    :: s_in
  
    self%c_l  = s_in%c_l
    self%s_lm = s_in%s_lm
    self%f_lm = s_in%f_lm
    self%scaled_f_lm = s_in%scaled_f_lm
    
  end subroutine camb_sample_set_equal


   subroutine sample_joint_Cl_theta_sampler(self, cpar, samp_group, handle, handle_noise)
    !
    ! 
    ! Routine for drawing new samples with Metropolis MCMC and a joint (C_l, s) accept-reject ratio;
    ! see Racine et al. (2016) for details.
    ! 
    ! 
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! old_sample: derived type (comm_camb_sample)
    !    Previous sample used to get new cosmological parameters and to
    !    scale fluctuation term f_lm
    ! L_mat: array
    !    Matrix Cholesky decomposition of covariance matrix used for to propose
    !    new cosmological parameters
    ! handle: derived type (planck_rng)
    !    Random number handle
    !
    ! Returns:
    ! --------
    ! accept: bool
    !    True if proposal is accepted
    ! new_sample: derived type (comm_camb_sample)
    !    New sample which includes proposed cosmological parameters
    !    with its corresponding mean field s_lm and (un-scaled) fluctuation f_lm
    ! 
    implicit none
    class(comm_camb),                   intent(inout) :: self
    type(comm_params),                  intent(in)    :: cpar
    type(planck_rng),                   intent(inout) :: handle, handle_noise
    integer(i4b),                       intent(in)    :: samp_group
    logical(lgt)                                      :: accept, finished
    type(comm_camb_sample), pointer                   :: new_sample
    type(comm_camb_sample)                            :: old_sample
    real(dp), dimension(:, :), allocatable            :: scaled_f_lm
    real(dp), dimension(4, 0: self%lmax)              :: new_sample_c_l
    integer(i4b)                                      :: ierr, comm, l
    real(dp), dimension(6)                            :: new_theta_proposal
    class(comm_comp), pointer                         :: c => null()

    allocate(new_sample)
    allocate(new_sample%c_l(self%nmaps, 0 : self%lmax))

    old_sample = self%sample_old
    finished = .false.

    !call self%cosmo_param_proposal(old_sample, self%L_mat, handle, new_sample)
    !call get_c_l_from_camb(new_sample%theta, new_sample_c_l)

    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
        comm = c%x%info%comm
        if (c%x%info%myid == 0) then
          ! Root core is making a proposal and getting CAMB C_\ell
          call self%cosmo_param_proposal(old_sample, handle, new_theta_proposal)
          call get_c_l_from_camb(new_theta_proposal, new_sample_c_l)
          
          call mpi_bcast(new_theta_proposal, size(new_theta_proposal),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
          call mpi_bcast(new_sample_c_l, size(new_sample_c_l),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)

          finished = .true.
          call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
        else
          loop: do while (.true.)
            call mpi_bcast(new_theta_proposal, size(new_theta_proposal),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
            call mpi_bcast(new_sample_c_l, size(new_sample_c_l),  MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)
            call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
            if (finished) exit loop
          end do loop
        end if
      end select
      c => c%next()
    end do

    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_diffuse_comp)
        ! Commander ordering
        !{TT,TE,TB,EE,EB,BB}

        ! CAMB ordering
        ! TT, EE, BB, and TE
        
        do l = 2, self%lmax
          c%Cl%Dl(l, 1) = new_sample_c_l(1, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 2) = new_sample_c_l(4, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 3) = 0
          c%Cl%Dl(l, 4) = new_sample_c_l(2, l) * l * (l+1) / (2*pi)
          c%Cl%Dl(l, 5) = 0
          c%Cl%Dl(l, 6) = new_sample_c_l(3, l) * l * (l+1) / (2*pi)
        end do
        call c%Cl%updateS
      end select
      c => c%next()
    end do
    
    new_sample%theta = new_theta_proposal
    new_sample%c_l = new_sample_c_l

    call self%get_s_lm_f_lm(cpar, samp_group, handle, handle_noise, new_sample, old_sample, do_scale_f_lm=.true.)
    
    accept = acceptance(self, new_sample, old_sample, handle)

    self%total_samples = self%total_samples + 1 
    if (accept) then
      ! Sample was accepted
      self%accepted_samples = self%accepted_samples + 1
      self%sample_old = new_sample
    end if

    ! Save sample
    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
        if (c%x%info%myid == 0) then
          write(*, *) 'accepted samples, total samples, acceptance rate', self%accepted_samples, self%total_samples, self%accepted_samples / self%total_samples 

          open(unit=1, file='cosmo_param_out.dat', position="append", action='write')
            write(1, '( 6(2X, ES14.6) )') self%sample_old%theta
          close(1)
        end if
      end select
      c => c%next()
    end do
    

  end subroutine sample_joint_Cl_theta_sampler

  subroutine compute_fluctuation_acceptance(self, c, comm, res)
    implicit none
    class(comm_camb),                         intent(inout) :: self
    class(comm_cmb_comp),                     intent(in)    :: c
    integer(i4b),                             intent(in)    :: comm
    real(dp),                                 intent(out)   :: res
    class(comm_map), pointer                                :: f_lm_map

    integer(i4b)                                            :: i, ierr, nmaps
    class(comm_mapinfo), pointer                            :: info 

    res = 0d0

    do i = 1, numband
      nmaps =  min(data(i)%info%nmaps, c%nmaps)
      info  => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, data(i)%info%lmax, nmaps, nmaps==3)
      f_lm_map   => comm_map(info)
      call c%x%alm_equal(f_lm_map)

      ! Beam convolution
      call data(i)%B(0)%p%conv(trans=.false., map=f_lm_map)
      call f_lm_map%Y()

      ! Divide by sqrt(N_ell)
      call data(i)%N%sqrtInvN(f_lm_map)

      ! Take the square
      res = res + sum(f_lm_map%map**2)

      call f_lm_map%dealloc(); deallocate(f_lm_map)
    end do

    call mpi_allreduce(MPI_IN_PLACE, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

  end subroutine compute_fluctuation_acceptance


  subroutine cosmo_param_proposal(self, old_sample, handle, new_theta_proposal)
    ! 
    ! Proposal function w. Finds new sample based on covariance
    ! matrix L*L^T. Proposal_multiplier = 0.3 to make sure the proposal theta
    ! becomes too large.
    !
    ! Arguments
    ! ---------
    ! old_sample: derived type (comm_camb_sample)
    !    Previous sample used to get new cosmological parameters and to
    !    scale fluctuation term f_lm
    ! L: array
    !    Matrix Cholesky decomposition of covariance matrix used for to propose
    !    new cosmological parameters
    ! handle: derived type (planck_rng)
    !    Random number handle
    !
    ! Returns:
    ! --------
    ! new_sample: derived type (comm_camb_sample)
    !    New sample which includes proposed cosmological parameters
    ! 
    implicit none
    class(comm_camb),                        intent(inout) :: self
    type(comm_camb_sample),                  intent(in)    :: old_sample
    type(planck_rng),                        intent(inout) :: handle
    real(dp), dimension(6),                  intent(out)   :: new_theta_proposal
    
    real(dp), dimension(6,6) :: L
    real(dp), dimension(6) :: z
    integer(i4b) :: i
    
    
    do i = 1, 6
       z(i) = self%proposal_multiplier * rand_gauss(handle)
    end do
    new_theta_proposal = old_sample%theta + matmul(self%L_mat, z)
    
  end subroutine cosmo_param_proposal


  function acceptance(self, new_sample, old_sample, handle)
    ! 
    ! This function determines if the new sample should be accepted or not.
    ! Assumes no priors and that the proposal is symmetric. Hence 
    ! A = min(1, pi(theta^{i+1})/pi(theta^i))
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! scaled_f_lm: array
    !    Scaled fluctuation term f_lm. See Racine et al. (2016) for details.
    ! new_sample: derived type (comm_camb_sample)
    !    Proposed sample with cosmological parameters, CAMB power spectras, s_lm and f_lm
    ! old_sample: derived type (comm_camb_sample)
    !    Previous sample with cosmological parameters, CAMB power spectras, s_lm and f_lm
    ! handle: derived type (planck_rng)
    !    Random number handle
    !
    ! Returns
    ! -------
    ! acceptance: boolean
    !    Returns true if sample is accepted
    ! 
    implicit none
    class(comm_camb),                    intent(inout) :: self
    type(comm_camb_sample),              intent(inout) :: old_sample, new_sample
    type(planck_rng),                    intent(inout) :: handle

    class(comm_comp), pointer                          :: c => null()
    real(dp), dimension(:, :), allocatable             :: temp, sigma_l, temp_flm
    real(dp), dimension(4, 0: self%lmax)               :: old_c_l, new_c_l
    real(dp)                                           :: ln_pi_ip1, ln_pi_i, cur_i, cur_ip1, chisq_ip1, chisq_i, log_probability, uni, fluc_ip1, fluc_i
    real(dp), dimension(3, 3)                          :: new_S, old_S
    integer(i4b)                                       :: k, i, l, m, comm, ierr
    logical(i4b)                                       :: acceptance, finished
    character(len=512)                                 :: filename

    old_c_l  = old_sample%c_l
    new_c_l  = new_sample%c_l
    finished = .false.

    ln_pi_ip1 = 0.d0
    ln_pi_i   = 0.d0

    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
        comm = c%x%info%comm
        do i = 0, c%x%info%nalm-1
          l = c%x%info%lm(1, i)
          if (l >= 2) then
            m = c%x%info%lm(2, i)
            
            new_S = reshape((/ new_c_l(1, l), new_c_l(4, l), 0d0, new_c_l(4, l), new_c_l(2, l), 0d0, 0d0, 0d0, new_c_l(3, l) /), shape(new_S))
            old_S = reshape((/ old_c_l(1, l), old_c_l(4, l), 0d0, old_c_l(4, l), old_c_l(2, l), 0d0, 0d0, 0d0, new_c_l(3, l) /), shape(old_S))
            
            call invert_matrix(new_S)
            call invert_matrix(old_S)

            cur_ip1 = dot_product(new_sample%s_lm(i, :), matmul(new_S, new_sample%s_lm(i, :)))
            cur_i = dot_product(old_sample%s_lm(i, :), matmul(old_S, old_sample%s_lm(i, :)))

            ln_pi_ip1 = ln_pi_ip1 + cur_ip1
            ln_pi_i   = ln_pi_i + cur_i
          end if
        end do

        
        allocate(temp(size(c%x%alm,1), size(c%x%alm,2)))
        temp = c%x%alm
        c%x%alm = new_sample%s_lm

        allocate(sigma_l(0:c%x%info%lmax,c%x%info%nspec))
        call c%x%getSigmaL(sigma_l)
        if (c%x%info%myid == 0) then
          write(*,*) 'writing to dat', shape(sigma_l), sigma_l(100, 1), sigma_l(100, 1) *100*101/(6.28d0) 
          filename = 'sigma_l_s_hat.dat'
          call write_sigma_l(filename, sigma_l)
        end if
        call compute_chisq(c%x%info%comm, chisq_fullsky=chisq_ip1, evalpol=.true.)

        c%x%alm = old_sample%s_lm
        call compute_chisq(c%x%info%comm, chisq_fullsky=chisq_i, evalpol=.true.)

        ! For plotting purposes
        c%x%alm = new_sample%s_lm+new_sample%f_lm
        call c%x%getSigmaL(sigma_l)
        if (c%x%info%myid == 0) then
          write(*,*) 'writing to dat', shape(sigma_l), sigma_l(100, 1), sigma_l(100, 1) *100*101/(6.28d0) 
          filename = 'sigma_l_s_hat_f_hat.dat'
          call write_sigma_l(filename, sigma_l)
        end if
        c%x%alm = new_sample%f_lm
        call c%x%getSigmaL(sigma_l)
        if (c%x%info%myid == 0) then
          write(*,*) 'writing to dat', shape(sigma_l), sigma_l(100, 1), sigma_l(100, 1) *100*101/(6.28d0) 
          filename = 'sigma_l_f_hat.dat'
          call write_sigma_l(filename, sigma_l)
        end if
      

        c%x%alm = new_sample%scaled_f_lm
        call compute_fluctuation_acceptance(self, c, comm, fluc_ip1)
        c%x%alm = old_sample%f_lm
        call compute_fluctuation_acceptance(self, c, comm, fluc_i)

        c%x%alm = temp

        call mpi_allreduce(MPI_IN_PLACE, ln_pi_ip1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        call mpi_allreduce(MPI_IN_PLACE, ln_pi_i, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        
        log_probability = -(ln_pi_ip1 - ln_pi_i) / 2.0d0
        if (c%x%info%myid == 0) then
          write(*,*) 's_lm cl', ln_pi_ip1, ln_pi_i, log_probability
          write(*,*) 'd-As_lm', chisq_ip1, chisq_i, -(chisq_ip1 - chisq_i)/2.0d0
          write(*,*) 'f_lm', fluc_ip1, fluc_i, -(fluc_ip1 - fluc_i)/2.0d0
          write(*, *) 'new chi^2', ln_pi_ip1 + chisq_ip1 + fluc_ip1
          write(*, *) 'old chi^2', ln_pi_i + chisq_i + fluc_i
          log_probability = log_probability - (chisq_ip1 - chisq_i)/2.0d0 - (fluc_ip1 - fluc_i)/2.0d0
          write(*,*) 'Total log prob:', log_probability
          uni = rand_uni(handle)
          if (log_probability > 0) then
            ! If log_p > 0 then exp(log_p) > 1 and will always be accepted
            ! We do this to avoid overflow
            acceptance = .true.
            write(*,*) '---------- ACCEPTED -----------'
          else if (log_probability < -20) then
            ! If log_p < -20 then exp(log_p) < 2*10^(-9) and we discard the sample
            ! We do this to avoid overflow
            acceptance = .false.
            write(*,*) '-------- NOT ACCEPTED ---------'
          else if (uni < exp(log_probability)) then
            acceptance = .true.
            write(*,*) '---------- ACCEPTED -----------'
          else
            acceptance = .false.
            write(*,*) '-------- NOT ACCEPTED ---------'
          end if
          call mpi_bcast(acceptance, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
          finished = .true.
          call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
        else
          loop: do while (.true.)
            call mpi_bcast(acceptance, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
            call mpi_bcast(finished, 1,  MPI_LOGICAL, 0, c%x%info%comm, ierr)
            if (finished) exit loop
          end do loop
        end if
      end select
      c => c%next()
    end do
    
  end function acceptance

  subroutine get_s_lm_f_lm(self, cpar, samp_group, handle, handle_noise, new_sample, old_sample, do_scale_f_lm)
    ! 
    ! Calculates mean field s_lm and fluctuation f_lm from c_l from new_sample
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! handle: derived type (planck_rng)
    !    Random number handle
    ! 
    ! Returns
    ! -------
    ! new_sample: derived type (comm_camb_sample)
    !    Returns s_lm and f_lm in new_sample
    !

    implicit none
    class(comm_camb),                     intent(inout) :: self
    type(comm_camb_sample),               intent(inout) :: new_sample
    type(comm_camb_sample),               intent(inout)    :: old_sample
    logical(lgt)                                        :: do_scale_f_lm
    type(comm_params)                                   :: cpar
    type(planck_rng)                                    :: handle, handle_noise
    integer(i4b)                                        :: samp_group, i, l, m
    logical(lgt)                                        :: include_mean, include_fluct
    class(comm_comp), pointer                           :: c => null()
    real(dp), dimension(3, 3)                           :: new_S, old_S
    real(dp), dimension(4, 0: self%lmax)                :: old_c_l, new_c_l

    old_c_l  = old_sample%c_l
    new_c_l  = new_sample%c_l

    ! Solve for mean-field map
    include_mean = .true.
    include_fluct = .false.
    call sample_amps_by_CG(cpar, samp_group, handle, handle_noise, include_mean, include_fluct)

    !s_lm = c
    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
          allocate(new_sample%s_lm(0:c%x%info%nalm-1, 3))
          new_sample%s_lm = c%x%alm
      end select
      c => c%next()
    end do

    ! Solve for fluctuation map
    include_mean = .false.
    include_fluct = .true.
    call sample_amps_by_CG(cpar, samp_group, handle, handle_noise, include_mean, include_fluct)
    !s_lm = c
    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
          allocate(new_sample%f_lm(0:c%x%info%nalm-1, 3))
          new_sample%f_lm = c%x%alm

          if (do_scale_f_lm == .true.) then
            allocate(new_sample%scaled_f_lm(0:c%x%info%nalm-1, 3))
            do i = 0, c%x%info%nalm-1
              l = c%x%info%lm(1, i)

              if (l >= 2) then
                new_S = reshape((/ new_c_l(1, l), new_c_l(4, l), 0d0, new_c_l(4, l), new_c_l(2, l), 0d0, 0d0, 0d0, new_c_l(3, l) /), shape(new_S))
                old_S = reshape((/ old_c_l(1, l), old_c_l(4, l), 0d0, old_c_l(4, l), old_c_l(2, l), 0d0, 0d0, 0d0, new_c_l(3, l) /), shape(old_S))
                
                call compute_hermitian_root(new_S, 0.5d0)
                call compute_hermitian_root(old_S, -0.5d0)

                new_sample%scaled_f_lm(i, :) = matmul(new_S, matmul(old_S, old_sample%f_lm(i, :)))
              end if
            end do
          end if
      end select
      c => c%next()
    end do

  end subroutine get_s_lm_f_lm

  subroutine init_covariance_matrix(self, L)
    ! 
    ! Caclulates a multivariate Gaussian for the proposal function w
    ! If there exists a cosmo_param_out.dat file in the parent folder it uses
    ! that to calculate a covariance matrix and then uses Cholesky decomposition
    ! to find the L (L*L^T=Covariance matrix), and so theta_new = theta_old + L*z
    ! where z is a random Gaussian vector
    ! 
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    !
    ! Returns
    ! ---------
    ! L: array
    !    Cholesky decomposition matrix L from covariance matrix (L*L^T)      
    implicit none
    class(comm_camb),                  intent(inout) :: self
    real(dp),         dimension(6, 6), intent(out)   :: L

    real(dp), dimension(6, 6) :: covariance_matrix
    integer(i4b) :: i, j, k, nlines
    real(dp)     :: tot
    logical(lgt) :: previous_sample
    real(dp),    dimension(:, :), allocatable :: old_samples 
    real(dp),    dimension(6)                 :: averages
    
    INQUIRE(FILE="cosmo_param_chain.dat", EXIST=previous_sample)

    if (previous_sample) then
       ! Caclulate covariance matrix
       print *, 'Found previous sample chain. Calculating covariance matrix and use Cholesky decomposition in proposals for this run'
       open(1, file="cosmo_param_chain.dat", status='old', action='read')
       
       nlines = 0
       DO
          READ (1, *, END=10)
          nlines = nlines + 1 
       END DO
10     close(1)
       
       allocate(old_samples(nlines, 6))
       
       open(1, file="cosmo_param_chain.dat", status='old', action='read')
       DO k = 1, nlines
          READ (1, *) old_samples(k, :)
       END DO
       close(1)
       
       averages = sum(old_samples, dim = 1) / nlines
       print *, averages
       DO j = 1,6
          DO k = 1,6
             covariance_matrix(j, k) = 0.d0
             DO i = 1, nlines
                covariance_matrix(j, k) = covariance_matrix(j, k) + (old_samples(i, j) - averages(j))*(old_samples(i, k) - averages(k))
             END DO
          END DO
       END DO
       covariance_matrix = covariance_matrix / (nlines - 1)
    else
       print *, 'Did not find previous sample chain. Using hard-coded diagonal covariance matrix.'
       covariance_matrix = 0.d0
       do i = 1, 6
          covariance_matrix(i, i) = self%sigma_cosmo_param(i)**2
       end do
    end if
    
    L = 0.d0
    DO i = 1, 6
       DO j = 1, i
          tot = 0.d0
          if (i == j) then
             DO k = 1, j
                tot = tot + L(j, k)**2
             END DO
             L(j, j) = sqrt(covariance_matrix(j, j) - tot)
             
          else
             DO k = 1, j
                tot = tot + L(i, k) * L(j, k)
             END DO
             L(i, j) = (covariance_matrix(i, j) - tot) / L(j, j)
          end if
       END DO
    END DO
  end subroutine init_covariance_matrix
end module comm_camb_mod

