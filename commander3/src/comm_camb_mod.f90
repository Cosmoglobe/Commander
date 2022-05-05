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
  use CAMB
  use comm_signal_mod

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
     type(comm_camb_sample), pointer :: sample_curr   ! Current sample 
     type(comm_camb_sample), pointer :: sample_prev   ! Previous sample 

     integer(i4b)     :: lmin               ! Minimum multipole
     integer(i4b)     :: lmax               ! Maximum multipole
     integer(i4b)     :: nalm               ! Number of alms = (l_lmax+1)**2
     integer(i4b)     :: nmaps              ! Number of polarization maps (typically {T,E,B})
     integer(i4b)     :: nspec              ! Number of power spectra (typically six, {TT,TE,TB,EE,EB,BB})
     integer(i4b)     :: nr_of_samples      
     character(len=4) :: dat_length
     real(dp)         :: proposal_multiplier
     real(dp),         dimension(2)   :: noise_l
     real(dp),         dimension(6)   :: correct_cosmo_param ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(6)   :: sigma_cosmo_param   ! [Ombh2, omch2, H0, tau, ln(10^10As), n_s]
     real(dp),         dimension(:, :), allocatable :: L_mat               ! Cholesky factor
     character(len=2), dimension(3)   :: spectra_list
     
   contains
     procedure sample_joint_Cl_theta_sampler
     procedure get_new_sample
     procedure init_covariance_matrix
     procedure get_c_l_from_a_lm
     procedure cosmo_param_proposal
     procedure acceptance
     procedure init_CMB_and_noise
     procedure get_s_lm_f_lm
     procedure get_scaled_f_lm
     procedure get_c_l_from_camb
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
    integer(i4b),              intent(in) :: lmin, lmax
    class(comm_camb), pointer             :: constructor
    type(comm_params)                     :: cpar
    class(comm_camb_sample), pointer      :: old_sample
    type(planck_rng), intent(inout)       :: handle, handle_noise

    integer(i4b) :: samp_group


    allocate(constructor)
    allocate(constructor%L_mat(6, 6))

    samp_group = 1
    constructor%lmin = lmin
    constructor%lmax = lmax
    constructor%nmaps = 3
    constructor%nspec = constructor%nmaps*(constructor%nmaps+1)/2
    constructor%nalm  = (lmax+1)**2
    !constructor%sample_curr => comm_camb_sample(lmin, lmax, constructor%nmaps)
    constructor%sample_curr => constructor_camb_sample(cpar, lmin, lmax, constructor%nmaps)
    constructor%sample_prev => constructor_camb_sample(cpar, lmin, lmax, constructor%nmaps)
    ! Static variables
    constructor%nr_of_samples       = 10
    constructor%dat_length          = '1501'
    constructor%proposal_multiplier = 0.3d0
    constructor%noise_l             = [3.04d-16, 6.08d-16] ! Noise for TT and EE power spectra
    constructor%correct_cosmo_param = [0.02202d0, 0.1153d0, 68.2d0, 0.066d0, 3.035d0, 0.960d0] ! Cosmo Params used to simulate CMB power spectra
    constructor%sigma_cosmo_param   = [0.001d0,   0.005d0,  1.d0,   0.005d0, 0.005d0, 0.005d0] ! Hard coded uncertainty in cosmo param proposal
    constructor%spectra_list        = ['TT', 'EE', 'TE']

    ! Initialize
    constructor%sample_prev%theta = constructor%correct_cosmo_param
    call constructor%get_c_l_from_camb(constructor%sample_prev)
    call constructor%init_covariance_matrix(constructor%L_mat)

    ! First sample, initial guess are the correct cosmological parameters
    constructor%sample_curr%theta = constructor%correct_cosmo_param
    call constructor%get_c_l_from_camb(constructor%sample_curr)

    call constructor%get_s_lm_f_lm(cpar, samp_group, handle, handle_noise, constructor%sample_curr, constructor%sample_curr, .false.)

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
    ! Routine for drawing new samples with Metropolis MCMC and a joint (C_l, s) accept-reject ratio;
    ! see Racine et al. (2016) for details.
    ! 
    ! Arguments
    ! ---------
    ! self:    derived type (comm_camb)
    !          CAMB object
    ! handle: derived type (planck_rng)
    !          Random number handle
    ! 

    implicit none
    class(comm_camb),   intent(inout)       :: self
    type(comm_params), intent(in)                       :: cpar
    type(planck_rng), intent(inout)                        :: handle, handle_noise
    integer(i4b), intent(in)                            :: samp_group
    class(comm_camb_sample), pointer        :: new_sample
    real(dp), dimension(:, :), allocatable  :: cur_sigma_l, hat_c_l, correct_c_l

    integer(i4b) :: nalm, nspec, nmaps, lmax
    logical(lgt) :: accept
    
    lmax  = self%lmax
    nmaps = self%nmaps
    nspec = nmaps*(nmaps+1)/2
    nalm  = (lmax+1)**2

    allocate(cur_sigma_l(nmaps, 0 : lmax))
    allocate(hat_c_l(nmaps, 0 : lmax))
    allocate(new_sample)
    allocate(new_sample%c_l(nmaps, 0 : lmax))

    ! Initialize
    !init_sample%theta = self%correct_cosmo_param
    !write(*, *) self%correct_cosmo_param
    !call self%get_c_l_from_camb(init_sample)
    !call self%init_covariance_matrix(self%L_mat)

    ! First sample, initial guess are the correct cosmological parameters
    !old_sample%theta = self%correct_cosmo_param
    !call self%get_c_l_from_camb(old_sample)
    !call self%get_s_lm_f_lm(cpar, samp_group, handle, handle_noise, d_lm, new_sample)

    call self%get_new_sample(cpar, samp_group, handle, handle_noise, self%sample_curr, self%L_mat, accept, new_sample)

  end subroutine sample_joint_Cl_theta_sampler


   subroutine get_new_sample(self, cpar, samp_group, handle, handle_noise, old_sample, L_mat, accept, new_sample)
   ! 
   ! Gets new cosmological parameter theta sample. Finds c_l from theta, and then s_lm
    ! and scaled f_lm
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
    class(comm_camb),                                 intent(inout) :: self
    logical(lgt),                                     intent(out)   :: accept
    type(comm_camb_sample),                           intent(out)   :: new_sample
    type(comm_camb_sample),                           intent(in)    :: old_sample
    real(dp),         dimension(6, 6),                intent(in)    :: L_mat
    type(comm_params), intent(in)   :: cpar
    type(planck_rng), intent(inout)    :: handle, handle_noise
    integer(i4b), intent(in)        :: samp_group
    integer(i4b) :: ierr
    real(dp), dimension(:, :), allocatable :: scaled_f_lm

    logical(lgt)   :: do_scale_f_lm
    do_scale_f_lm = .true.

    call self%cosmo_param_proposal(old_sample, L_mat, handle, new_sample)
    call self%get_c_l_from_camb(new_sample)
    call self%get_s_lm_f_lm(cpar, samp_group, handle, handle_noise, new_sample, old_sample, do_scale_f_lm)
    call mpi_finalize(ierr)
    ! Here we have to gather all alms. Since we need to scale the fluctuation terms for different ell
    write(*,*) '2 -'
    
    !call self%get_scaled_f_lm(new_sample, old_sample, scaled_f_lm) 
    accept = acceptance(self, new_sample, old_sample, handle)
  end subroutine get_new_sample

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
    
    INQUIRE(FILE="../cosmo_param_out.dat", EXIST=previous_sample)

    if (previous_sample) then
       ! Caclulate covariance matrix
       print *, 'Found previous sample chain. Calculating covariance matrix and use Cholesky decomposition in proposals for this run'
       open(1, file="../cosmo_param_out.dat", status='old', action='read')
       
       nlines = 0
       DO
          READ (1, *, END=10)
          nlines = nlines + 1 
       END DO
10     close(1)
       
       allocate(old_samples(nlines, 6))
       
       open(1, file="../cosmo_param_out.dat", status='old', action='read')
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

  subroutine get_c_l_from_a_lm(self, a_lm, c_l)  
    ! 
    ! Calculates power spectra from a_lm. TT, EE, and TE
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! a_lm: array
    !    alms
    ! Returns
    ! -------
    ! c_l: array
    !    Power spectra caculated from a_lm
    implicit none

    class(comm_camb),                          intent(inout) :: self
    real(dp), dimension(2, (self%lmax+1)**2),  intent(in) :: a_lm
    real(dp), dimension(4, 0:self%lmax),       intent(out) :: c_l

    integer(i4b) :: k, l, m, index
    real(dp)     :: cur_c_l

    c_l = 0.d0
    ! Do TT and EE
    DO k = 1, 2
       DO l = self%lmin, self%lmax
          cur_c_l = a_lm(k, l**2 + l + 1)**2
          DO m = 1, l
             index = l**2 + l + m + 1
             cur_c_l = cur_c_l + 2*a_lm(k, index)**2
          END DO
          c_l(k, l) = cur_c_l / (2*l+1)
       END DO
    END DO
    
    ! Do TE
    DO l = self%lmin, self%lmax
       cur_c_l = a_lm(1, l**2 + l + 1)*a_lm(2, l**2 + l + 1)
       DO m = 1, l
          index = l**2 + l + m + 1
          cur_c_l = cur_c_l + 2*a_lm(1, index)*a_lm(2, index)
       END DO
       c_l(3, l) = cur_c_l / (2*l+1)
    END DO

  end subroutine get_c_l_from_a_lm


  subroutine cosmo_param_proposal(self, old_sample, L, handle, new_sample)
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
    real(dp),               dimension(6, 6), intent(in)    :: L
    type(planck_rng),                        intent(inout) :: handle
    type(comm_camb_sample),                  intent(out)   :: new_sample
    
    real(dp), dimension(6) :: z
    integer(i4b) :: i
    
    do i = 1, 6
       z(i) = self%proposal_multiplier * rand_gauss(handle)
    end do
    new_sample%theta = old_sample%theta + matmul(L, z)
    
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
    !    Scaled fluctutaion term f_lm. See Racine et al. (2016) for details.
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
    class(comm_camb),                                   intent(inout) :: self
    type(comm_camb_sample),                             intent(in)    :: old_sample, new_sample
    type(planck_rng),                                   intent(inout) :: handle

    class(comm_comp), pointer                           :: c => null()
    real(dp), dimension(:, :), allocatable :: old_f_lm, scaled_f_lm
    real(dp), dimension(4, 0: self%lmax) :: old_c_l, new_c_l
    real(dp) :: ln_pi_ip1, ln_pi_i, probability, uni
    real(dp), dimension(2, 2) :: new_S, old_S
    integer(i4b) :: k, i, l, m, comm, ierr
    logical(i4b) :: acceptance

    allocate(old_f_lm(size(old_sample%f_lm,1), size(old_sample%f_lm,2)))
    allocate(scaled_f_lm(size(new_sample%scaled_f_lm,1), size(new_sample%scaled_f_lm,2)))

    old_f_lm = old_sample%f_lm
    old_c_l  = old_sample%c_l
    new_c_l  = new_sample%c_l
    scaled_f_lm = new_sample%scaled_f_lm

    ln_pi_ip1 = 0.d0
    ln_pi_i   = 0.d0

    c => compList
    do while (associated(c))
      select type (c)
      class is (comm_cmb_comp)
        comm = c%x%info%comm
        do i = 0, c%x%info%nalm-1
          l = c%x%info%lm(1, i)
          m = c%x%info%lm(2, i)

          new_S = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
          old_S = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
          call invert_matrix(new_S)
          call invert_matrix(old_S)

          ln_pi_ip1 = ln_pi_ip1 + dot_product(new_sample%s_lm(i, :), matmul(new_S, new_sample%s_lm(i, :)))
          ln_pi_i   = ln_pi_i + dot_product(old_sample%s_lm(i, :), matmul(old_S, old_sample%s_lm(i, :)))
        end do
      end select
      c => c%next()
    end do
  
    !DO l = self%lmin, self%lmax
    !   new_S = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
    !   old_S = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
    !   call invert_matrix(new_S)
    !   call invert_matrix(old_S)
    !   DO m = 0, l   ! HKE: Shouldn't this sum run from -m to m?
    !      i = l**2 + l + m + 1
          
          ! This part is a bit ugly. Everything is diagonal except c_l (because of
          ! C^TE) and so that is done after the k loop. k=1 is a^T_lm and k=2 is
          ! a^E_lm

          !DO k = 1, 2   
          !   ln_pi_ip1 = ln_pi_ip1 + (d_lm(k, i) - new_s_lm(k, i))**2 / self%noise_l(k) + scaled_f_lm(k, i)**2 / self%noise_l(k)
          !   ln_pi_i   = ln_pi_i   + (d_lm(k, i) - old_s_lm(k, i))**2 / self%noise_l(k) + old_f_lm(k, i)**2 / self%noise_l(k)
          !END DO
    !      ln_pi_ip1 = ln_pi_ip1 + dot_product(new_s_lm(:, i), matmul(new_S, new_s_lm(:, i)))
    !      ln_pi_i   = ln_pi_i + dot_product(old_s_lm(:, i), matmul(old_S, old_s_lm(:, i)))
    !   END DO
    !END DO
  
    probability = exp(-(ln_pi_ip1 - ln_pi_i) / 2.0d0)
    
    call mpi_allreduce(MPI_IN_PLACE, probability, 1, MPI_DOUBLE_PRECISION, MPI_PROD, comm, ierr)
    print *, 'prob:', probability
    
    uni = rand_uni(handle) 
    if (uni < probability) then
       acceptance = .true.
       print *, '---------- ACCEPTED -----------'
    else
       acceptance = .false.
       print *, '-------- NOT ACCEPTED ---------'
    end if
  end function acceptance

  subroutine init_CMB_and_noise(self, cur_sample, handle, d_lm)
    ! 
    ! Initializes simulated d_lm which is a_lm from CMB plus noise.
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! cur_sample: derived type (comm_camb_sample)
    !    Sample from which the power spectras are used to simulate d_lm
    ! handle: derived type (planck_rng)
    !    Random number handle
    !
    ! Returns
    ! -------
    ! d_lm: array
    !    Simulated observed alms which are CMB + noise
    ! 
    implicit none
    class(comm_camb),                           intent(inout) :: self
    real(dp),               dimension(2, (self%lmax+1)**2), intent(out)   :: d_lm
    type(comm_camb_sample),                             intent(in)    :: cur_sample
    type(planck_rng),                                   intent(inout) :: handle  
    
    integer(i4b) :: index, i, l, m, k
    real(dp), dimension(4) :: z
    real(dp), dimension(4, 0: self%lmax) :: c_l  
    
    ! Do a^T_lm and a^E_lm
    c_l = cur_sample%c_l
    d_lm = 0.d0
    do l = self%lmin, self%lmax
       do m = 0, l
          index = l**2 + l + m + 1
          z = (/ rand_gauss(handle), rand_gauss(handle), rand_gauss(handle), rand_gauss(handle) /)  
          
          d_lm(1, index) = sqrt(c_l(1, l)) * z(1) + sqrt(self%noise_l(1)) * z(2)
          d_lm(2, index) = c_l(3, l) / sqrt(c_l(1, l)) * z(1) + sqrt(c_l(2, l) - c_l(3, l)**2 / c_l(1, l)) * z(3) + sqrt(self%noise_l(2)) * z(4)
       end do
    end do
    do i = 1, 3
       open(unit=1, file='lcdm_c_'//self%spectra_list(i)//'_l_out.dat', status='replace', action='write')
       write(1, '( '//self%dat_length//'(2X, ES14.6) )') c_l(i, :)
       close(1)
    end do

  end subroutine init_CMB_and_noise

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
    type(comm_camb_sample),               intent(in)    :: old_sample
    logical(lgt)                                        :: do_scale_f_lm
    type(comm_params)                                   :: cpar
    type(planck_rng)                                    :: handle, handle_noise
    integer(i4b)                                        :: samp_group, i, l
    logical(lgt)                                        :: include_mean, include_fluct
    class(comm_comp), pointer                           :: c => null()
    real(dp), dimension(2, 2)                           :: new_S, old_S
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
          allocate(new_sample%s_lm(0:c%x%info%nalm-1, 2))
          new_sample%s_lm = c%x%alm(:, 1:2)
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
          allocate(new_sample%f_lm(0:c%x%info%nalm-1, 2))
          
          new_sample%f_lm = c%x%alm(:, 1:2)

          if (do_scale_f_lm == .true.) then
            allocate(new_sample%scaled_f_lm(0:c%x%info%nalm-1, 2))
            do i = 0, c%x%info%nalm-1
              l = c%x%info%lm(1, i)

              if (l >= 2) then
                new_S = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
                call compute_hermitian_root(new_S, 0.5d0)
                old_S = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
                call compute_hermitian_root(old_S, -0.5d0)

                new_sample%scaled_f_lm(i, :) = matmul(new_S, matmul(old_S, new_sample%f_lm(i, :)))
              end if
            end do
            !new_sample%scaled_f_lm = c%x%alm
          end if
      end select
      c => c%next()
    end do

    !do l = self%lmin, self%lmax
    !   S_inv      = reshape((/ c_l(1, l), c_l(3, l), c_l(3, l), c_l(2, l) /), shape(S_mat))
    !   S_sqrt_inv = S_inv
    !   call invert_matrix(S_inv)
    !   call compute_hermitian_root(S_sqrt_inv, -0.5d0)
    !   N_inv      = reshape((/ self%noise_l(1), 0.d0, 0.d0, self%noise_l(2) /), shape(N))
    !   N_sqrt_inv = N_inv
    !   call invert_matrix(N_inv)
    !   call compute_hermitian_root(N_sqrt_inv, -0.5d0)
    !   common_matrix = S_inv + N_inv
    !   call invert_matrix(common_matrix)
    !   do m = 0, l
    !      index          = l**2 + l + m + 1
    !      d_vector       = [d_lm(1, index), d_lm(2, index)]
    !      omega_1        = [rand_gauss(handle), rand_gauss(handle)]
    !      omega_2        = [rand_gauss(handle), rand_gauss(handle)]
    !      s_lm(:, index) = matmul(common_matrix, matmul(N_inv, d_vector))
    !      f_lm(:, index) = matmul(common_matrix, matmul(S_sqrt_inv, omega_1)) + matmul(common_matrix, matmul(N_sqrt_inv, omega_2))
    !   end do
    !end do
    
  end subroutine get_s_lm_f_lm
  
  subroutine get_scaled_f_lm(self, new_sample, old_sample, scaled_f_lm)
    ! 
    ! Scaled f_lm from new_sample. f_scaled = sqrt(c_l^{i+1}/c_l^i)f^{i+1}
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    ! new_sample: derived type (comm_camb_sample)
    !    Proposed sample with cosmological parameters, CAMB power spectras, s_lm and f_lm
    ! old_sample: derived type (comm_camb_sample)
    !    Previous sample with cosmological parameters, CAMB power spectras, s_lm and f_lm
    !
    ! Returns
    ! -------
    ! scaled_f_lm: array
    !    Scaled fluctutaion term f_lm. See Racine et al. (2016) for details.
    ! 
    implicit none

    class(comm_camb),       intent(inout) :: self
    type(comm_camb_sample), intent(in)    :: new_sample, old_sample

    integer(i4b) :: index, l, m, k, ierr
    real(dp) :: prefactor
    real(dp), dimension(4, 0: self%lmax)      :: old_c_l, new_c_l
    real(dp), dimension(:, :), allocatable    :: old_f_lm
    real(dp), dimension(:, :), allocatable    :: scaled_f_lm
    real(dp), dimension(2, 2)                 :: new_S, old_S, inv_old_S
    
    allocate(old_f_lm(size(old_sample%f_lm,1), size(old_sample%f_lm,2)))
    allocate(scaled_f_lm(size(old_sample%f_lm,1), size(old_sample%f_lm,2)))
    old_f_lm = old_sample%f_lm
    old_c_l  = old_sample%c_l
    new_c_l  = new_sample%c_l

    !do i = 0, info%nalm-1
    !  l = info%lm(i,1)
    !  m = info%lm(i,2)
    !  flm_new%alm(i,1) = sqrt(Cl_new(l)/Cl_old(l)) * flm_old%alm(i,1)
    !end do


    do l = self%lmin, self%lmax
       !Scaling is non-trivial when C^TE != 0, then we need to do matrix operations
       new_S = reshape((/ new_c_l(1, l), new_c_l(3, l), new_c_l(3, l), new_c_l(2, l) /), shape(new_S))
       call compute_hermitian_root(new_S, 0.5d0)
       old_S = reshape((/ old_c_l(1, l), old_c_l(3, l), old_c_l(3, l), old_c_l(2, l) /), shape(old_S))
       call compute_hermitian_root(old_S, -0.5d0)
       write(*,*) 'herher', shape(new_S), shape(old_S), shape(old_f_lm)
       DO m = 0, l ! HKE: Should this loop run from 0 to l?
          index = l**2 + l + m + 1
          scaled_f_lm(:, index) = matmul(new_S, matmul(old_S, old_f_lm(:, index)))
       END DO
    END DO
  end subroutine get_scaled_f_lm
  
  subroutine get_c_l_from_camb(self, new_sample)
    ! 
    ! Gets TT, EE, and TE power spectra from camb using the cosmological
    ! parameters in theta.
    !
    ! Arguments
    ! ---------
    ! self: derived type (comm_camb)
    !    CAMB object
    !
    ! Returns
    ! -------
    ! new_sample: derived type (comm_camb_sample)
    !    Proposed sample with power spectras from CAMB
    ! 
    implicit none
    class(comm_camb),                           intent(inout) :: self
    type(comm_camb_sample) :: new_sample 
    real(dp), dimension(6) :: cosmo_param
    integer(i4b) :: l, k
    real(dp), dimension(4, 0: self%lmax) :: c_l
    
    type(CAMBparams) P
    type(CAMBdata) camb_data
    real(dp) :: CMB_outputscale
    cosmo_param = new_sample%theta
    call CAMB_SetDefParams(P)
    
    P%ombh2 = cosmo_param(1)
    P%omch2 = cosmo_param(2)
    P%omk = 0.d0
    P%H0= cosmo_param(3)
    !call P%set_H0_for_theta(0.0104d0)!cosmo_param(4)
    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = exp(cosmo_param(5))*1e-10
       InitPower%ns = cosmo_param(6)
       InitPower%r = 0.0
    end select
    
    select type(Reion=>P%Reion)
    class is (TTanhReionization)
       Reion%use_optical_depth = .true.
       Reion%optical_depth = cosmo_param(4)
    end select
    
    P%WantScalars           = .true.
    P%WantTensors           = .true.
    P%Accuracy%AccurateBB   = .true.
    P%WantCls               = .true.
    P%DoLensing             = .false.

    P%Max_l=self%lmax+2
    P%Max_eta_k=6000
    P%Max_l_tensor=self%lmax+2
    P%Max_eta_k_tensor=6000

    ! From the CAMB documentation you need this to get micro K^2.
    ! This is (2.726 K * 10^6)^2
    CMB_outputscale = 7.4311e12
    
    call CAMB_GetResults(camb_data, P)
    
    ! Set TT, EE, and TE
    c_l = 0.d0
    !write(*,*) 'camb', shape(camb_data%CLData%Cl_scalar)
    !write(*,*) 'camb k=1', camb_data%CLData%Cl_scalar(2:4, 1)*CMB_outputscale
    !write(*,*) 'camb k=2', camb_data%CLData%Cl_scalar(2:4, 2)*CMB_outputscale
    !write(*,*) 'camb k=3', camb_data%CLData%Cl_scalar(2:4, 3)*CMB_outputscale
    !write(*, *) 'my cl', shape(c_l)
    DO k = 1, 3
       c_l(k, 0) = 0.d0
       c_l(k, 1) = 0.d0
       DO l = 2, self%lmax
          c_l(k, l) = 2.d0 * pi / (l * (l + 1)) * camb_data%CLData%Cl_scalar(l, k) * CMB_outputscale
       END DO
    END DO
    new_sample%c_l = c_l
    
  end subroutine get_c_l_from_camb

end module comm_camb_mod

