module comm_like_model_mod
  use comm_like_utils
  use alm_tools
  use rngmod
  implicit none

  ! General parameters
  integer(i4b)                                       :: lmax, nspec, n_h
  integer(i4b)                                       :: nmaps, numcomp
  integer(i4b)                                       :: verbosity
  real(dp),     allocatable, dimension(:,:), private :: cl_fid
  logical(lgt), allocatable, dimension(:),   private :: inc_spec
  integer(i4b), allocatable, dimension(:)            :: indmap

  ! General model parameters
  character(len=128) :: model
  integer(i4b)       :: npar, lmin_fit, lmax_fit, ind_min, ind_max

  real(dp)                                       :: dir_prop_radius
  character(len=15), allocatable, dimension(:)   :: par_label
  real(dp),          allocatable, dimension(:,:) :: prior_uni, prior_gauss
  real(dp),          allocatable, dimension(:)   :: par_modulo, par_rms
  real(dp),          allocatable, dimension(:,:) :: L_prop_in

  ! Specific model parameters
  real(dp),          private :: l_pivot
  

contains

  subroutine initialize_model_mod(parfile)
    implicit none

    character(len=*), intent(in) :: parfile

    integer(i4b)       :: ind1, ind2, i, l, m, n, unit
    logical(lgt)       :: polarization, exist
    character(len=256) :: clfile, propmat_file, likefile

    unit = comm_getlun()
    call comm_get_parameter(parfile,  'LIKELIHOOD_PARAMETER_FILE', par_string=likefile)
    call comm_get_parameter(parfile,  'LMAX',                      par_int=lmax)
    call comm_get_parameter(likefile, 'LMIN',                      par_int=lmin_fit)
    call comm_get_parameter(likefile, 'LMAX',                      par_int=lmax_fit)
    call comm_get_parameter(parfile,  'MODEL',                     par_string=model)
    call comm_get_parameter(parfile,  'FIDUCIAL_CL_FILE',          par_string=clfile)
    call comm_get_parameter(parfile,  'POLARIZATION',              par_lgt=polarization)
    call comm_get_parameter(parfile,  'INPUT_PROPOSAL_MATRIX',     par_string=propmat_file)
    call comm_get_parameter(parfile, 'VERBOSITY',                  par_int=verbosity)
    ind_min = lmin_fit**2+1
    ind_max = (lmax_fit+1)**2

    if (polarization) then
       nmaps = 3
       nspec = 6
       allocate(inc_spec(nspec))
       call comm_get_parameter(parfile, 'INCLUDE_TT', par_lgt=inc_spec(1))       
       call comm_get_parameter(parfile, 'INCLUDE_TE', par_lgt=inc_spec(2))       
       call comm_get_parameter(parfile, 'INCLUDE_TB', par_lgt=inc_spec(3))       
       call comm_get_parameter(parfile, 'INCLUDE_EE', par_lgt=inc_spec(4))       
       call comm_get_parameter(parfile, 'INCLUDE_EB', par_lgt=inc_spec(5))       
       call comm_get_parameter(parfile, 'INCLUDE_BB', par_lgt=inc_spec(6))       
    else
       nmaps = 1
       nspec = 1
       allocate(inc_spec(nspec))
       inc_spec = .true.
    end if
    numcomp = (lmax+1)**2
    n_h     = numcomp*nmaps

    ! Set up index map, ie., list of active basis vectors
    allocate(indmap(nmaps*(ind_max-ind_min+1)))
    ind1 = 1
    ind2 = 1
    do i = 1, nmaps
       do l = 0, lmax
          do m = -l, l
             if (l >= lmin_fit .and. l <= lmax_fit) then
                indmap(ind2) = ind1
                ind2         = ind2+1
             end if
             ind1 = ind1+1
          end do
       end do
    end do
    
    if (trim(model) == 'qn') then

       npar = 5
       allocate(par_label(npar), prior_uni(npar,2), prior_gauss(npar,2), &
            & par_modulo(npar), par_rms(npar))
       par_label        = ['phi  ',   'theta',   'psi  ',    'q    ', 'n    ']
       prior_uni(:,1)   = [-1000.d0,   0.d0,     -1000.d0,    0.01d0, -10.d0]
       prior_uni(:,2)   = [ 1000.d0,   pi,        1000.d0,   10.d0,    10.d0]
       prior_gauss(:,1) = [  0.d0,     0.d0,      0.d0,       0.d0,     0.d0]
       prior_gauss(:,2) = [  0.d0,     0.d0,      0.d0,       0.d0,     0.d0]
       par_modulo       = [ 2.d0*pi,   pi,        2.d0*pi,   10.d0,    10.d0]
       par_rms          = [    0.d0,   0.d0,         0.d0,    0.1d0,    0.1d0]
       call comm_get_parameter(parfile, 'Q_GAUSS_PRIOR_MEAN', par_dp=prior_gauss(4,1))
       call comm_get_parameter(parfile, 'Q_GAUSS_PRIOR_RMS',  par_dp=prior_gauss(4,2))
       call comm_get_parameter(parfile, 'Q_PROP_RMS',         par_dp=par_rms(4))
       call comm_get_parameter(parfile, 'N_GAUSS_PRIOR_MEAN', par_dp=prior_gauss(5,1))
       call comm_get_parameter(parfile, 'N_GAUSS_PRIOR_RMS',  par_dp=prior_gauss(5,2))
       call comm_get_parameter(parfile, 'N_PROP_RMS',         par_dp=par_rms(5))
       call comm_get_parameter(parfile, 'L_PIVOT',            par_dp=l_pivot)
       dir_prop_radius = 0.d0

    else if (trim(model) == 'power_asymmetry') then

       npar = 6
       allocate(par_label(npar), prior_uni(npar,2), prior_gauss(npar,2), &
            & par_modulo(npar), par_rms(npar))
       par_label        = ['phi  ',   'theta',    'psi  ',    'q    ', 'n    ', 'alpha']
       prior_uni(:,1)   = [-1000.d0,   0.d0,     -1000.d0,    0.01d0, -10.d0,   0.d0]
       prior_uni(:,2)   = [ 1000.d0,   pi,        1000.d0,   10.d0,    10.d0,   1.d0]
       prior_gauss(:,1) = [  0.d0,     0.d0,      0.d0,       0.d0,     0.d0,   0.d0]
       prior_gauss(:,2) = [ -1.d0,    -1.d0,      0.d0,       0.d0,     0.d0,   0.d0]
       par_modulo       = [ 2.d0*pi,   pi,        2.d0*pi,   10.d0,    10.d0,   1.d0]
       par_rms          = [  0.d0,     0.d0,      0.d0,      0.1d0,    0.1d0,  0.01d0]
       call comm_get_parameter(parfile, 'Q_GAUSS_PRIOR_MEAN',     par_dp=prior_gauss(4,1))
       call comm_get_parameter(parfile, 'Q_GAUSS_PRIOR_RMS',      par_dp=prior_gauss(4,2))
       call comm_get_parameter(parfile, 'Q_PROP_RMS',             par_dp=par_rms(4))
       call comm_get_parameter(parfile, 'N_GAUSS_PRIOR_MEAN',     par_dp=prior_gauss(5,1))
       call comm_get_parameter(parfile, 'N_GAUSS_PRIOR_RMS',      par_dp=prior_gauss(5,2))
       call comm_get_parameter(parfile, 'N_PROP_RMS',             par_dp=par_rms(5))
       call comm_get_parameter(parfile, 'ALPHA_GAUSS_PRIOR_MEAN', par_dp=prior_gauss(6,1))
       call comm_get_parameter(parfile, 'ALPHA_GAUSS_PRIOR_RMS',  par_dp=prior_gauss(6,2))
       call comm_get_parameter(parfile, 'ALPHA_PROP_RMS',         par_dp=par_rms(6))
       call comm_get_parameter(parfile, 'DIR_PROP_RADIUS',        par_dp=dir_prop_radius)
       call comm_get_parameter(parfile, 'L_PIVOT',                par_dp=l_pivot)
       dir_prop_radius = dir_prop_radius * pi/180.d0

    else
       write(*,*) 'Unknown model: ', trim(model)
       stop
    end if

    ! Read in fiducial spectrum, to be conditioned upon outside range of interest
    call read_fiducial_spectrum(clfile, cl_fid)

    ! Read initial proposal matrix if it exists
    inquire(file=trim(propmat_file), exist=exist)
    if (exist) then
       open(unit,file=trim(propmat_file),form='unformatted')
       read(unit) n
       if (n == npar) then
          allocate(L_prop_in(npar,npar))
          read(unit) L_prop_in
       end if
       close(unit)
    end if

  end subroutine initialize_model_mod

  ! S_hat wrapper
  subroutine compute_S_hat(p, S_hat)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: p
    real(dp), dimension(1:,1:), intent(out) :: S_hat

    if (trim(model) == 'qn') then
       call compute_S_hat_qn(p, S_hat)
    else if (trim(model) == 'power_asymmetry') then
       call compute_S_hat_power_asymmetry(p, S_hat)
    end if

  end subroutine compute_S_hat

  !==================================================================

  subroutine compute_S_hat_qn(p, S_hat)
    implicit none

    real(dp), dimension(:),   intent(in)  :: p
    real(dp), dimension(:,:), intent(out) :: S_hat

    integer(i4b) :: ind, l, m, i
    real(dp) :: q, n, scale

    q = p(1) ! Power spectrum amplitude at l=l_pivot
    n = p(2) ! Power spectrum tilt

    S_hat = 0.d0
    ind   = lmin_fit**2+1
    do l = lmin_fit, lmax_fit
       if (l <= 0) then
          scale = 0.d0
       else
          scale = q*(l/l_pivot)**n / (l*(l+1.d0)/(2.d0*pi))
       end if
       do m = -l, l
          if (inc_spec(TT)) S_hat(0*numcomp+ind,0*numcomp+ind) = scale * cl_fid(l,TT)
          if (inc_spec(TE)) then
             S_hat(1*numcomp+ind,0*numcomp+ind) = scale * cl_fid(l,TE)
             S_hat(0*numcomp+ind,1*numcomp+ind) = scale * cl_fid(l,TE)
          end if
          if (inc_spec(TB)) then
             S_hat(2*numcomp+ind,0*numcomp+ind) = scale * cl_fid(l,TB)
             S_hat(0*numcomp+ind,2*numcomp+ind) = scale * cl_fid(l,TB)
          end if
          if (inc_spec(EE)) S_hat(1*numcomp+ind,1*numcomp+ind) = scale * cl_fid(l,EE)
          if (inc_spec(EB)) then
             S_hat(2*numcomp+ind,1*numcomp+ind) = scale * cl_fid(l,EB)
             S_hat(1*numcomp+ind,2*numcomp+ind) = scale * cl_fid(l,EB)
          end if
          if (inc_spec(BB)) S_hat(2*numcomp+ind,2*numcomp+ind) = scale * cl_fid(l,BB)
          ind = ind+1
       end do
    end do

  end subroutine compute_S_hat_qn


  subroutine compute_S_hat_power_asymmetry(p, S_hat)
    implicit none

    real(dp), dimension(:),   intent(in)  :: p
    real(dp), dimension(:,:), intent(out) :: S_hat

    integer(i4b) :: i, j, k, l1, l2, m1, m2, ind1, ind2, n, nside, npix
    real(dp)     :: alpha, lstar, theta, phi
    real(dp),     allocatable, dimension(:,:) :: M, S, S_b
    integer(i4b), allocatable, dimension(:)   :: ind

    alpha = p(3)

    ! Get the unmodulated covariance matrix
    allocate(S(n_h,n_h))
    S = 0.d0
    call compute_S_hat_qn(p(1:2), S)

    ! Construct harmonic space modulation matrix
    allocate(M(numcomp,numcomp))
    ind1 = 1
    M    = 0.d0
    do l1 = 0, lmax
       if (l1 < lmin_fit .or. l1 > lmax_fit) then
          ind1 = ind1 + 2*l1+1
          cycle
       end if
       do m1 = -l1, l1

          ind2 = 1
          do l2 = 0, lmax
             if (l2 < lmin_fit .or. l2 > lmax_fit) then
                ind2 = ind2 + 2*l2+1
                cycle
             end if
             do m2 = -l2, l2
                if (l1 == l2 .and. m1 == m2) then
                   M(ind1,ind2) = 1.d0
                else if (abs(l1-l2) == 1 .and. m1 == m2) then
                   lstar = min(l1,l2)
                   M(ind1,ind2) = alpha * sqrt((lstar+m1+1)*(lstar-m1+1)/(2*lstar+1)/(2*lstar+3))
                   M(ind2,ind1) = M(ind1,ind2)
                end if
                ind2 = ind2+1
             end do
          end do

          ind1 = ind1+1
       end do
    end do

    ! Multiply with M from the left
    n = size(indmap)/nmaps
    allocate(ind(n))
    ind = indmap(1:n)
    ind1 = 0
    allocate(S_b(n,n))
    do i = 1, nmaps
       do j = i, nmaps
          ind1 = ind1+1
          if (.not. inc_spec(ind1)) cycle
          call dgemm('N', 'N', n, n, n, 1.d0, M(ind,ind), n, &
               & S((i-1)*numcomp+ind,(j-1)*numcomp+ind), n, 0.d0, S_b, n)
          S_hat((i-1)*numcomp+ind,(j-1)*numcomp+ind) = S_b
          if (i /= j) S_hat((j-1)*numcomp+ind,(i-1)*numcomp+ind) = S_b
       end do
    end do

    ! Multiply with M^t from the right
    ind1 = 0
    do i = 1, nmaps
       do j = i, nmaps
          ind1 = ind1+1
          if (.not. inc_spec(ind1)) cycle
          call dgemm('N', 'T', n, n, n, 1.d0, S_hat((i-1)*numcomp+ind,(j-1)*numcomp+ind), n, &
               & M(ind,ind), n, 0.d0, S_b, n)
          S_hat((i-1)*numcomp+ind,(j-1)*numcomp+ind) = S_b
          if (i /= j) S_hat((j-1)*numcomp+ind,(i-1)*numcomp+ind) = S_b
       end do
    end do
    deallocate(M, S, S_b, ind)

  end subroutine compute_S_hat_power_asymmetry


  subroutine initialize_model_parameters(rng_handle, p, L_prop)
    implicit none

    type(planck_rng),                 intent(inout) :: rng_handle
    real(dp),         dimension(:),   intent(out)   :: p
    real(dp),         dimension(:,:), intent(out)   :: L_prop
    
    integer(i4b) :: i

    L_prop = 0.d0
    do i = 1, npar
       if (prior_gauss(i,2) == 0) then
          p(i)        = modulo(prior_gauss(i,1), par_modulo(i))
          L_prop(i,i) = 0.d0
       else
          p(i)        = prior_uni(i,1) + rand_uni(rng_handle) * (prior_uni(i,2)-prior_uni(i,1))
          p(i)        = modulo(p(i), par_modulo(i))
          L_prop(i,i) = par_rms(i)
       end if
    end do

    if (allocated(L_prop_in)) L_prop = L_prop_in

  end subroutine initialize_model_parameters

  subroutine get_sqrt_S(p, sqrt_S, ierr)
    implicit none

    real(dp),     dimension(:),   intent(in)  :: p
    real(dp),     dimension(:,:), intent(out) :: sqrt_S
    integer(i4b),                 intent(out) :: ierr
    
    integer(i4b) :: i, j, n, ind1
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:)         :: ind
    real(dp),     allocatable, dimension(:),   save :: p_prev
    real(dp),     allocatable, dimension(:,:), save :: R, sqrt_S_hat
    real(dp),     allocatable, dimension(:,:)       :: R_b, S_b, S_hat_b

    ierr = 0

    if (.not. allocated(p_prev)) then
       allocate(p_prev(npar), R(numcomp,numcomp), sqrt_S_hat(n_h,n_h))
       p_prev     = -1.d30
       R          = 0.d0
       sqrt_S_hat = 0.d0
    end if

    ! Check Euler parameters, and update rotation matrix R if necessary
    if (any(p(1:3) /= p_prev(1:3))) then
       p_prev(1:3) = p(1:3)
       call compute_rotation_matrix(p(1:3), R)
    end if

    ! Check model parameters, and update model matrix S if necessary
    if (any(p_prev(4:npar) /= p(4:npar))) then
       call compute_S_hat(p(4:npar), sqrt_S_hat)
       call cholesky_decompose_with_mask_dp(sqrt_S_hat, ierr=ierr)
       p_prev(4:npar) = p(4:npar)
    end if

    ! Multiply with rotation matrix
    if (all(prior_gauss(1:3,2) == 0.d0)) then
       ! No rotation requested
       sqrt_S = sqrt_S_hat 
    else
       n = size(indmap)/nmaps
       allocate(ind(n))
       ind = indmap(1:n)
       ind1 = 0
       allocate(S_b(n,n))
       do i = 1, nmaps
          do j = i, nmaps
             ind1 = ind1+1
             call dgemm('N', 'N', n, n, n, 1.d0, R(ind,ind), n, &
                  & sqrt_S_hat((i-1)*numcomp+ind,(j-1)*numcomp+ind), n, 0.d0, S_b, n)
             sqrt_S((i-1)*numcomp+ind,(j-1)*numcomp+ind) = S_b
          end do
       end do
       deallocate(S_b)

    end if

  end subroutine get_sqrt_S

  subroutine compute_rotation_matrix(euler, R)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: euler
    real(dp), dimension(1:,1:), intent(out) :: R

    real(dp)     :: psi, theta, phi, t1, t2
    integer(i4b) :: ind1, ind2, i, j, l, m, mp
    real(dp),     allocatable, dimension(:,:)   :: alms
    complex(dpc), allocatable, dimension(:,:,:) :: alms_cmplx

    phi   = euler(1)
    theta = euler(2)
    psi   = euler(3)

    allocate(alms(numcomp,1), alms_cmplx(1,0:lmax,0:lmax))

    R    = 0.d0
    do m = -lmax_fit, lmax_fit
       alms = 0.d0
       do l = abs(m), lmax_fit
          ind1 = l**2 + l + m + 1
          alms(ind1,1) = 1.d0
       end do
       call convert_real_to_complex_alms_dp(alms, alms_cmplx)
       call rotate_alm(lmax, alms_cmplx, psi, theta, phi)
       call convert_complex_to_real_alms(alms_cmplx,alms)
       do l = abs(m), lmax_fit
          ind1 = l**2 + l + m + 1
          do mp = -l, l
             ind2 = l**2 + l + mp + 1
             R(ind2,ind1) = alms(ind2,1)
          end do
       end do
    end do

    deallocate(alms, alms_cmplx)

  end subroutine compute_rotation_matrix

  subroutine output_prop_matrix(filename, scale, samples)
    implicit none

    real(dp),                         intent(in) :: scale
    character(len=*),                 intent(in) :: filename
    real(dp),         dimension(:,:), intent(in) :: samples

    integer(i4b) :: i, j, unit, n, nsamp
    real(dp), allocatable, dimension(:)   :: mu
    real(dp), allocatable, dimension(:,:) :: L

    unit  = comm_getlun()
    nsamp = size(samples,1)
    n     = size(samples,2)

    allocate(L(n,n), mu(n))
    L = 0.d0
    do i = 1, n
       mu = sum(samples(:,i))/nsamp
       do j = 1, i
          L(i,j) = sum((samples(:,i)-mu(i))*(samples(:,j)-mu(j))) / (nsamp-1)
       end do
    end do
    call cholesky_decompose_with_mask_dp(L)
    L = scale * L

    open(unit, file=trim(filename), form='unformatted')
    write(unit) n
    write(unit) L
    close(unit)

    deallocate(mu, L)

  end subroutine output_prop_matrix





end module comm_like_model_mod
