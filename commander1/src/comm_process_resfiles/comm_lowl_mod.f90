module comm_lowl_mod
  use comm_like_utils
  use comm_proc_utils
  implicit none

  ! *********************************************************************
  ! *  comm_lowl_mod -- An F90 module for computing a Gaussian          *
  ! *                   low-l likelihood by brute force                 *
  ! *                                                                   *
  ! *           H. K. Eriksen, E. Gjerløw (University of Oslo)          *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *                         I. K. Wehus (JPL)                         *   
  ! *                                                                   *
  ! *                                                                   *
  ! *  Please cite the following papers when using this code:           *
  ! *                                                                   *
  ! *  Gjerloew, Eriksen and Wehus 2014, ApJS, in preparation          *
  ! *                                                                   *
  ! *  History:                                                         *
  ! *      January 10th, 2014 -- First fully functional version         *
  ! *                                                                   *
  ! *********************************************************************

  ! ==========================================================================
  ! User routines:
  !
  !  ** subroutine comm_lowl_initialize_object(paramfile, handle)
  !
  !    Input parameters:
  !      paramfile   (character=*) :: Low-l parameter file
  ! 
  !    Output parameters:
  !      handle (i4b)  :: this code supports multiple objects (optional; default=1),
  !                       and the handle identifies specific data sets
  !
  ! 
  !  ** function comm_lowl_compute_lnL(cls, ierr, handle)
  !
  !    Input parameters:
  !      cls     (dp)  :: array containing at least cls(2:lmax) in units of l(l+1)/2pi
  !      ierr    (i4b) :: error flag
  !      handle  (i4b) :: handle selecting which data set to use (optional; default = 1)
  !    
  !  
  !  ** subroutine comm_lowl_deallocate_object(handle)
  ! 
  !    Input parameters:
  !      handle  (i4b) :: handle of object to deallocate (optional; default=remove all)
  !
  ! ==============================================================================

  integer(i4b), parameter :: MAX_N_LOWL = 5

  type comm_lowl_data
     logical(lgt) :: initialized=.false.
     integer(i4b) :: n, n_h, n_d, nmaps, llow, lhigh, lmax
     real(dp)     :: loglike_weight, cond_threshold
     real(dp)     :: lnL_recent, chisq_recent, red_chisq_recent
     real(dp)     :: cond_number_recent
     real(dp), allocatable, dimension(:,:)   :: d
     real(dp), allocatable, dimension(:,:)   :: cl_fid
     real(dp), allocatable, dimension(:,:)   :: N_cov, P_harm
     real(dp), allocatable, dimension(:,:,:) :: w
     real(dp), allocatable, dimension(:,:)   :: beam
  end type comm_lowl_data

  type(comm_lowl_data), dimension(MAX_N_LOWL) :: comm_lowl
  

contains

  ! Initialization routines
  subroutine comm_lowl_initialize_object(paramfile, handle)
    implicit none

    character(len=*),  intent(in)  :: paramfile
    integer(i4b),      intent(out), optional :: handle

    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep
    integer(i4b)       :: b, i, j, k, l, m, n, q, c1, c2, f, ind, d, n_h, n_d
    integer(i4b)       :: unit, numsamples, numchains, lmax_chain, id, bin(2), n_p, n_g, nmode
    integer(i4b)       :: col, p, nmaps, nspec, nsamp, lmax_cl, lmin_bin, lmax_bin, ncomp
    logical(lgt)       :: pattern(3,3), polarization, exist
    real(dp)           :: lnL, t1, t2
    character(len=512) :: line, s1, s2, datafile
    character(len=128) :: sigmafile, clfile
    real(dp),     allocatable, dimension(:,:)     :: cls
    integer(i4b), allocatable, dimension(:)       :: i2p
    real(dp),     allocatable, dimension(:,:,:,:) :: sigma
    real(dp),     allocatable, dimension(:,:,:)   :: sigma_1D
    real(dp),     allocatable, dimension(:,:)     :: S, sqrtS_P

    id = 1
    unit = comm_getlun()
    if (present(handle)) then
       do while (comm_lowl(id)%initialized .and. id < MAX_N_LOWL)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > MAX_N_LOWL) stop 'Error -- planck_br_mod: lowl id out of range'

    ! Initialize all distributions
    comm_lowl(id)%initialized = .true.
    call comm_get_parameter(paramfile, 'DATAFILE',         par_string=datafile)
    call comm_get_parameter(paramfile, 'FIDUCIAL_CL_FILE', par_string=clfile)
    call comm_get_parameter(paramfile, 'LOGLIKE_WEIGHT',   par_dp=comm_lowl(id)%loglike_weight)
    call comm_get_parameter(paramfile, 'LMIN',             par_int=comm_lowl(id)%llow)
    call comm_get_parameter(paramfile, 'LMAX',             par_int=comm_lowl(id)%lhigh)
    call comm_get_parameter(paramfile, 'CONDITION_NUMBER_THRESHOLD', &
         & par_dp=comm_lowl(id)%cond_threshold)

    ! Read data file
    inquire(file=trim(datafile), exist=exist)
    if (.not. exist) then
       write(*,*) 'Error -- low-l datafile = ', trim(datafile), ' does not exist'
       stop
    end if

    call read_lowl_datafile(datafile, comm_lowl(id)%n, comm_lowl(id)%n_d, comm_lowl(id)%n_h, &
         & comm_lowl(id)%lmax, comm_lowl(id)%nmaps, comm_lowl(id)%d, comm_lowl(id)%N_cov, &
         & comm_lowl(id)%beam, comm_lowl(id)%P_harm, comm_lowl(id)%w)

    !comm_lowl(id)%N_cov = comm_lowl(id)%N_cov * 0.75d0**2

    ncomp = (comm_lowl(id)%lmax + 1) ** 2
    nmode = comm_lowl(id)%n
    nmaps = comm_lowl(id)%nmaps
    n_h   = comm_lowl(id)%n_h

    if (comm_lowl(id)%lhigh > comm_lowl(id)%lmax) then
       write(*,*) 'comm_lowl_mod -- Error: Requested LMAX greater than maximum '
       write(*,*) '                        multipole in datafile'
       stop
    end if

    ! Read in fiducial spectrum, to be conditioned upon outside range of interest
    call read_fiducial_spectrum(clfile, comm_lowl(id)%cl_fid)

    ! Add basis vectors with fixed C_l's directly into noise covariance
    allocate(sqrtS_P(n_h,nmode), S(nmaps,nmaps))
    sqrtS_P = 0.d0
    ind     = 5
    do l = 2, comm_lowl(id)%lmax
       if (l >= comm_lowl(id)%llow .and. l <= comm_lowl(id)%lhigh) then
          ind = ind + 2*l+1
          cycle
       end if
       call cl2s(comm_lowl(id)%cl_fid(l,:), S) ! Fix to fiducial
       S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       call cholesky_decompose_with_mask_dp(S)
       do m = -l, l
          call dgemm('N', 'N', nmaps, nmode, nmaps, 1.d0, S, nmaps, &
               & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
               & nmaps, 0.d0, sqrtS_P(ind:n_h:ncomp,:), nmaps)
          ind = ind+1
       end do
    end do
    call dsyrk('L','T',nmode,n_h,1.d0,sqrtS_P,n_h,1.d0,comm_lowl(id)%N_cov,nmode)
    deallocate(sqrtS_P, S)

    ! Symmetrize matrix
    do i = 1, nmode
       do j = i+1, nmode
          comm_lowl(id)%N_cov(i,j) = comm_lowl(id)%N_cov(j,i)
       end do
    end do

  end subroutine comm_lowl_initialize_object

  subroutine comm_lowl_initialize_object_direct(lmin_like, lmax_like, lmax, n_d, &
       & n_h, nmode, nmaps, d, N_cov, w, beam, P_harm, cl_fid, loglike_weight, &
       & cond_threshold, handle)
    implicit none

    integer(i4b),       intent(in)      :: lmin_like, lmax_like, lmax, n_d
    integer(i4b),       intent(in)      :: n_h, nmode, nmaps
    real(dp), intent(in)        :: loglike_weight, cond_threshold
    real(dp), dimension(:, :), intent(in)  :: d, N_cov, P_harm
    real(dp), dimension(0:, :), intent(in)  :: beam, cl_fid
    real(dp), dimension(0:, :, :), intent(in)  :: w
    integer(i4b),      intent(out), optional :: handle

    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep
    integer(i4b)       :: b, i, j, k, l, m, n, q, c1, c2, f, ind
    integer(i4b)       :: unit, numsamples, numchains, lmax_chain, id, bin(2), n_p, n_g
    integer(i4b)       :: col, p, nspec, nsamp, lmax_cl, lmin_bin, lmax_bin, ncomp
    logical(lgt)       :: pattern(3,3), polarization, exist
    real(dp)           :: lnL, t1, t2
    character(len=512) :: line, s1, s2, datafile
    character(len=128) :: sigmafile, clfile
    real(dp),     allocatable, dimension(:,:)     :: cls
    integer(i4b), allocatable, dimension(:)       :: i2p
    real(dp),     allocatable, dimension(:,:,:,:) :: sigma
    real(dp),     allocatable, dimension(:,:,:)   :: sigma_1D
    real(dp),     allocatable, dimension(:,:)     :: S, sqrtS_P

    id = 1
    unit = comm_getlun()
    if (present(handle)) then
       do while (comm_lowl(id)%initialized .and. id < MAX_N_LOWL)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > MAX_N_LOWL) stop 'Error -- planck_br_mod: lowl id out of range'

    ! Initialize all distributions
    comm_lowl(id)%initialized = .true.

   comm_lowl(id)%loglike_weight = loglike_weight
   comm_lowl(id)%llow = lmin_like
   comm_lowl(id)%lhigh = lmax_like
   comm_lowl(id)%cond_threshold = cond_threshold
   comm_lowl(id)%lmax = lmax
   ncomp = (comm_lowl(id)%lmax+1) ** 2
   comm_lowl(id)%n = nmode
   comm_lowl(id)%nmaps = nmaps
   comm_lowl(id)%n_h = n_h
   comm_lowl(id)%n_d = n_d
   allocate(comm_lowl(id)%w(0:lmax,nmaps,nmaps))
   allocate(comm_lowl(id)%beam(0:lmax,nmaps))
   allocate(comm_lowl(id)%d(nmode,n_d), comm_lowl(id)%N_cov(nmode,nmode))
   allocate(comm_lowl(id)%P_harm(n_h,nmode))
   comm_lowl(id)%w = w(0:lmax, 1:nmaps, 1:nmaps)
   comm_lowl(id)%beam = beam(0:lmax, 1:nmaps)
   comm_lowl(id)%d = d
   comm_lowl(id)%N_cov = N_cov
   comm_lowl(id)%P_harm = P_harm

    if (comm_lowl(id)%lhigh > comm_lowl(id)%lmax) then
       write(*,*) 'comm_lowl_mod -- Error: Requested LMAX greater than maximum '
       write(*,*) '                        multipole in datafile'
       stop
    end if

   allocate(comm_lowl(id)%cl_fid(0:comm_lowl(id)%lmax,6))
   comm_lowl(id)%cl_fid = cl_fid(0:comm_lowl(id)%lmax, :)

    ! Add basis vectors with fixed C_l's directly into noise covariance
    allocate(sqrtS_P(n_h,nmode), S(nmaps,nmaps))
    sqrtS_P = 0.d0
    ind     = 5
    do l = 2, comm_lowl(id)%lmax
       if (l >= comm_lowl(id)%llow .and. l <= comm_lowl(id)%lhigh) then
          ind = ind + 2*l+1
          cycle
       end if
       call cl2s(comm_lowl(id)%cl_fid(l,:), S) ! Fix to fiducial
       S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       call cholesky_decompose_with_mask_dp(S)
       do m = -l, l
          call dgemm('N', 'N', nmaps, nmode, nmaps, 1.d0, S, nmaps, &
               & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
               & nmaps, 0.d0, sqrtS_P(ind:n_h:ncomp,:), nmaps)
          ind = ind+1
       end do
    end do
    call dsyrk('L','T',nmode,n_h,1.d0,sqrtS_P,n_h,1.d0,comm_lowl(id)%N_cov,nmode)
    deallocate(sqrtS_P, S)

    ! Symmetrize matrix
    do i = 1, nmode
       do j = i+1, nmode
          comm_lowl(id)%N_cov(i,j) = comm_lowl(id)%N_cov(j,i)
       end do
    end do

  end subroutine comm_lowl_initialize_object_direct


  subroutine comm_lowl_deallocate_object(handle)
    implicit none
    
    integer(i4b), optional :: handle

    integer(i4b) :: i, j, k, id, id_min, id_max

    id_min = 1; id_max = MAX_N_LOWL; 
    if (present(handle)) then 
       id_min = handle
       id_max = handle
    end if

    do id = id_min, id_max
       if (comm_lowl(id)%initialized) then
          deallocate(comm_lowl(id)%d, comm_lowl(id)%cl_fid, comm_lowl(id)%N_cov)
          deallocate(comm_lowl(id)%P_harm, comm_lowl(id)%w, comm_lowl(id)%beam)
       end if
      comm_lowl(id)%initialized = .false.
    end do

  end subroutine comm_lowl_deallocate_object


  ! Base computation routine
  function comm_lowl_compute_lnL(cls, sqrt_S, chisq, red_chisq, handle, enforce_pos_def, &
       & lnL_multi, cond_number, ierr)
    implicit none

    real(dp),     dimension(0:,1:), intent(in),  optional :: cls
    integer(i4b),                   intent(out), optional :: ierr
    integer(i4b),                   intent(in),  optional :: handle
    logical(lgt),                   intent(in),  optional :: enforce_pos_def
    real(dp),     dimension(1:,1:), intent(in),  optional :: sqrt_S
    real(dp),                       intent(out), optional :: chisq, red_chisq, cond_number
    real(dp),     dimension(1:),    intent(out), optional :: lnL_multi
    real(dp)                                              :: comm_lowl_compute_lnL

    integer(i4b) :: i, j, l, m, n, id, k, stat, n_d
    logical(lgt) :: posdef, enf_pos_def
    real(dp)     :: chi2, logdet, t1, t2, L_max, L_min, cond
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: C, invC_d, V, map

    if (present(ierr)) ierr = 0
    id = 1; if (present(handle)) id = handle
    enf_pos_def = .true.; if (present(enforce_pos_def)) enf_pos_def = enforce_pos_def

    ! Check that likelihood structure is initialized
    if (.not. comm_lowl(id)%initialized) then
       write(*,*) 'Error -- comm_lowl_mod: Requested handle ', id, ' is not initialized'
       stop
    end if

    ! Compute covariance matrix
    n   = comm_lowl(id)%n
    n_d = comm_lowl(id)%n_d
    allocate(C(n,n))
    if (present(sqrt_S)) then
       call comm_lowl_getC(enf_pos_def, C, stat, sqrt_S=sqrt_S)
    else
       call comm_lowl_getC(enf_pos_def, C, stat, cls=cls)
    end if

    if (stat /= 0) then
       if (present(ierr)) ierr = 1
       comm_lowl_compute_lnL = -1.d30
       deallocate(C)
       return
    end if

    ! Cholesky decompose matrix
    call dpotrf('L', n, C, n, stat)
    if (stat /= 0) then
       if (present(ierr)) ierr = ierr + 1
       comm_lowl_compute_lnL = -1.d30
       deallocate(C)
       return
    end if

    ! Compute log-determinant
    logdet =  0.d0
    L_max  = -1.d30
    L_min  =  1.d30
    do i = 1, n
       logdet = logdet + 2.d0 * log(C(i,i))
       L_max  = max(L_max, C(i,i))
       L_min  = min(L_min, C(i,i))
    end do
    cond = (L_max/L_min)**2
    if (present(cond_number)) cond_number = cond

    ! Compute chi-square term
    allocate(invC_d(n,n_d))
    invC_d = comm_lowl(id)%d
    call dpotrs('L', n, n_d, C, n, invC_d, n, stat)

    ! Return log-like value
    if (stat == 0 .and. cond < comm_lowl(id)%cond_threshold) then
       do i = 1, n_d
          chi2 = sum(comm_lowl(id)%d(:,i)*invC_d(:,i))
          if (present(lnL_multi)) lnL_multi(i) = -0.5d0 * (chi2 + logdet)
          if (i == 1) then
             if (present(chisq))     chisq     = chi2
             if (present(red_chisq)) red_chisq = chi2 / n
             comm_lowl_compute_lnL = -0.5d0 * (chi2 + logdet)
          end if
       end do
    else
       comm_lowl_compute_lnL = -1.d30
       if (present(ierr))      ierr      = 1
       if (present(lnL_multi)) lnL_multi = -1.d30
    end if

    ! Update data structures for quick look-up 
    comm_lowl(id)%lnL_recent         = comm_lowl_compute_lnL
    comm_lowl(id)%chisq_recent       = sum(comm_lowl(id)%d(:,1)*invC_d(:,1))
    comm_lowl(id)%red_chisq_recent   = comm_lowl(id)%chisq_recent / n
    comm_lowl(id)%cond_number_recent = cond

    deallocate(C, invC_d)

  end function comm_lowl_compute_lnL

  subroutine comm_lowl_getC(enforce_pos_def, C, ierr, handle, cls, sqrt_S)
    implicit none

    real(dp),              dimension(0:,1:), intent(in), optional :: cls
    real(dp),              dimension(1:,1:), intent(in), optional :: sqrt_S
    integer(i4b),                            intent(in), optional :: handle
    logical(lgt),                            intent(in)           :: enforce_pos_def
    real(dp), allocatable, dimension(:,:),   intent(out)          :: C
    integer(i4b),                            intent(out)          :: ierr

    real(dp)     :: t1, t2, var
    integer(i4b) :: i, j, l, m, ind, ind1, ind2, id, stat
    integer(i4b) :: n_h, n, nmaps, lmax, ncomp, ind_min, ind_max, nmode
    real(dp), allocatable, dimension(:,:) :: S, sqrtS_P, V, P_sub, map

    id = 1; if (present(handle)) id = handle
    n_h     = comm_lowl(id)%n_h
    n       = comm_lowl(id)%n
    nmaps   = comm_lowl(id)%nmaps
    lmax    = comm_lowl(id)%lmax
    ncomp   = comm_lowl(id)%n_h / comm_lowl(id)%nmaps
    ind_min = comm_lowl(id)%llow**2+1
    ind_max = (comm_lowl(id)%lhigh+1)**2
    nmode   = (ind_max-ind_min+1)*nmaps
    ierr    = 0

    allocate(C(n,n), sqrtS_P(nmode,n), S(nmaps,nmaps))

    if (present(sqrt_S)) then
       ! Use user-supplied harmonic space covariance

       ! Extract vectors to be multiplied with sqrt_S, and multiply with beam
       allocate(P_sub(nmode,n))
       ind2 = 1
       do i = 1, nmaps
          ind1 = (i-1)*ncomp + comm_lowl(id)%llow**2 + 1
          do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
             do m = -l, l
                P_sub(ind2,:) = comm_lowl(id)%beam(l,i)*comm_lowl(id)%P_harm(ind1,:)
                ind2 = ind2+1
                ind1 = ind1+1
             end do
          end do
       end do

       ! Multiply basis vectors with sqrt(S), ie., Cholesky factor of S
       call dgemm('T', 'N', nmode, n, nmode, 1.d0, sqrt_S, nmode, &
            & P_sub, nmode, 0.d0, sqrtS_P, nmode)
       call dsyrk('L','T',n,nmode,1.d0,sqrtS_P,nmode,0.d0,C,n)
       do i = 1, n
          do j = i+1, n
             C(i,j) = C(j,i)
          end do
       end do
       deallocate(P_sub)

    else if (enforce_pos_def) then
       ! Use user-supplied power spectrum for S

       ! Check that spectrum is positive definite
       if (.not. comm_cls_posdef(cls(2:lmax,:))) then
          ierr = 1
          C    = 0.d0
          deallocate(sqrtS_P, S)
          return
       end if

       ! Compute signal covariance matrix
       sqrtS_P        = 0.d0
       ind            = comm_lowl(id)%llow**2+1
       ind2           = 1
       do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
          call cl2s(cls(l,:), S)      
          S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
          call cholesky_decompose_with_mask_dp(S,ierr=stat)
          if (stat /= 0) then
             ierr = 1       
             C    = 0.d0
             deallocate(sqrtS_P, S)
             return
          end if
          do m = -l, l
             call dgemm('N', 'N', nmaps, n, nmaps, 1.d0, S, nmaps, &
                  & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
                  & nmaps, 0.d0, sqrtS_P(ind2:ind2+nmaps-1,:), nmaps)
             ind  = ind  + 1
             ind2 = ind2 + nmaps
          end do
       end do
       call dsyrk('L','T',n,nmode,1.d0,sqrtS_P,nmode,0.d0,C,n)
       do i = 1, n
          do j = i+1, n
             C(i,j) = C(j,i)
          end do
       end do
       
    else 
       ! Use user-supplied power spectrum for S; do not enforce positive definiteness

       ! Compute signal covariance matrix
       allocate(P_sub(nmode,n))
       sqrtS_P        = 0.d0
       ind            = comm_lowl(id)%llow**2+1
       ind2           = 1
       do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
          call cl2s(cls(l,:), S)      
          S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
          do m = -l, l
             call dgemm('N', 'N', nmaps, n, nmaps, 1.d0, S, nmaps, &
                  & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
                  & nmaps, 0.d0, sqrtS_P(ind2:ind2+nmaps-1,:), nmaps)
             P_sub(ind2:ind2+nmaps-1,:) = comm_lowl(id)%P_harm(ind:n_h:ncomp,:)
             ind  = ind  + 1
             ind2 = ind2 + nmaps
          end do
       end do
       call dgemm('T', 'N', n, n, nmode, 1.d0, P_sub, nmode, &
            & sqrtS_P, nmode, 0.d0, C, n)

       deallocate(P_sub)
    end if

    ! Add noise
    C = C + comm_lowl(id)%N_cov

    ! Validate against correlation function; only works for pixel space basis
    if (.false.) then
       var = 0.d0
       do l = 2, lmax
          var = var + (2*l+1)/(4.d0*pi) * cls(l,1) * &
               & 2.d0*pi/real(l*(l+1),dp) * comm_lowl(id)%beam(l,1)**2
       end do
       write(*,*) C(1,1), var
       stop
    end if

    deallocate(sqrtS_P, S)

  end subroutine comm_lowl_getC

  subroutine comm_lowl_get_recent_value(handle, cond_number, chisq, red_chisq, lnL)
    implicit none

    integer(i4b), intent(in),  optional :: handle
    real(dp),     intent(out), optional :: cond_number, chisq, red_chisq, lnL

    integer(i4b) :: id

    id = 1; if (present(handle)) id = handle
    if (present(cond_number)) cond_number = comm_lowl(id)%cond_number_recent
    if (present(chisq))       chisq       = comm_lowl(id)%chisq_recent
    if (present(red_chisq))   red_chisq   = comm_lowl(id)%red_chisq_recent
    if (present(lnL))         lnL         = comm_lowl(id)%lnL_recent

  end subroutine comm_lowl_get_recent_value

end module comm_lowl_mod
