module comm_diffuse_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  use comm_data_mod
  use comm_F_int_mod
  use comm_Cl_mod
  use math_tools
  use comm_cr_utils
  implicit none

  private
  public comm_diffuse_comp, precondDiff, precond_diff_comps, add_to_npre
  
  !**************************************************
  !            Diffuse component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nx, x0
     logical(lgt)       :: pol, output_mixmat
     integer(i4b)       :: lmax_amp, lmax_ind, lpiv
     real(dp)           :: cg_scale
     real(dp), allocatable, dimension(:,:) :: cls
     real(dp), allocatable, dimension(:,:) :: F_mean

     class(comm_map),               pointer     :: mask
     class(comm_map),               pointer     :: x      ! Spatial parameters
     class(comm_map),               pointer     :: mu     ! Spatial prior mean
     class(comm_B),                 pointer     :: B_out  ! Output beam
     class(comm_Cl),                pointer     :: Cl     ! Power spectrum
     class(map_ptr),  dimension(:), allocatable :: theta  ! Spectral parameters
     type(map_ptr),   dimension(:), allocatable :: F      ! Mixing matrix
     logical(lgt),    dimension(:), allocatable :: F_null ! Don't allocate space for null mixmat's
     type(F_int_ptr), dimension(:), allocatable :: F_int  ! SED integrator
   contains
     procedure :: initDiffuse
     procedure :: updateMixmat
!!$     procedure :: dumpHDF  => dumpDiffuseToHDF
     procedure :: getBand     => evalDiffuseBand
     procedure :: projectBand => projectDiffuseBand
     procedure :: dumpFITS    => dumpDiffuseToFITS
  end type comm_diffuse_comp

  type diff_ptr
     class(comm_diffuse_comp), pointer :: p
  end type diff_ptr
  
  type invM_lm
     integer(i4b)                              :: n
     integer(i4b), allocatable, dimension(:)   :: ind, comp2ind
     real(dp),     allocatable, dimension(:,:) :: M0, M
  end type invM_lm

  type precondDiff
     real(dp)     :: Nscale
     type(invM_lm), allocatable, dimension(:,:)      :: invM_ ! (0:nalm-1,nmaps)
   contains
     procedure :: update  => updateDiffPrecond
  end type precondDiff

  interface precondDiff
     procedure constPreDiff
  end interface precondDiff
  
  integer(i4b) :: npre      =  0
  integer(i4b) :: lmax_pre  = -1
  integer(i4b) :: nside_pre = 1000000
  integer(i4b) :: nmaps_pre = -1
  logical(lgt) :: output_cg_eigenvals
  character(len=512) :: outdir
  class(comm_mapinfo), pointer                   :: info_pre
  class(diff_ptr),     allocatable, dimension(:) :: diffComps
  
contains

  subroutine initDiffuse(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info
    
    call self%initComp(cpar, id, id_abs)

    ! Initialize variables specific to diffuse source type
    self%pol           = cpar%cs_polarization(id_abs)
    self%nside         = cpar%cs_nside(id_abs)
    self%lmax_amp      = cpar%cs_lmax_amp(id_abs)
    self%lmax_ind      = cpar%cs_lmax_ind(id_abs)
    self%cltype        = cpar%cs_cltype(id_abs)
    self%cg_scale      = cpar%cs_cg_scale(id_abs)
    self%nmaps         = 1; if (self%pol) self%nmaps = 3
    self%output_mixmat = cpar%output_mixmat
    output_cg_eigenvals = cpar%output_cg_eigenvals
    outdir              = cpar%outdir
    info               => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)

    ! Diffuse preconditioner variables
    npre      = npre+1
    lmax_pre  = max(lmax_pre,  self%lmax_amp)
    nside_pre = min(nside_pre, self%nside)
    nmaps_pre = max(nmaps_pre, self%nmaps)

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id_abs)) == 'zero' .or. trim(cpar%cs_input_amp(id_abs)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       self%x => comm_map(info, cpar%cs_input_amp(id_abs))
       call self%x%YtW
    end if
    self%ncr = size(self%x%alm)

    ! Initialize prior mean
    if (trim(cpar%cs_prior_amp(id_abs)) /= 'none') then
       self%mu => comm_map(info, cpar%cs_prior_amp(id_abs))
       call self%mu%YtW
    end if

    ! Allocate mixing matrix
    call update_status(status, "init_pre_mix")
    allocate(self%F(numband), self%F_mean(numband,self%nmaps), self%F_null(numband))
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & self%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
       self%F(i)%p    => comm_map(info)
       self%F_null(i) =  .false.
    end do
    call update_status(status, "init_postmix")

    ! Initialize output beam
    self%B_out => comm_B_bl(cpar, self%x%info, 0, fwhm=cpar%cs_fwhm(id_abs))

    ! Initialize power spectrum
    self%Cl => comm_Cl(cpar, self%x%info, id, id_abs)
    
  end subroutine initDiffuse

  function constPreDiff(comm, Nscale)
    implicit none
    integer(i4b),                intent(in) :: comm
    real(dp),                    intent(in) :: Nscale
    class(precondDiff), pointer             :: constPreDiff

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c
    class(comm_diffuse_comp), pointer :: p1, p2
    real(dp),     allocatable, dimension(:,:) :: mat

    allocate(constPreDiff)

    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%next()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if

    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(constPreDiff%invM_(0:info_pre%nalm-1,info_pre%nmaps))
    !!$OMP PARALLEL PRIVATE(mat, ind, j, i1, l, m, q, i2, k1, p1, k2, n)
    allocate(mat(npre,npre), ind(npre))
    do j = 1, info_pre%nmaps
       !!$OMP DO SCHEDULE(guided)
       do i1 = 0, info_pre%nalm-1
          call info_pre%i2lm(i1, l, m)
          mat = 0.d0
          do q = 1, numband
             call data(q)%info%lm2i(l,m,i2)
             if (i2 == -1) cycle
             do k1 = 1, npre
                p1 => diffComps(k1)%p
                if (l > p1%lmax_amp) cycle
                do k2 = 1, npre
                   p2 => diffComps(k2)%p
                   if (l > p2%lmax_amp) cycle
                   mat(k1,k2) = mat(k1,k2) + &
                        & data(q)%N%invN_diag%alm(i2,j) * &  ! invN_{lm,lm}
                        & data(q)%B%b_l(l,j)**2 * &          ! b_l^2
                        & p1%F_mean(q,j) * p2%F_mean(q,j)    ! F(c1)*F(c2)
                end do
             end do
          end do

          n = 0
          allocate(constPreDiff%invM_(i1,j)%comp2ind(npre))
          constPreDiff%invM_(i1,j)%comp2ind = -1
          do k1 = 1, npre
             if (mat(k1,k1) > 0.d0) then
                n = n+1
                ind(n) = k1
                constPreDiff%invM_(i1,j)%comp2ind(k1) = n
             end if
          end do
          constPreDiff%invM_(i1,j)%n = n
          allocate(constPreDiff%invM_(i1,j)%ind(n))
          allocate(constPreDiff%invM_(i1,j)%M0(n,n), constPreDiff%invM_(i1,j)%M(n,n))
          constPreDiff%invM_(i1,j)%ind = ind(1:n)
          constPreDiff%invM_(i1,j)%M0   = mat(ind(1:n),ind(1:n))

       end do
       !!$OMP END DO
    end do
    deallocate(ind, mat)
    !!$OMP END PARALLEL
    call wall_time(t2)
    !if (info_pre%myid == 0) write(*,*) 'constPreDiff = ', t2-t1

!!$    if (info_pre%myid == 0) then
!!$       write(*,*) 'wall time = ', t2-t1
!!$       do i = 0, 10
!!$          write(*,*) info_pre%myid, i, sum(abs(constPreDiff%invM_(i,1)%M0))
!!$       end do
!!$    end if
!!$    call mpi_finalize(i)
!!$    stop

  end function constPreDiff

  subroutine updateDiffPrecond(self)
    implicit none
    class(precondDiff), intent(inout) :: self

    integer(i4b) :: i, j, k, k1, k2, l, m, n, ierr, unit, p, q
    real(dp)     :: W_ref, t1, t2
    logical(lgt), save :: first_call = .true.
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: W
    real(dp),     allocatable, dimension(:) :: cond, W_min
    real(dp),     allocatable, dimension(:,:) :: alm

    ! Initialize current preconditioner to F^t * B^t * invN * B * F
    !self%invM    = self%invM0
    do j = 1, info_pre%nmaps
       do i = 0, info_pre%nalm-1
          self%invM_(i,j)%M = self%invM_(i,j)%M0
       end do
    end do

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W0 = ', real(W,sp)
!!$       deallocate(W)
!!$       write(*,*) 
!!$       do i = 1, 3
!!$          write(*,*) real(self%invM_(0,1)%M(i,:),sp)
!!$       end do
!!$       write(*,*)
!!$    end if

    ! Right-multiply with sqrt(Cl)
    call wall_time(t1)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                p = self%invM_(i,j)%comp2ind(k1)
                q = self%invM_(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = self%invM_(i,j)%M(q,p)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre)
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                p = self%invM_(i,j)%comp2ind(k1)
                q = self%invM_(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) self%invM_(i,j)%M(q,p) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do

    ! Left-multiply with sqrt(Cl)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                p = self%invM_(i,j)%comp2ind(k1)
                q = self%invM_(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = self%invM_(i,j)%M(p,q)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre)
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                p = self%invM_(i,j)%comp2ind(k1)
                q = self%invM_(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) self%invM_(i,j)%M(p,q) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do
    !call wall_time(t2)
    !write(*,*) 'sqrtS = ', t2-t1

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W1 = ', real(W,sp)
!!$       deallocate(W)
!!$       write(*,*) 
!!$       do i = 1, 3
!!$          write(*,*) real(self%invM_(0,1)%M(i,:),sp)
!!$       end do
!!$       write(*,*)
!!$    end if


    ! Add unity 
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       !$OMP PARALLEL PRIVATE(i,l,m,j,p)
       !$OMP DO SCHEDULE(guided)
       do i = 0, info_pre%nalm-1
          call info_pre%i2lm(i, l, m)
          if (l <= diffComps(k1)%p%lmax_amp) then
             do j = 1, info_pre%nmaps
                p = self%invM_(i,j)%comp2ind(k1)
                self%invM_(i,j)%M(p,p) = self%invM_(i,j)%M(p,p) + 1.d0
             end do
          end if
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end do

!!$    if (info_pre%myid == 0) then
!!$       n = self%invM_(0,1)%n
!!$       allocate(W(n))
!!$       call get_eigenvalues(self%invM_(0,1)%M, W)
!!$       write(*,*) 'W2 = ', real(W,sp)
!!$       deallocate(W)
!!$    end if

!!$    call mpi_finalize(i)
!!$    stop

    ! Print out worst condition number in first call
    call wall_time(t1)
    if (first_call .and. output_cg_eigenvals) then
       allocate(cond(0:lmax_pre), W_min(0:lmax_pre))
       cond  = 0.d0
       W_min = 1.d30
       do j = 1, nmaps_pre
          do i = 0, info_pre%nalm-1
             call info_pre%i2lm(i, l, m)
             n = self%invM_(i,j)%n
             if (n > 0) then
                allocate(W(n), ind(n))
                call get_eigenvalues(self%invM_(i,j)%M, W)
                W_ref    = minval(abs(W))
                cond(l)  = max(cond(l), maxval(abs(W/W_ref)))
                W_min(l) = min(W_min(l), minval(W))
                deallocate(W, ind)
             end if
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, cond,  lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MAX, info_pre%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, W_min, lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MIN, info_pre%comm, ierr)
       if (info_pre%myid == 0) then
          write(*,*) 'Precond -- largest condition number = ', real(maxval(cond),sp)
          write(*,*) 'Precond -- smallest eigenvalue      = ', real(minval(W_min),sp)
          unit = getlun()
          open(unit, file=trim(outdir)//'/precond_eigenvals.dat', recl=1024)
          do l = 0, lmax_pre
             write(unit,fmt='(i8,2e16.8)') l, cond(l), W_min(l)
          end do
          close(unit)
       end if

       first_call = .false.
       deallocate(cond, W_min)
    end if
    call wall_time(t2)
    !write(*,*) 'eigen = ', t2-t1

    ! Invert matrix
    call wall_time(t1)
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          call invert_matrix_with_mask(self%invM_(i,j)%M)
       end do
    end do
    call wall_time(t2)
    !write(*,*) 'invert = ', t2-t1
       
  end subroutine updateDiffPrecond
    
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  subroutine updateMixmat(self, theta)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta

    integer(i4b) :: i, j, k, n, nmaps, ierr
    real(dp),        allocatable, dimension(:,:) :: theta_p
    real(dp),        allocatable, dimension(:)   :: nu, s, buffer
    class(comm_mapinfo),          pointer        :: info
    class(comm_map),              pointer        :: t, t0
    
    ! Copy over alms from input structure, and compute pixel-space parameter maps
    if (present(theta)) then
       do i = 1, self%npar
          self%theta(i)%p%alm = theta(i)%alm
       end do
    end if
    
    ! Compute mixing matrix
    allocate(theta_p(self%npar,self%nmaps))
    do i = 1, numband

       ! Don't update null mixing matrices
       if (self%F_null(i)) cycle

       ! Compute spectral parameters at the correct resolution for this channel
       if (self%npar > 0) then
          info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, self%theta(1)%p%info%lmax, &
               & data(i)%info%nmaps, data(i)%info%pol)
          t    => comm_map(info)
          nmaps            = min(data(i)%info%nmaps, self%theta(1)%p%info%nmaps)
          t%alm(:,1:nmaps) = self%theta(1)%p%alm(:,1:nmaps)
          call t%Y
          nullify(info)
          do j = 2, self%npar
             info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                  & self%theta(j)%p%info%lmax, data(i)%info%nmaps, data(i)%info%pol)
             t0    => comm_map(info)
             nmaps            = min(data(i)%info%nmaps, self%theta(j)%p%info%nmaps)
             t0%alm(:,1:nmaps) = self%theta(j)%p%alm(:,1:nmaps)
             call t0%Y
             call t%add(t0)
             nullify(info)
          end do
       end if
       
       ! Loop over all pixels, computing mixing matrix for each
       do j = 0, self%F(i)%p%info%np-1
          if (self%npar > 0) then
             ! Collect all parameters
             t0 => t
             theta_p(1,:) = t0%map(j,:)
             do k = 2, self%npar
                t0 => t0%next()
                theta_p(k,:) = t0%map(j,:)
             end do

             ! Check polarization type
             if (self%nmaps == 3) then
                do k = 1, self%npar
                   if (self%poltype(k) < 2) theta_p(k,2) = theta_p(k,1)
                   if (self%poltype(k) < 3) theta_p(k,3) = theta_p(k,2)
                end do
             end if
          end if

          ! Temperature
          self%F(i)%p%map(j,1) = self%F_int(i)%p%eval(theta_p(:,1)) * data(i)%gain * self%cg_scale

          ! Polarization
          if (self%nmaps == 3) then
             ! Stokes Q
             if (all(self%poltype < 2)) then
                self%F(i)%p%map(j,2) = self%F(i)%p%map(j,1) * data(i)%gain * self%cg_scale
             else
                self%F(i)%p%map(j,2) = self%F_int(i)%p%eval(theta_p(:,2)) * data(i)%gain * self%cg_scale
             end if
       
             ! Stokes U
             if (all(self%poltype < 3)) then
                self%F(i)%p%map(j,3) = self%F(i)%p%map(j,2) * data(i)%gain * self%cg_scale
             else
                self%F(i)%p%map(j,3) = self%F_int(i)%p%eval(theta_p(:,3)) * data(i)%gain * self%cg_scale
             end if
          end if
          
       end do

       ! Compute mixing matrix average; for preconditioning only
       do j = 1, self%nmaps
          self%F_mean(i,j) = sum(self%F(i)%p%map(:,j))
       end do
       call mpi_allreduce(MPI_IN_PLACE, self%F_mean(i,:), self%nmaps, &
            & MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       self%F_mean(i,:) = self%F_mean(i,:) / self%F(i)%p%info%npix

    end do
    nullify(t, t0)
    deallocate(theta_p)

  end subroutine updateMixmat

  function evalDiffuseBand(self, band, amp_in, pix, alm_out)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    real(dp),        dimension(:,:), allocatable                        :: evalDiffuseBand

    integer(i4b) :: i, j, np, nmaps, lmax, nmaps_comp
    logical(lgt) :: alm_out_
    real(dp)     :: t1, t2
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: m

    alm_out_ = .false.; if (present(alm_out)) alm_out_ = alm_out

    if (self%F_null(band)) then
       if (alm_out_) then
          if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
       else
          if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       end if
       evalDiffuseBand = 0.d0
       return
    end if

    ! Initialize amplitude map
    nmaps =  min(data(band)%info%nmaps, self%nmaps)
    info  => comm_mapinfo(data(band)%info%comm, data(band)%info%nside, self%lmax_amp, &
         & data(band)%info%nmaps, data(band)%info%pol)
    m     => comm_map(info)
    if (present(amp_in)) then
       m%alm(:,1:nmaps) = amp_in
    else
       m%alm(:,1:nmaps) = self%x%alm(:,1:nmaps)
    end if
    if (self%lmax_amp > data(band)%map%info%lmax) then
       ! Nullify elements above band-specific lmax to avoid aliasing during projection
       do i = 0, m%info%nalm-1
          if (m%info%lm(1,i) > data(band)%info%lmax) m%alm(i,:) = 0.d0
       end do
    end if

    ! Scale to correct frequency through multiplication with mixing matrix
    if (self%lmax_ind == 0) then
       do i = 1, m%info%nmaps
          m%alm(:,i) = m%alm(:,i) * self%F_mean(band,i)
       end do
    else
       call m%Y()
       m%map = m%map * self%F(band)%p%map
       call m%YtW()
    end if
       
    ! Convolve with band-specific beam
    call data(band)%B%conv(alm_in=(self%lmax_ind==0), alm_out=.false., trans=.false., map=m)
    if (.not. alm_out_) call m%Y()

    ! Return correct data product
    if (alm_out_) then
       if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
       evalDiffuseBand = m%alm
    else
       if (.not. allocated(evalDiffuseBand)) allocate(evalDiffuseBand(0:data(band)%info%np-1,data(band)%info%nmaps))
       evalDiffuseBand = m%map
    end if
       

    ! Clean up
    deallocate(m, info)

  end function evalDiffuseBand

  ! Return component projected from map
  function projectDiffuseBand(self, band, map, alm_in)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    real(dp),        dimension(:,:), allocatable                        :: projectDiffuseBand

    integer(i4b) :: i, nmaps
    logical(lgt) :: alm_in_
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: m, m_out

    nmaps     =  min(self%x%info%nmaps, map%info%nmaps)
    info      => comm_mapinfo(self%x%info%comm, map%info%nside, self%lmax_amp, &
         & self%x%info%nmaps, self%x%info%nmaps==3)
    m_out     => comm_map(info)
    m         => comm_map(map)
    alm_in_ = .false.; if (present(alm_in)) alm_in_ = alm_in

    ! Scale to correct frequency through multiplication with mixing matrix
    if (self%F_null(band)) then
       m_out%alm = 0.d0
    else
       ! Convolve with band-specific beam
       if (.not. alm_in) call m%Yt()
       call data(band)%B%conv(alm_in=.false., alm_out=(self%lmax_ind==0), trans=.true., map=m)

       if (self%lmax_ind == 0) then
          call m%alm_equal(m_out)
          do i = 1, nmaps
             m_out%alm(:,i) = m_out%alm(:,i) * self%F_mean(band,i)
          end do
       else
          call m%Y()
          m_out%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band)%p%map(:,1:nmaps)
          call m_out%YtW()
       end if
    end if

    if (.not. allocated(projectDiffuseBand)) allocate(projectDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
    projectDiffuseBand = m_out%alm

    deallocate(m, m_out, info)
    
  end function projectDiffuseBand

  subroutine precond_diff_comps(P, x)
    implicit none
    class(precondDiff),               intent(in)    :: P
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y

    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       call cr_extract_comp(diffComps(i)%p%id, x, alm)
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Multiply with preconditioner
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          !y(:,i,j) = matmul(P%invM(:,:,i,j), y(:,i,j))
          y(P%invM_(i,j)%ind,i,j) = matmul(P%invM_(i,j)%M, y(P%invM_(i,j)%ind,i,j))
       end do
    end do

    ! Reformat y(npre,nalm,nmaps) structure into linear array
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       allocate(alm(0:diffComps(i)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = y(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(i)%p%id, .false., alm, x)
       deallocate(alm)
    end do
    
    deallocate(y)

  end subroutine precond_diff_comps
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpDiffuseToFITS(self, postfix, dir)
    implicit none
    character(len=*),                        intent(in)           :: postfix
    class(comm_diffuse_comp),                intent(in)           :: self
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, m, ierr, unit
    real(dp)           :: vals(10)
    logical(lgt)       :: exist, first_call = .true.
    character(len=512) :: filename
    class(comm_map), pointer :: map

    if (trim(self%type) == 'md') then
       filename = trim(dir)//'/md_' // trim(postfix) // '.dat'
       vals = 0.d0
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l == 0) then                 ! Monopole
             vals(1)  = 1.d0/sqrt(4.d0*pi) * self%x%alm(i,1)             * self%RJ2unit_
             vals(5)  = 1.d0/sqrt(4.d0*pi) * self%mu%alm(i,1)            * self%RJ2unit_
             vals(9)  = 1.d0/sqrt(4.d0*pi) * sqrt(self%Cl%Dl(0,1))       * self%RJ2unit_
          end if
          if (l == 1 .and. m == -1) then   ! Y dipole
             vals(3)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_
             vals(7)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_
          end if
          if (l == 1 .and. m ==  0) then   ! Z dipole
             vals(4)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_
             vals(8)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_
          end if
          if (l == 1 .and. m ==  1) then   ! X dipole
             vals(2)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_
             vals(6)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_
             vals(10) =  1.d0/sqrt(4.d0*pi/3.d0) * sqrt(self%Cl%Dl(1,1))  * self%RJ2unit_
          end if
       end do
       call mpi_allreduce(MPI_IN_PLACE, vals, 10, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       if (self%x%info%myid == 0) then
          inquire(file=trim(filename), exist=exist)
          unit = getlun()
          if (first_call) then
             open(unit, file=trim(filename), recl=1024)
             write(unit,*) '# Band          M_def       X_def       Y_def       Z_def     M_mean      X_mean      Y_mean      Z_mean     M_rms       D_rms'
             first_call = .false.
          else
             open(unit, file=trim(filename), recl=1024, position='append')
          end if
          write(unit,'(a12,10e12.3)') trim(self%label), vals
          close(unit)
       end if
    else

       ! Write amplitude
       map => comm_map(self%x)
       map%alm = map%alm * self%RJ2unit_ * self%cg_scale  ! Output in requested units
       call self%B_out%conv(alm_in=.true., alm_out=.false., trans=.false., map=map)
       
       filename = trim(self%label) // '_' // trim(postfix) // '.fits'
       call map%Y
       call map%writeFITS(trim(dir)//'/'//trim(filename))
       deallocate(map)
       
       ! Write spectral index maps
       do i = 1, self%npar
          filename = trim(self%label) // '_' // trim(self%indlabel(i)) // '_' // &
               & trim(postfix) // '.fits'
          call self%theta(i)%p%Y
          call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
       end do
       
       ! Write mixing matrices
       if (self%output_mixmat) then
          do i = 1, numband
             if (self%F_null(i)) cycle
             filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
                  & trim(postfix) // '.fits'
             call self%F(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
          end do
       end if
    end if
        
  end subroutine dumpDiffuseToFITS

  subroutine add_to_npre(n, nside, lmax, nmaps)
    implicit none
    integer(i4b), intent(in) :: n, nside, lmax, nmaps
    npre      = npre + n
    nside_pre = min(lmax_pre, nside)
    lmax_pre  = max(lmax_pre, lmax)
    nmaps_pre = max(nmaps_pre, nmaps)
  end subroutine add_to_npre
  
end module comm_diffuse_comp_mod
