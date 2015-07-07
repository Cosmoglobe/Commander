module comm_diffuse_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  use comm_data_mod
  use comm_F_int_mod
  use comm_Cl_mod
  use math_tools
  implicit none

  private
  public comm_diffuse_comp
  
  !**************************************************
  !            Diffuse component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nx, x0
     logical(lgt)       :: pol
     integer(i4b)       :: lmax_amp, lmax_ind, lpiv
     real(dp), allocatable, dimension(:,:) :: cls
     real(dp), allocatable, dimension(:,:) :: F_mean

     class(comm_map),               pointer     :: mask
     class(comm_map),               pointer     :: x      ! Spatial parameters
     class(comm_map),               pointer     :: inv_M  ! Preconditioner
     class(comm_B),                 pointer     :: B_out  ! Output beam
     class(comm_Cl),                pointer     :: Cl     ! Power spectrum
     class(map_ptr),  dimension(:), allocatable :: theta  ! Spectral parameters
     type(map_ptr),   dimension(:), allocatable :: F      ! Mixing matrix
     type(F_int_ptr), dimension(:), allocatable :: F_int  ! SED integrator
   contains
     procedure :: initDiffuse
     procedure :: updateMixmat
     procedure :: invM        => applyInvM
!!$     procedure :: dumpHDF  => dumpDiffuseToHDF
     procedure :: getBand     => evalDiffuseBand
     procedure :: projectBand => projectDiffuseBand
     procedure :: dumpFITS    => dumpDiffuseToFITS
     procedure :: initPrecond 
  end type comm_diffuse_comp

contains

  subroutine initDiffuse(self, cpar, id)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info
    
    call self%initComp(cpar, id)

    ! Initialize variables specific to diffuse source type
    self%pol      = cpar%cs_polarization(id)
    self%nside    = cpar%cs_nside(id)
    self%lmax_amp = cpar%cs_lmax_amp(id)
    self%lmax_ind = cpar%cs_lmax_ind(id)
    self%cltype   = cpar%cs_cltype(id)
    self%nmaps    = 1; if (self%pol) self%nmaps = 3
    info          => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id)) == 'zero' .or. trim(cpar%cs_input_amp(id)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       self%x => comm_map(info, cpar%cs_input_amp(id))
       call self%x%YtW
    end if
    self%ncr = size(self%x%alm)

    ! Initialize preconditioner structure
    self%inv_M => comm_map(info)

    ! Allocate mixing matrix
    allocate(self%F(numband), self%F_mean(numband,self%nmaps))
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & self%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
       self%F(i)%p => comm_map(info)
    end do

    ! Initialize output beam
    self%B_out => comm_B_bl(cpar, self%x%info, 0, fwhm=cpar%cs_fwhm(id))

    ! Initialize power spectrum
    self%Cl => comm_Cl(cpar, self%x%info, id)
    
  end subroutine initDiffuse

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
    allocate(theta_p(self%nmaps,self%npar))
    do i = 1, numband

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
          self%F(i)%p%map(j,1) = self%F_int(i)%p%eval(theta_p(:,1)) * data(i)%RJ2data()

          ! Polarization
          if (self%nmaps == 3) then
             ! Stokes Q
             if (all(self%poltype < 2)) then
                self%F(i)%p%map(j,2) = self%F(i)%p%map(j,1)
             else
                self%F(i)%p%map(j,2) = self%F_int(i)%p%eval(theta_p(:,2)) * self%RJ2unit()
             end if
       
             ! Stokes U
             if (all(self%poltype < 3)) then
                self%F(i)%p%map(j,3) = self%F(i)%p%map(j,2)
             else
                self%F(i)%p%map(j,3) = self%F_int(i)%p%eval(theta_p(:,3)) * self%RJ2unit()
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

  end subroutine updateMixmat

  function evalDiffuseBand(self, band, amp_in, pix)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    real(dp),        dimension(:,:), allocatable                        :: evalDiffuseBand

    integer(i4b) :: i, j, np, nmaps, lmax, nmaps_comp
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: m

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
    call m%Y

    ! Scale to correct frequency through multiplication with mixing matrix
    m%map = m%map * self%F(band)%p%map

    ! Convolve with band-specific beam
    call data(band)%B%conv(alm_in=.false., alm_out=.false., trans=.false., map=m)

    ! Return pixelized map
    allocate(evalDiffuseBand(0:data(band)%info%np-1,data(band)%info%nmaps))
    evalDiffuseBand = m%map

    ! Clean up
    deallocate(m, info)

  end function evalDiffuseBand

  ! Return component projected from map
  function projectDiffuseBand(self, band, map)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    real(dp),        dimension(:,:), allocatable                        :: projectDiffuseBand

    class(comm_map), pointer :: m

    m => comm_map(map)

    ! Convolve with band-specific beam
    call data(band)%B%conv(alm_in=.false., alm_out=.false., trans=.true., map=m)

    ! Scale to correct frequency through multiplication with mixing matrix
    m%map = m%map * self%F(band)%p%map

    ! Extract spherical harmonics coefficients
    allocate(projectDiffuseBand(self%x%info%nalm,self%x%info%nmaps))
    call m%Yt
    projectDiffuseBand = m%alm

    deallocate(m)
    
  end function projectDiffuseBand

  ! Return component projected from map
  subroutine initPrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout) :: self

    integer(i4b) :: i, j, k, l, m, ind
    real(dp)     :: alm(self%inv_M%info%nmaps)

    self%inv_M%map = 0.d0
    do j = 0, self%x%info%nalm-1
       call self%x%info%i2lm(j,l,m)
       ! Compute inv_M = sum_nu F^t * B^t * invN * B * F
       do i = 1, numband
          call data(i)%N%invN_diag%info%lm2i(l,m,k)
          if (k > -1) then
             alm = data(i)%N%invN_diag%alm(k,:) * data(i)%B%b_l(l,:)**2 * self%F_mean(i,:)**2
             if (trim(self%cltype) /= 'none') then
                ! inv_M = C^(1/2) * F^t * B^t * invN * B * F * C^(1/2)
             end if
             self%inv_M%alm(j,:) = self%inv_M%alm(j,:) + alm
          end if
       end do
    end do

  end subroutine initPrecond

  subroutine applyInvM(self, Nscale, alm)
    implicit none

    class(comm_diffuse_comp),                   intent(in)    :: self
    real(dp),                                   intent(in)    :: Nscale
    real(dp),                 dimension(0:,1:), intent(inout) :: alm

    integer(i4b) :: i, l, j, ierr
    real(dp)     :: M(self%nmaps,self%nmaps), a(self%nmaps,1)
    
    if (trim(self%cltype) /= 'none') then
       do i = 0, self%x%info%nalm-1
          l      = self%x%info%lm(1,i)
          a(:,1) = matmul(self%Cl%sqrtS_mat(:,:,l), self%inv_M%alm(i,:))
          M      = Nscale * matmul(a,transpose(a))
          do j = 1, self%nmaps
             M(j,j) = M(j,j) + 1.d0
          end do
          call cholesky_decompose_single(M)
          call cholesky_solve(M, alm(i,:), a(:,1))
          alm(i,:) = a(:,1)
       end do
    else
       alm = alm / (Nscale * self%inv_M%alm)
    end if

  end subroutine applyInvM
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpDiffuseToFITS(self, postfix, dir)
    implicit none
    character(len=*),                        intent(in)           :: postfix
    class(comm_diffuse_comp),                intent(in)           :: self
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i
    character(len=512) :: filename
    class(comm_map), pointer :: map

    ! Write amplitude
    map => comm_map(self%x)
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
    do i = 1, numband
       filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
            & trim(postfix) // '.fits'
       call self%F(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
    end do
        
  end subroutine dumpDiffuseToFITS
  
end module comm_diffuse_comp_mod
