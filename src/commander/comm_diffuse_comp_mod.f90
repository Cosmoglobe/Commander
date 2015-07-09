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
  public comm_diffuse_comp, precondDiff, precond_diff_comps
  
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
!!$     procedure :: dumpHDF  => dumpDiffuseToHDF
     procedure :: getBand     => evalDiffuseBand
     procedure :: projectBand => projectDiffuseBand
     procedure :: dumpFITS    => dumpDiffuseToFITS
     procedure :: initPrecond 
  end type comm_diffuse_comp

  type diff_ptr
     class(comm_diffuse_comp), pointer :: p
  end type diff_ptr
  
  type precondDiff
     integer(i4b) :: lmax, nside, nmaps
     real(dp)     :: Nscale
     real(dp), allocatable, dimension(:,:,:,:) :: invM0 ! (npre,npre,0:nalm-1,nmaps)
     real(dp), allocatable, dimension(:,:,:,:) :: invM  ! (npre,npre,0:nalm-1,nmaps)
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
  class(comm_mapinfo), pointer                   :: info_pre
  class(diff_ptr),     allocatable, dimension(:) :: diffComps
  
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

    ! Diffuse preconditioner variables
    npre      = npre+1
    lmax_pre  = max(lmax_pre,  self%lmax_amp)
    nside_pre = min(nside_pre, self%nside)
    nmaps_pre = max(nmaps_pre, self%nmaps)

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

  function constPreDiff(comm, Nscale)
    implicit none
    integer(i4b),                intent(in) :: comm
    real(dp),                    intent(in) :: Nscale
    class(precondDiff), pointer             :: constPreDiff

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m
    class(comm_comp),         pointer :: c
    class(comm_diffuse_comp), pointer :: p1, p2

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
    allocate(constPreDiff%invM0(npre,npre,0:info_pre%nalm-1,info_pre%nmaps))
    allocate(constPreDiff%invM(npre,npre,0:info_pre%nalm-1,info_pre%nmaps))
    constPreDiff%invM0 = 0.d0
    do j = 1, info_pre%nmaps
       do i1 = 0, info_pre%nalm-1
          call info_pre%i2lm(i1, l, m)
          do q = 1, numband
             call data(q)%info%lm2i(l,m,i2)
             if (i2 == -1) cycle
             do k1 = 1, npre
                p1 => diffComps(k1)%p
                do k2 = 1, npre
                   p2 => diffComps(k2)%p
                   constPreDiff%invM0(k1,k2,i1,j) = constPreDiff%invM0(k1,k2,i1,j) + &
                        & data(q)%N%invN_diag%alm(i2,j) * &                         ! invN_{lm,lm}
                        & data(q)%B%b_l(l,j)**2 * &                                 ! b_l^2
                        & p1%F_mean(q,j) * p2%F_mean(q,j) ! F(c1)*F(c2)
                end do
             end do
          end do
       end do
    end do

  end function constPreDiff

  subroutine updateDiffPrecond(self)
    implicit none
    class(precondDiff), intent(inout) :: self

    integer(i4b) :: i, j, k1, k2
    
    self%invM = self%invM0
    
    ! Right-multiply with sqrt(Cl)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       do k2 = 1, npre
          call diffComps(k1)%p%Cl%sqrtS(alm=self%invM(k2,k1,:,:))
       end do
    end do

    ! Left-multiply with sqrt(Cl)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       do k2 = 1, npre
          call diffComps(k1)%p%Cl%sqrtS(alm=self%invM(k1,k2,:,:))
       end do
    end do

    ! Add unity 
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       self%invM(k1,k1,:,:) = self%invM(k1,k1,:,:) + 1.d0
    end do

    ! Invert
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          call invert_matrix_with_mask(self%invM(:,:,i,j))
       end do
    end do    

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

    ! Scale to correct frequency through multiplication with mixing matrix
    if (self%lmax_ind == 0) then
       do i = 1, m%info%nmaps
          m%alm(:,i) = m%alm(:,i) * self%F_mean(band,i)
       end do
    else
       call m%Y
       m%map = m%map * self%F(band)%p%map
    end if

    ! Convolve with band-specific beam
    call data(band)%B%conv(alm_in=(self%lmax_ind==0), alm_out=.false., trans=.false., map=m)

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

    integer(i4b) :: i, nmaps
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: m, m_out

    nmaps     =  min(self%x%info%nmaps, map%info%nmaps)
    info      => comm_mapinfo(self%x%info%comm, map%info%nside, self%lmax_amp, &
         & self%x%info%nmaps, self%x%info%nmaps==3)
    m_out     => comm_map(info)
    m         => comm_map(map)

    ! Convolve with band-specific beam
    call data(band)%B%conv(alm_in=.false., alm_out=(self%lmax_ind==0), trans=.true., map=m)

    ! Scale to correct frequency through multiplication with mixing matrix
    if (self%lmax_ind == 0) then
       call m%alm_equal(m_out)
       do i = 1, nmaps
          m_out%alm(:,i) = m_out%alm(:,i) * self%F_mean(band,i)
       end do
       m_out%alm = m_out%alm * self%x%info%npix/(4.d0*pi) ! From beam convolution
    else
       m_out%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band)%p%map(:,1:nmaps)
       call m_out%Yt
    end if

    allocate(projectDiffuseBand(0:self%x%info%nalm-1,self%x%info%nmaps))
    projectDiffuseBand = m_out%alm

    deallocate(m, m_out, info)
    
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
          y(:,i,j) = matmul(P%invM(:,:,i,j), y(:,i,j))
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
