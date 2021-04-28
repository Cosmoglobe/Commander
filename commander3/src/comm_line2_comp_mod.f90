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
module comm_line2_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_line_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_line2_comp

  !**************************************************
  !           Line emission component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_line2_comp
     integer(i4b)                            :: ref_band       ! band number of the reference band
     real(dp)                                :: line2RJ_ref    ! reference band K km/s to uK_RJ ratio 
     integer(i4b), allocatable, dimension(:) :: ind2band       !array connecting index and band number
     real(dp),     allocatable, dimension(:) :: line2RJ        ! K km/s to uK_RJ ratio 
     real(dp),     allocatable, dimension(:) :: lineratio      !line ratio parameter
     real(dp),     allocatable, dimension(:) :: tilt           !bandpass tilt parameter
     real(dp),     allocatable, dimension(:) :: vmap_uni       !limits on velocity map
     class(comm_map),            pointer     :: vmap => null()  !velocity shift map
     class(comm_map),            pointer     :: numap => null() !frequency shift map
   contains
     procedure :: S    => evalSED
     procedure :: sampleSpecInd => sampleLineRatios
  end type comm_line2_comp

  interface comm_line2_comp
     procedure constructor
  end interface comm_line2_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_line2_comp), pointer   :: constructor

    integer(i4b) :: i, j, k, l, m, nline, b, n, ierr
    real(dp)     :: f
    real(dp)     :: light_speed
    logical(lgt) :: ref_exist
    logical(lgt) :: no_spatial !flag to check if there are to be spatial structure in theta
    character(len=512), allocatable, dimension(:) :: label    !band label
    real(dp),           allocatable, dimension(:) :: line2RJ  ! K km/s to uK_RJ ratio
    real(dp),           allocatable, dimension(:) :: lr_init, lr_mu, lr_sigma !line ratio values
    real(dp),           allocatable, dimension(:) :: tilt_init, tilt_mu, tilt_sigma !bandpass tilt values
    integer(i4b),       allocatable, dimension(:) :: poltype !polarization type
    type(comm_mapinfo), pointer :: info
    
    ! General parameters
    allocate(constructor)
    constructor%npar = 0 !temporary value to skip spectral parameter initialization in initDiffuse. This is done here instead
    call constructor%initDiffuse(cpar, id, id_abs)

    constructor%lmax_ind = -1 !assume spatial structure

    ! Read line template file
    call read_line2_template(trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), &
         & nline, label, lr_init, lr_mu, lr_sigma, tilt_init, tilt_mu, tilt_sigma, line2RJ, poltype)

    ! Check how many lines are included in current run
    n         = 0
    ref_exist = .false.
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) n = n+1
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) ref_exist = .true.
    end do
    if (.not. ref_exist) call report_error("Line component reference band does not exist")

    allocate(constructor%ind2band(n))
    constructor%npar = n

    allocate(constructor%theta_def(2*n), constructor%p_gauss(2,2*n), constructor%p_uni(2,2*n))
    allocate(constructor%poltype(n), constructor%indlabel(n), constructor%line2RJ(n))
    n         = 0
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) then
             n = n+1
             constructor%ind2band(n)  = i
             constructor%theta_def(n) = lr_init(j)
             constructor%p_gauss(1,n) = lr_mu(j)
             constructor%p_gauss(2,n) = lr_sigma(j)
             constructor%theta_def(n+constructor%npar) = tilt_init(j)
             constructor%p_gauss(1,n+constructor%npar) = tilt_mu(j)
             constructor%p_gauss(2,n+constructor%npar) = tilt_sigma(j)
             constructor%p_uni(1,n)   = 0.d0     !lower limit for line ratio
             constructor%p_uni(2,n)   = 100.d0   !upper limit for line ratio
             constructor%p_uni(1,n+constructor%npar)   = -100.d0 !lower limit for tilt
             constructor%p_uni(2,n+constructor%npar)   = 100.d0  !upper limit for tilt
             constructor%poltype(n)   = poltype(j)
             constructor%indlabel(n)  = label(j)
             constructor%line2RJ(n)   = line2RJ(j)
             exit
          end if
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) then
          constructor%ref_band    = i
          constructor%line2RJ_ref = constructor%line2RJ(n)
       end if
    end do

    !################################!
    ! Initialize spectral index maps
    !################################!
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    no_spatial=.false.
    !read velocity map
    allocate(constructor%vmap_uni(2))
    constructor%vmap_uni(1)=-3000.d0 !lower limit of vmap in km/s, 1% of c
    constructor%vmap_uni(2)=3000.d0  !upper limit of vmap in km/s, 1% of c
    if (trim(cpar%cs_vmap(id_abs)) == 'none') then
       no_spatial=.true.
       constructor%vmap => comm_map(info)
       constructor%vmap%map = 0.d0
    else
       constructor%vmap => comm_map(info, trim(cpar%datadir) // '/' // &
            & trim(cpar%cs_vmap(id_abs)))
       !truncate vmap on global limits
       constructor%vmap%map = min(constructor%vmap_uni(2),max(constructor%vmap_uni(1), &
            & constructor%vmap%map))
    end if
    
    !transform velocity map from km/s to total frequency shift
    constructor%numap => comm_map(info)
    light_speed = 299792.458d0 !speed of light in km/s

    !loop all maps and transform velocity into frequency shift [km/s -> GHz]. 
    do i = 1, info%nmaps
       constructor%numap%map(:,i) = ( sqrt((1.d0-constructor%vmap%map(:,i)/light_speed)/ &
            & (1+constructor%vmap%map(:,i)/light_speed))-1.d0 )*(constructor%nu_ref(i)/1.d9)
    end do

    allocate(constructor%lmax_ind_pol(3,constructor%npar))
    constructor%lmax_ind_pol = -1 !Default is pixel domain for line2 component
    if (allocated(constructor%lmax_ind_mix)) deallocate(constructor%lmax_ind_mix)
    allocate(constructor%lmax_ind_mix(3,constructor%npar))
    constructor%lmax_ind_mix = -1 !Default is pixel domain for line2 component
    allocate(constructor%pol_pixreg_type(3,constructor%npar))
    constructor%pol_pixreg_type = 0 !to not write unnecessary/non-allocated maps during FITS-dump
    !check if tilts are zero and not to be sampled. If so there is no spatial structure
    if (all( abs( constructor%theta_def(constructor%npar+1:) ) < 1.d-8 ) .and. &
         all( constructor%p_gauss(2,constructor%npar+1:) == 0.d0 )) then
       no_spatial = .true.
    end if

    if (no_spatial) then
       !make the lmax = 0 to speed up CG search
       constructor%lmax_ind_pol = 0 
       constructor%lmax_ind_mix = 0 
       constructor%lmax_ind = 0 

       if (info%myid == 0) write(*,*) 'No non-zero tilt parameter or the velocity map is zero. Setting lmax = 0'

       !update the info pointer with the new lmax
       info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
            & constructor%nmaps, constructor%pol)
    end if


    allocate(constructor%theta(n))
    do i = 1, n
       constructor%theta(i)%p     => comm_map(info)
       constructor%theta(i)%p%map = constructor%theta_def(i) + constructor%theta_def(i+constructor%npar) * constructor%numap%map 
       if (constructor%lmax_ind >= 0) call constructor%theta(i)%p%YtW_scalar
    end do

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    j = 1
    do l = 1, 3
       do i = 1, numband
          if (l > 1) then
             do m = 0,constructor%ndet
                constructor%F_int(l,i,m)%p => constructor%F_int(l-1,i,m)%p
             end do
             cycle
          end if
          if (any(constructor%ind2band == i)) then
             do k = 0, data(i)%ndet
                constructor%F_int(l,i,k)%p => comm_F_line(constructor, data(i)%bp(k)%p, .true., &
                     & constructor%line2RJ(j) / constructor%line2RJ_ref * data(i)%RJ2data(k), j)
             end do
             j = j+1
          else
             do k = 0, data(i)%ndet
                constructor%F_int(l,i,k)%p => comm_F_line(constructor, data(i)%bp(k)%p, .false., 0.d0, j)
             end do
          end if
       end do
    end do

    ! Initialize mixing matrix
    if (trim(cpar%init_chain_prefix) == 'none' &
         & .or. trim(constructor%init_from_HDF) == 'none') &
         & call constructor%updateMixmat

    deallocate(label, lr_init, lr_mu, lr_sigma, tilt_init, tilt_mu, tilt_sigma, line2RJ, poltype)

  end function constructor

  ! Definition:
  !    SED  = delta_{band,
  function evalSED(self, nu, band, pol, theta)
    class(comm_line2_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, ind

    if (band == self%ref_band) then
       evalSED = 1.d0
    else 
       do i = 1, self%npar
          if (band == self%ind2band(i)) exit
       end do
       if (i > self%npar) then
          evalSED = 0.d0
       else
          evalSED = theta(i) * self%line2RJ(i) / self%line2RJ_ref
       end if
    end if

  end function evalSED

  subroutine read_line2_template(filename, nline, label, lr_init, lr_mu, lr_sigma, tilt_init, tilt_mu, tilt_sigma, line2RJ, poltype)
    implicit none
    character(len=*),                              intent(in)  :: filename
    integer(i4b),                                  intent(out) :: nline
    character(len=512), allocatable, dimension(:), intent(out) :: label
    real(dp),           allocatable, dimension(:), intent(out) :: line2RJ
    real(dp),           allocatable, dimension(:), intent(out) :: lr_init, lr_mu, lr_sigma
    real(dp),           allocatable, dimension(:), intent(out) :: tilt_init, tilt_mu, tilt_sigma
    integer(i4b),       allocatable, dimension(:), intent(out) :: poltype

    integer(i4b)        :: i, unit
    character(len=1024) :: line

    unit  = getlun()

    ! Find number of lines
    nline = 0
    open(unit, file=trim(filename), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
    end do
1   close(unit)
    
    allocate(label(nline), line2RJ(nline), poltype(nline))
    allocate(lr_init(nline), lr_mu(nline), lr_sigma(nline))
    allocate(tilt_init(nline), tilt_mu(nline), tilt_sigma(nline))
    open(unit, file=trim(filename), recl=1024)
    nline = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
       read(line,*) label(nline), lr_init(nline), lr_mu(nline), lr_sigma(nline), &
            & tilt_init(nline), tilt_mu(nline), tilt_sigma(nline), & 
            & line2RJ(nline), poltype(nline)
    end do
2   close(unit)

  end subroutine read_line2_template

  ! Sample line ratios
  subroutine sampleLineRatios(self, cpar, handle, id, iter)
    implicit none
    class(comm_line2_comp),                   intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b)    :: i, j, l, n, m, band, ierr
    real(dp)        :: A, b, mu, sigma, par, sigma_p, scale, w
    class(comm_map), pointer :: invN_amp, amp, mask
    character(len=2) :: id_text
    
    band = self%ind2band(id)
    !if (band == self%ref_band) return

    ! Construct mask
!!$    if (associated(self%indmask)) then
!!$       if (data(band)%info%nside /= self%indmask%info%nside) then
!!$          call report_error("Mask udgrade in line_comp not yet supported.")
!!$       else
!!$          mask => self%indmask
!!$       end if
!!$    end if
    
    ! Compute likelihood term
    w            = self%theta(id)%p%map(1,1)
    amp          => comm_map(data(band)%info)
    invN_amp     => comm_map(data(band)%info)
    amp%map      =  self%getBand(band)/w
    invN_amp%map = amp%map
    call data(band)%N%invN(invN_amp)     ! Inverse noise variance weighted amplitude map
    
!!$    call int2string(id, id_text)
!!$    call mask%writeFITS('co_mask'//id_text//'.fits')
!!$    call amp%writeFITS('co_amp'//id_text//'.fits')
!!$    call data(band)%res%writeFITS('co_res'//id_text//'.fits')
!!$    call data(band)%N%invN_diag%writeFITS('co_invN'//id_text//'.fits')

    ! Reduce across processors
    if (associated(mask)) then
       A = sum(invN_amp%map * mask%map * amp%map)
       b = sum(invN_amp%map * mask%map * data(band)%res%map)
    else
       A = sum(invN_amp%map * amp%map)
       b = sum(invN_amp%map * data(band)%res%map)
    end if
    call mpi_allreduce(MPI_IN_PLACE, A, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    
    call amp%dealloc(); deallocate(amp)
    call invN_amp%dealloc(); deallocate(invN_amp)
    
    ! Compute new line ratio; just root processor
    if (self%x%info%myid == 0) then

!       write(*,*) 'A,b = ', A, b

       if (A > 0.d0) then
          mu    = b / A
          sigma = sqrt(1.d0 / A)
       else if (self%p_gauss(2,id) > 0.d0) then
          mu    = 0.d0
          sigma = 1.d30
       else
          mu    = self%p_uni(1,id) + (self%p_uni(2,id)-self%p_uni(1,id))*rand_uni(handle)
          sigma = 0.d0
       end if

!       write(*,*) '  mu, sigma = ', mu, sigma
       
       ! Add prior
       if (self%p_gauss(2,id) > 0.d0) then
          sigma_p = self%p_gauss(2,id) !/ sqrt(real(npix_reg,dp))
          mu      = (mu*sigma_p**2 + self%p_gauss(1,id) * sigma**2) / (sigma_p**2 + sigma**2)
          sigma   = sqrt(sigma**2 * sigma_p**2 / (sigma**2 + sigma_p**2))
       end if

!       write(*,*) '  mu_prior, sigma_prior = ', mu, sigma
       
       ! Draw sample
       par = -1.d30
       if (trim(self%operation) == 'optimize') then
          if (mu < self%p_uni(1,id)) then
             par = self%p_uni(1,id)
          else if (mu > self%p_uni(2,id)) then
             par = self%p_uni(2,id)
          else
             par = mu
          end if
       else
          do while (par < self%p_uni(1,id) .or. par > self%p_uni(2,id))
             if (mu < self%p_uni(1,id)) then
                par = rand_trunc_gauss(handle, mu, self%p_uni(1,id), sigma)
             else if (mu > self%p_uni(2,id)) then
                par = 2.d0*mu-rand_trunc_gauss(handle, mu, 2.d0*mu-self%p_uni(2,id), sigma)
             else
                par = mu + sigma * rand_gauss(handle)
             end if
          end do
       end if
       
       write(*,*) '  Line ratio i = ', id, ' = ', par
    end if
    
    ! Distribute new relative line ratio, and update
    call mpi_bcast(par, 1, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)

    if (band == self%ref_band) then
       self%x%map = self%x%map * par  ! Rescale amplitude map, but leave mixing matrix
       self%x%alm = self%x%alm * par  
       do i = 1, self%npar            ! Rescale line ratios at other frequencies
          if (self%ind2band(i) == self%ref_band) cycle
          self%theta(i)%p%map = self%theta(i)%p%map / par
          if (self%lmax_ind >= 0) then
             self%theta(i)%p%alm(:,1) = self%theta(i)%p%alm(:,1) / par
          end if
       end do
    else
       self%theta(id)%p%map = par
       if (self%lmax_ind >= 0) then
          self%theta(id)%p%alm(:,1) = par * sqrt(4.d0*pi)
       end if
    end if
    call self%updateMixmat()

  end subroutine sampleLineRatios
  
end module comm_line2_comp_mod
