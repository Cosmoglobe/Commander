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
module comm_md_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_F_line_mod
  use comm_data_mod
  use comm_bp_utils
  implicit none

  private
  public comm_md_comp, initialize_md_comps

  !**************************************************
  !           Monopole/dipole component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_md_comp
     !integer(i4b)                            :: ref_band
     logical(lgt) :: mono_from_prior  ! true if the band used as zero-level prior
     real(dp)     :: mono_alm         ! alm value of monopole for when we CG-sample,
                                      ! can revert to pre CG sampling value afterwards
   contains
     procedure :: S    => evalSED
  end type comm_md_comp

  interface comm_md_comp
     procedure constructor
  end interface comm_md_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs, band, label, mu, rms, def)
    implicit none
    type(comm_params),              intent(in) :: cpar
    integer(i4b),                   intent(in) :: id, id_abs, band
    character(len=*),               intent(in) :: label
    real(dp),         dimension(4), intent(in) :: mu, def
    real(dp),         dimension(2), intent(in) :: rms
    class(comm_md_comp), pointer               :: constructor

    integer(i4b) :: i, j, k, l, m, n
    character(len=16), dimension(1000) :: comp_label
    type(comm_mapinfo), pointer :: info => null()

!    write(*,*) 'mu', trim(label), real(mu,sp)
!    write(*,*) 'rms', trim(label), real(rms,sp)
    ! General parameters
    write(*,*) 's1'
    allocate(constructor)

    ! Initialize comm_comp_mod parameters
    constructor%id              = id
    constructor%active          = cpar%cs_include(id_abs)
    constructor%label           = data(band)%label
    constructor%type            = cpar%cs_type(id_abs)
    constructor%class           = cpar%cs_class(id_abs)
    constructor%unit            = data(band)%unit
    constructor%nu_ref          = data(band)%bp(0)%p%nu_c
    constructor%cg_scale        = 1.d0
    constructor%myid            = cpar%myid_chain
    constructor%comm            = cpar%comm_chain
    constructor%numprocs        = cpar%numprocs_chain
    constructor%init_from_HDF   = cpar%cs_initHDF(id_abs)
    constructor%lmax_pre_lowl   = -1
    constructor%lmax_def        = -1
    constructor%nside_def       = 0
    constructor%fwhm_def        = 0.d0
    constructor%mono_prior_type = 'none'
    precond_type                = cpar%cg_precond

    call get_tokens(cpar%output_comps, ",", comp_label, n)
    constructor%output = .false.
    do i = 1, n
       if (trim(comp_label(i)) == trim(constructor%label) .or. trim(comp_label(i)) == 'all') then
          constructor%output = .true.
          exit
       end if
    end do

    write(*,*) 's2'

    !constructor%ref_band = band

    ! Set up conversion factor between RJ and native component unit
    select case (trim(constructor%unit))
    case ('uK_cmb')
       constructor%RJ2unit_ = data(band)%bp(0)%p%a2t
    case ('mK_cmb')
       constructor%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-3
    case ('K_cmb')
       constructor%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-6
    case ('MJy/sr') 
      !  constructor%RJ2unit_ = data(band)%bp(0)%p%a2t / data(band)%bp(0)%p%f2t
       constructor%RJ2unit_ = data(band)%bp(0)%p%a2f
    case ('uK_RJ') 
       constructor%RJ2unit_ = 1.d0
    case ('Kkm/s') 
       constructor%RJ2unit_ = 1.d0
       ! constructor%RJ2unit_ = 1.d0 / (constructor%nu_ref/c * 1d9)
!       write(*,*) 'Kkm/s not yet supported in md_mod'
!       call mpi_finalize(i)
!       stop
    case default
       call report_error('Unsupported unit: ' // trim(constructor%unit))
    end select

    ! Initialize comm_diffuse_comp_mod parameters
    constructor%pol      = .false.
    constructor%nside    = data(band)%map%info%nside
    constructor%lmax_amp = 1
    constructor%l_apod   = 2
    constructor%lmax_ind = 0
    constructor%cltype   = 'binned'
    constructor%nmaps    = 1
    allocate(constructor%lmax_ind_mix(3,1))
    constructor%lmax_ind_mix = 0

    !info          => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_amp, &
    !     & constructor%nmaps, constructor%pol)
    info          => comm_mapinfo(cpar%comm_chain, 128, constructor%lmax_amp, &
         & constructor%nmaps, constructor%pol)

    write(*,*) 's3'

    ! Diffuse preconditioner variables
    call add_to_npre(1,constructor%nside,1,1)

    ! Initialize amplitude and prior maps
    constructor%x   => comm_map(info)
    constructor%mu  => comm_map(info)
    constructor%ncr =  size(constructor%x%alm)
    do i = 0, constructor%x%info%nalm-1
       call constructor%x%info%i2lm(i,l,m)
       if (l == 0) then ! Monopole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi) * def(1) / constructor%RJ2unit_(1)
          constructor%mu%alm(i,1) = sqrt(4.d0*pi) * mu(1)  / constructor%RJ2unit_(1)
       end if
       if (l == 1 .and. m == -1) then ! Y dipole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(3) / constructor%RJ2unit_(1)
          constructor%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(3)  / constructor%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  0) then ! Z dipole
          constructor%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(4) / constructor%RJ2unit_(1)
          constructor%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(4)  / constructor%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  1) then ! X dipole
          constructor%x%alm(i,1)  = -sqrt(4.d0*pi/3.d0) * def(2) / constructor%RJ2unit_(1)
          constructor%mu%alm(i,1) = -sqrt(4.d0*pi/3.d0) * mu(2)  / constructor%RJ2unit_(1)
       end if
    end do


    ! Allocate mixing matrix
    constructor%ndet = maxval(data%ndet)
    allocate(constructor%F(numband,0:constructor%ndet), constructor%F_mean(numband,0:constructor%ndet,constructor%nmaps))
    allocate(constructor%F_null(numband,0:constructor%ndet))
    do i = 1, numband
       if (i == band) then
          info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
               & constructor%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
          constructor%F(i,0)%p      => comm_map(info)
          constructor%F(i,0)%p%map  = constructor%RJ2unit_(1)
          constructor%F(i,0)%p%alm  = constructor%RJ2unit_(1) * sqrt(4.d0*pi)
          constructor%F_null(i,:)   = .false.
          constructor%F_mean(i,:,:) = constructor%RJ2unit_(1)
          do j = 1, data(i)%ndet
             constructor%F(i,j)%p      => constructor%F(i,0)%p
          end do
       else
          do j = 0, data(i)%ndet
             constructor%F_null(i,j)   = .true.
             !constructor%F(i)%p%map  = 0.d0
             !constructor%F(i)%p%alm  = 0.d0
             constructor%F_mean(i,j,:) = 0.d0
          end do
       end if
    end do

    write(*,*) 's4'

    ! Initialize output beam
    constructor%B_out => comm_B_bl(cpar, constructor%x%info, 0, 0, fwhm=0.d0, init_realspace=.false.)

    ! Initialize power spectrum
    allocate(constructor%Cl)
    constructor%Cl%type   = 'binned'
    constructor%Cl%info   => constructor%x%info
    constructor%Cl%label  = 'md_'//trim(data(band)%label)
    constructor%Cl%lmax   = 1
    constructor%Cl%nmaps  = 1
    constructor%Cl%nspec  = 1
    constructor%Cl%outdir = 'none'
    allocate(constructor%Cl%Dl(0:1,1), constructor%Cl%sqrtS_mat(1,1,0:1))
    allocate(constructor%Cl%sqrtInvS_mat(1,1,0:1), constructor%Cl%S_mat(1,1,0:1))
    constructor%Cl%Dl(0,1) = 4.d0*pi      * rms(1)**2 / constructor%RJ2unit_(1)**2
    constructor%Cl%Dl(1,1) = 4.d0*pi/3.d0 * rms(2)**2 / constructor%RJ2unit_(1)**2
    constructor%Cl%S_mat(1,1,0:1)        = rms**2     / constructor%RJ2unit_(1)**2
    constructor%Cl%sqrtS_mat(1,1,0:1)    = rms        / constructor%RJ2unit_(1)
    if (rms(1) > 0.d0) constructor%Cl%sqrtInvS_mat(1,1,0) = 1.d0/rms(1) * constructor%RJ2unit_(1)
    if (rms(2) > 0.d0) constructor%Cl%sqrtInvS_mat(1,1,1) = 1.d0/rms(2) * constructor%RJ2unit_(1)
    
    ! Initialize md_mod specific parameters
    constructor%npar = 0

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    do k = 1, 3
       do i = 1, numband
          if (i == band) then
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                      constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
                constructor%F_int(k,i,j)%p => comm_F_line(constructor, data(i)%bp(j)%p, .true., 1.d0, -1)
             end do
          else
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                      constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
                constructor%F_int(k,i,j)%p => comm_F_line(constructor, data(i)%bp(j)%p, .true., 0.d0, -1)
             end do
          end if
       end do
    end do

    ! Set up CG sampling groups                                                                                                                                             
    allocate(constructor%active_samp_group(cpar%cg_num_samp_groups))
    constructor%active_samp_group = .false.
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(constructor%label) == trim(comp_label(j))) then
             constructor%active_samp_group(i) = .true.
             if (n == 1) constructor%cg_unique_sampgroup = i ! Dedicated sampling group for this component  
             exit
          end if
       end do
    end do


    ! Set up default values for prior sampling (to be potentially changed at end of init)  
    constructor%mono_from_prior=.false.
    constructor%mono_alm = 0.d0

    write(*,*) 's5', trim(constructor%label)
    
  end function constructor

  ! Definition:
  !    SED  = delta_{band,ref_band}
  function evalSED(self, nu, band, pol, theta)
    class(comm_md_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, ind

    evalSED = 1.d0

  end function evalSED

  function initialize_md_comps(cpar, id, id_abs, n)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_md_comp), pointer   :: initialize_md_comps

    integer(i4b)        :: i, unit
    real(dp)            :: mu(4), rms(2), def(4)
    character(len=1024) :: line, label
    class(comm_comp), pointer :: c => null()

    unit  = getlun()
    ! Find number of lines
    n = 0
    open(unit, file=trim(cpar%cs_SED_template(1,id_abs)), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       read(line,*) label, def, mu, rms
       do i = 1, numband
          if (trim(label) == trim(data(i)%label)) exit
       end do
       if (i > numband) cycle

       if (n == 0) then
          initialize_md_comps => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
       else
          c => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
          call initialize_md_comps%add(c)
       end if
!       write(*,*) 'cc', trim(c%label)
       n = n+1
    end do
1   close(unit)
  
    if (n < numband .and. cpar%myid == 0) then
       write(*,'(a,i6)') ' | Warning: Number of channels without a monopole/dipole definition = ', numband-n
       stop
    end if

  end function initialize_md_comps
  
end module comm_md_comp_mod
