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
  function constructor(cpar, id, id_abs, band, label, mu, rms, def) result(comm_md_ptr)
    implicit none
    type(comm_params),              intent(in) :: cpar
    integer(i4b),                   intent(in) :: id, id_abs, band
    character(len=*),               intent(in) :: label
    real(dp),         dimension(4), intent(in) :: mu, def
    real(dp),         dimension(2), intent(in) :: rms
    class(comm_md_comp), pointer               :: comm_md_ptr

    integer(i4b) :: i, j, k, l, m, n
    character(len=16), dimension(1000) :: comp_label
    type(comm_mapinfo), pointer :: info => null()

!    write(*,*) 'mu', trim(label), real(mu,sp)
!    write(*,*) 'rms', trim(label), real(rms,sp)

    ! General parameters
    allocate(comm_md_ptr)

    ! Initialize comm_comp_mod parameters
    comm_md_ptr%id              = id
    comm_md_ptr%active          = cpar%cs_include(id_abs)
    comm_md_ptr%label           = data(band)%label
    comm_md_ptr%type            = cpar%cs_type(id_abs)
    comm_md_ptr%class           = cpar%cs_class(id_abs)
    comm_md_ptr%unit            = data(band)%unit
    comm_md_ptr%nu_ref          = data(band)%bp(0)%p%nu_c
    comm_md_ptr%cg_scale        = 1.d0
    comm_md_ptr%myid            = cpar%myid_chain
    comm_md_ptr%comm            = cpar%comm_chain
    comm_md_ptr%numprocs        = cpar%numprocs_chain
    comm_md_ptr%init_from_HDF   = cpar%cs_initHDF(id_abs)
    comm_md_ptr%lmax_pre_lowl   = -1
    comm_md_ptr%lmax_def        = -1
    comm_md_ptr%nside_def       = 0
    comm_md_ptr%fwhm_def        = 0.d0
    comm_md_ptr%mono_prior_type = 'none'

    call get_tokens(cpar%output_comps, ",", comp_label, n)
    comm_md_ptr%output = .false.
    do i = 1, n
       if (trim(comp_label(i)) == trim(comm_md_ptr%label) .or. trim(comp_label(i)) == 'all') then
          comm_md_ptr%output = .true.
          exit
       end if
    end do

    !comm_md_ptr%ref_band = band

    ! Set up conversion factor between RJ and native component unit
    select case (trim(comm_md_ptr%unit))
    case ('uK_cmb')
       comm_md_ptr%RJ2unit_ = data(band)%bp(0)%p%a2t
    case ('mK_cmb')
       comm_md_ptr%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-3
    case ('K_cmb')
       comm_md_ptr%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-6
    case ('MJy/sr') 
       comm_md_ptr%RJ2unit_ = data(band)%bp(0)%p%a2t / data(band)%bp(0)%p%f2t
    case ('uK_RJ') 
       comm_md_ptr%RJ2unit_ = 1.d0
    case ('K km/s') 
       comm_md_ptr%RJ2unit_ = 1.d0
    case default
       call report_error('Unsupported unit: ' // trim(comm_md_ptr%unit))
    end select

    ! Initialize comm_diffuse_comp_mod parameters
    comm_md_ptr%pol      = .false.
    comm_md_ptr%nside    = data(band)%map%info%nside
    comm_md_ptr%lmax_amp = 1
    comm_md_ptr%l_apod   = 2
    comm_md_ptr%lmax_ind = 0
    comm_md_ptr%cltype   = 'binned'
    comm_md_ptr%nmaps    = 1
    allocate(comm_md_ptr%lmax_ind_mix(3,1))
    comm_md_ptr%lmax_ind_mix = 0

    !info          => comm_mapinfo(cpar%comm_chain, comm_md_ptr%nside, comm_md_ptr%lmax_amp, &
    !     & comm_md_ptr%nmaps, comm_md_ptr%pol)
    info          => comm_mapinfo(cpar%comm_chain, 128, comm_md_ptr%lmax_amp, &
         & comm_md_ptr%nmaps, comm_md_ptr%pol)

    ! Diffuse preconditioner variables
    call add_to_npre(1,comm_md_ptr%nside,1,1)

    ! Initialize amplitude and prior maps
    comm_md_ptr%x   => comm_map(info)
    comm_md_ptr%mu  => comm_map(info)
    comm_md_ptr%ncr =  size(comm_md_ptr%x%alm)
    do i = 0, comm_md_ptr%x%info%nalm-1
       call comm_md_ptr%x%info%i2lm(i,l,m)
       if (l == 0) then ! Monopole
          comm_md_ptr%x%alm(i,1)  = sqrt(4.d0*pi) * def(1) / comm_md_ptr%RJ2unit_(1)
          comm_md_ptr%mu%alm(i,1) = sqrt(4.d0*pi) * mu(1)  / comm_md_ptr%RJ2unit_(1)
       end if
       if (l == 1 .and. m == -1) then ! Y dipole
          comm_md_ptr%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(3) / comm_md_ptr%RJ2unit_(1)
          comm_md_ptr%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(3)  / comm_md_ptr%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  0) then ! Z dipole
          comm_md_ptr%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(4) / comm_md_ptr%RJ2unit_(1)
          comm_md_ptr%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(4)  / comm_md_ptr%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  1) then ! X dipole
          comm_md_ptr%x%alm(i,1)  = -sqrt(4.d0*pi/3.d0) * def(2) / comm_md_ptr%RJ2unit_(1)
          comm_md_ptr%mu%alm(i,1) = -sqrt(4.d0*pi/3.d0) * mu(2)  / comm_md_ptr%RJ2unit_(1)
       end if
    end do


    ! Allocate mixing matrix
    comm_md_ptr%ndet = maxval(data%ndet)
    allocate(comm_md_ptr%F(numband,0:comm_md_ptr%ndet), comm_md_ptr%F_mean(numband,0:comm_md_ptr%ndet,comm_md_ptr%nmaps))
    allocate(comm_md_ptr%F_null(numband,0:comm_md_ptr%ndet))
    do i = 1, numband
       if (i == band) then
          info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
               & comm_md_ptr%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
          comm_md_ptr%F(i,0)%p      => comm_map(info)
          comm_md_ptr%F(i,0)%p%map  = comm_md_ptr%RJ2unit_(1)
          comm_md_ptr%F(i,0)%p%alm  = comm_md_ptr%RJ2unit_(1) * sqrt(4.d0*pi)
          comm_md_ptr%F_null(i,:)   = .false.
          comm_md_ptr%F_mean(i,:,:) = comm_md_ptr%RJ2unit_(1)
          do j = 1, data(i)%ndet
             comm_md_ptr%F(i,j)%p      => comm_md_ptr%F(i,0)%p
          end do
       else
          do j = 0, data(i)%ndet
             comm_md_ptr%F_null(i,j)   = .true.
             !comm_md_ptr%F(i)%p%map  = 0.d0
             !comm_md_ptr%F(i)%p%alm  = 0.d0
             comm_md_ptr%F_mean(i,j,:) = 0.d0
          end do
       end if
    end do

    ! Initialize output beam
    comm_md_ptr%B_out => comm_B_bl(cpar, comm_md_ptr%x%info, 0, 0, fwhm=0.d0, init_realspace=.false.)

    ! Initialize power spectrum
    allocate(comm_md_ptr%Cl)
    comm_md_ptr%Cl%type   = 'binned'
    comm_md_ptr%Cl%info   => comm_md_ptr%x%info
    comm_md_ptr%Cl%label  = 'md_'//trim(data(band)%label)
    comm_md_ptr%Cl%lmax   = 1
    comm_md_ptr%Cl%nmaps  = 1
    comm_md_ptr%Cl%nspec  = 1
    comm_md_ptr%Cl%outdir = 'none'
    allocate(comm_md_ptr%Cl%Dl(0:1,1), comm_md_ptr%Cl%sqrtS_mat(1,1,0:1))
    allocate(comm_md_ptr%Cl%sqrtInvS_mat(1,1,0:1), comm_md_ptr%Cl%S_mat(1,1,0:1))
    comm_md_ptr%Cl%Dl(0,1) = 4.d0*pi      * rms(1)**2 / comm_md_ptr%RJ2unit_(1)**2
    comm_md_ptr%Cl%Dl(1,1) = 4.d0*pi/3.d0 * rms(2)**2 / comm_md_ptr%RJ2unit_(1)**2
    comm_md_ptr%Cl%S_mat(1,1,0:1)        = rms**2     / comm_md_ptr%RJ2unit_(1)**2
    comm_md_ptr%Cl%sqrtS_mat(1,1,0:1)    = rms        / comm_md_ptr%RJ2unit_(1)
    if (rms(1) > 0.d0) comm_md_ptr%Cl%sqrtInvS_mat(1,1,0) = 1.d0/rms(1) * comm_md_ptr%RJ2unit_(1)
    if (rms(2) > 0.d0) comm_md_ptr%Cl%sqrtInvS_mat(1,1,1) = 1.d0/rms(2) * comm_md_ptr%RJ2unit_(1)
    
    ! Initialize md_mod specific parameters
    comm_md_ptr%npar = 0

    ! Precompute mixmat integrator for each band
    allocate(comm_md_ptr%F_int(3,numband,0:comm_md_ptr%ndet))
    do k = 1, 3
       do i = 1, numband
          if (i == band) then
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (comm_md_ptr%nu_ref(k) == comm_md_ptr%nu_ref(k-1)) then
                      comm_md_ptr%F_int(k,i,j)%p => comm_md_ptr%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
                comm_md_ptr%F_int(k,i,j)%p => comm_F_line(comm_md_ptr, data(i)%bp(j)%p, .true., 1.d0, -1)
             end do
          else
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (comm_md_ptr%nu_ref(k) == comm_md_ptr%nu_ref(k-1)) then
                      comm_md_ptr%F_int(k,i,j)%p => comm_md_ptr%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
                comm_md_ptr%F_int(k,i,j)%p => comm_F_line(comm_md_ptr, data(i)%bp(j)%p, .true., 0.d0, -1)
             end do
          end if
       end do
    end do

    ! Set up CG sampling groups                                                                                                                                             
    allocate(comm_md_ptr%active_samp_group(cpar%cg_num_samp_groups))
    comm_md_ptr%active_samp_group = .false.
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(comm_md_ptr%label) == trim(comp_label(j))) then
             comm_md_ptr%active_samp_group(i) = .true.
             if (n == 1) comm_md_ptr%cg_unique_sampgroup = i ! Dedicated sampling group for this component  
             exit
          end if
       end do
    end do


    ! Set up default values for prior sampling (to be potentially changed at end of init)  
    comm_md_ptr%mono_from_prior=.false.
    comm_md_ptr%mono_alm = 0.d0
    
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

  function initialize_md_comps(cpar, id, id_abs, n) result(md_comps)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_md_comp), pointer   :: md_comps

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
          md_comps => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
       else
          c => comm_md_comp(cpar, id+n, id_abs, i, label, mu, rms, def)
          call md_comps%add(c)
       end if
       n = n+1
    end do
1   close(unit)
  
    if (n < numband .and. cpar%myid == 0) then
       write(*,'(a,i6)') ' | Warning: Number of channels without a monopole/dipole definition = ', numband-n
    end if
  
  end function initialize_md_comps
  
end module comm_md_comp_mod
