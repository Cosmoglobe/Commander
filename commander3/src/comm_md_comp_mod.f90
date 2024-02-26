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
  use comm_comp_interface_mod
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
     procedure :: S    => evalSED_md
  end type comm_md_comp

  interface comm_md_comp
     procedure constructor_md
  end interface comm_md_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_md(cpar, id, id_abs, band, label, mu, rms, def) result(c)
    implicit none
    type(comm_params),              intent(in) :: cpar
    integer(i4b),                   intent(in) :: id, id_abs, band
    character(len=*),               intent(in) :: label
    real(dp),         dimension(4), intent(in) :: mu, def
    real(dp),         dimension(2), intent(in) :: rms
    class(comm_md_comp), pointer               :: c

    integer(i4b) :: i, j, k, l, m, n
    character(len=16), dimension(1000) :: comp_label
    type(comm_mapinfo), pointer :: info => null()

    ! General parameters
    allocate(c)

    ! Initialize comm_comp_mod parameters
    c%id              = id
    c%active          = cpar%cs_include(id_abs)
    c%label           = data(band)%label
    c%type            = cpar%cs_type(id_abs)
    c%class           = cpar%cs_class(id_abs)
    c%unit            = data(band)%unit
    c%nu_ref          = data(band)%bp(0)%p%nu_c
    c%cg_scale        = 1.d0
    c%myid            = cpar%myid_chain
    c%comm            = cpar%comm_chain
    c%numprocs        = cpar%numprocs_chain
    c%init_from_HDF   = cpar%cs_initHDF(id_abs)
    c%lmax_pre_lowl   = -1
    c%lmax_def        = -1
    c%nside_def       = 0
    c%fwhm_def        = 0.d0
    c%mono_prior_type = 'none'
    precond_type                = cpar%cg_precond

    call get_tokens(cpar%output_comps, ",", comp_label, n)
    c%output = .false.
    do i = 1, n
       if (trim(comp_label(i)) == trim(c%label) .or. trim(comp_label(i)) == 'all') then
          c%output = .true.
          exit
       end if
    end do

    !c%ref_band = band

    ! Set up conversion factor between RJ and native component unit
    select case (trim(c%unit))
    case ('uK_cmb')
       c%RJ2unit_ = data(band)%bp(0)%p%a2t
    case ('mK_cmb')
       c%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-3
    case ('K_cmb')
       c%RJ2unit_ = data(band)%bp(0)%p%a2t * 1d-6
    case ('MJy/sr') 
      !  c%RJ2unit_ = data(band)%bp(0)%p%a2t / data(band)%bp(0)%p%f2t
       c%RJ2unit_ = data(band)%bp(0)%p%a2f
    case ('uK_RJ') 
       c%RJ2unit_ = 1.d0
    case ('Kkm/s') 
       c%RJ2unit_ = 1.d0
       ! c%RJ2unit_ = 1.d0 / (c%nu_ref/c * 1d9)
!       write(*,*) 'Kkm/s not yet supported in md_mod'
!       call mpi_finalize(i)
!       stop
    case default
       call report_error('Unsupported unit: ' // trim(c%unit))
    end select

    ! Initialize comm_diffuse_comp_mod parameters
    c%pol      = .false.
    c%nside    = data(band)%map%info%nside
    c%lmax_amp = 1
    c%l_apod   = 2
    c%lmax_ind = 0
    c%cltype   = 'binned'
    c%nmaps    = 1
    allocate(c%lmax_ind_mix(3,1))
    c%lmax_ind_mix = 0

    !info          => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_amp, &
    !     & c%nmaps, c%pol)
    info          => comm_mapinfo(cpar%comm_chain, 128, c%lmax_amp, &
         & c%nmaps, c%pol)


    ! Diffuse preconditioner variables
    call add_to_npre(1,c%nside,1,1)

    ! Initialize amplitude and prior maps
    c%x   => comm_map(info)
    c%mu  => comm_map(info)
    c%ncr =  size(c%x%alm)
    do i = 0, c%x%info%nalm-1
       call c%x%info%i2lm(i,l,m)
       if (l == 0) then ! Monopole
          c%x%alm(i,1)  = sqrt(4.d0*pi) * def(1) / c%RJ2unit_(1)
          c%mu%alm(i,1) = sqrt(4.d0*pi) * mu(1)  / c%RJ2unit_(1)
       end if
       if (l == 1 .and. m == -1) then ! Y dipole
          c%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(3) / c%RJ2unit_(1)
          c%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(3)  / c%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  0) then ! Z dipole
          c%x%alm(i,1)  = sqrt(4.d0*pi/3.d0) * def(4) / c%RJ2unit_(1)
          c%mu%alm(i,1) = sqrt(4.d0*pi/3.d0) * mu(4)  / c%RJ2unit_(1)
       end if
       if (l == 1 .and. m ==  1) then ! X dipole
          c%x%alm(i,1)  = -sqrt(4.d0*pi/3.d0) * def(2) / c%RJ2unit_(1)
          c%mu%alm(i,1) = -sqrt(4.d0*pi/3.d0) * mu(2)  / c%RJ2unit_(1)
       end if
    end do


    ! Allocate mixing matrix
    c%ndet = maxval(data%ndet)
    allocate(c%F(numband,0:c%ndet), c%F_mean(numband,0:c%ndet,c%nmaps))
    allocate(c%F_null(numband,0:c%ndet))
    do i = 1, numband
       if (i == band) then
          info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
               & c%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
          c%F(i,0)%p      => comm_map(info)
          c%F(i,0)%p%map  = c%RJ2unit_(1)
          c%F(i,0)%p%alm  = c%RJ2unit_(1) * sqrt(4.d0*pi)
          c%F_null(i,:)   = .false.
          c%F_mean(i,:,:) = c%RJ2unit_(1)
          do j = 1, data(i)%ndet
             c%F(i,j)%p      => c%F(i,0)%p
          end do
       else
          do j = 0, data(i)%ndet
             c%F_null(i,j)   = .true.
             !c%F(i)%p%map  = 0.d0
             !c%F(i)%p%alm  = 0.d0
             c%F_mean(i,j,:) = 0.d0
          end do
       end if
    end do


    ! Initialize output beam
    c%B_out => comm_B_bl(cpar, c%x%info, 0, 0, fwhm=0.d0, init_realspace=.false.)

    ! Initialize power spectrum
    allocate(c%Cl)
    c%Cl%type   = 'binned'
    c%Cl%info   => c%x%info
    c%Cl%label  = 'md_'//trim(data(band)%label)
    c%Cl%lmax   = 1
    c%Cl%nmaps  = 1
    c%Cl%nspec  = 1
    c%Cl%outdir = 'none'
    allocate(c%Cl%Dl(0:1,1), c%Cl%sqrtS_mat(1,1,0:1))
    allocate(c%Cl%sqrtInvS_mat(1,1,0:1), c%Cl%S_mat(1,1,0:1))
    c%Cl%Dl(0,1) = 4.d0*pi      * rms(1)**2 / c%RJ2unit_(1)**2
    c%Cl%Dl(1,1) = 4.d0*pi/3.d0 * rms(2)**2 / c%RJ2unit_(1)**2
    c%Cl%S_mat(1,1,0:1)        = rms**2     / c%RJ2unit_(1)**2
    c%Cl%sqrtS_mat(1,1,0:1)    = rms        / c%RJ2unit_(1)
    if (rms(1) > 0.d0) c%Cl%sqrtInvS_mat(1,1,0) = 1.d0/rms(1) * c%RJ2unit_(1)
    if (rms(2) > 0.d0) c%Cl%sqrtInvS_mat(1,1,1) = 1.d0/rms(2) * c%RJ2unit_(1)
    
    ! Initialize md_mod specific parameters
    c%npar = 0

    ! Precompute mixmat integrator for each band
    allocate(c%F_int(3,numband,0:c%ndet))
    do k = 1, 3
       do i = 1, numband
          if (i == band) then
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (c%nu_ref(k) == c%nu_ref(k-1)) then
                      c%F_int(k,i,j)%p => c%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
!                write(*,*) 'mono didabled'
                c%F_int(k,i,j)%p => comm_F_line(c, data(i)%bp(j)%p, .true., 1.d0, -1)
             end do
          else
             do j = 0, data(i)%ndet
                if (k > 1) then
                   if (c%nu_ref(k) == c%nu_ref(k-1)) then
                      c%F_int(k,i,j)%p => c%F_int(k-1,i,j)%p
                      cycle
                   end if
                end if
!                                write(*,*) 'mono didabled'
                c%F_int(k,i,j)%p => comm_F_line(c, data(i)%bp(j)%p, .true., 0.d0, -1)
             end do
          end if
       end do
    end do

    ! Set up CG sampling groups                                                                                                                                             
    allocate(c%active_samp_group(cpar%cg_num_samp_groups))
    c%active_samp_group = .false.
    do i = 1, cpar%cg_num_samp_groups
       call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
       do j = 1, n
          if (trim(c%label) == trim(comp_label(j))) then
             c%active_samp_group(i) = .true.
             if (n == 1) c%cg_unique_sampgroup = i ! Dedicated sampling group for this component  
             exit
          end if
       end do
    end do

    ! Set up default values for prior sampling (to be potentially changed at end of init)  
    c%mono_from_prior=.false.
    c%mono_alm = 0.d0
    
  end function constructor_md

  ! Definition:
  !    SED  = delta_{band,ref_band}
  function evalSED_md(self, nu, band, pol, theta)
    class(comm_md_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_md

    integer(i4b) :: i, ind

    evalSED_md = 1.d0

  end function evalSED_md

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
          call initialize_md_comps%addComp(c)
       end if
       n = n+1
    end do
1   close(unit)
  
    if (n < numband .and. cpar%myid == 0) then
       write(*,'(a,i6)') ' | Warning: Number of channels without a monopole/dipole definition = ', numband-n
       stop
    end if

  end function initialize_md_comps
  
end module comm_md_comp_mod
