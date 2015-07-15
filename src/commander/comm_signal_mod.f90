module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_template_comp_mod
  use comm_ptsrc_comp_mod
  use comm_cr_mod
  use comm_cr_utils
  implicit none

contains

  subroutine initialize_signal_mod(cpar)
    implicit none

    type(comm_params), intent(in) :: cpar

    integer(i4b) :: i
    class(comm_comp), pointer :: c
    
    ncomp = cpar%cs_ncomp
    do i = 1, ncomp
       ! Initialize object
       select case (trim(cpar%cs_class(i)))
       case ("diffuse")
          ! Diffuse components
          select case (trim(cpar%cs_type(i)))
          case ("cmb")
             c => comm_cmb_comp(cpar, i)
          case ("power_law")
             c => comm_powlaw_comp(cpar, i)
          case default
             call report_error("Unknown component type: "//trim(cpar%cs_type(i)))
          end select
          call add_to_complist(c)
       case ("ptsrc")
          call initialize_ptsrc_comp(cpar, i)
       case ("template")
          ! Point source components
          select case (trim(cpar%cs_type(i)))
          case ("monopole")
!             c => comm_cmb_comp(cpar, i)
          case ("dipole")
!             c => comm_powlaw_comp(cpar, i)
          case ("rel_quadrupole")
!             c => comm_powlaw_comp(cpar, i)
          case ("file")
!             c => comm_powlaw_comp(cpar, i)
          case default
             call report_error("Unknown component type: "//trim(cpar%cs_type(i)))
          end select
       case default
          call report_error("Unknown component class: "//trim(cpar%cs_class(i)))
       end select

    end do

    ! Compute position and length of each component in parameter array
    allocate(ind_comp(ncomp,3))
    ncr = 0
    i   = 1
    ind_comp(i,:) = [1,0]
    c => compList
    do while (associated(c))
       ind_comp(i,1) = ncr+1
       ind_comp(i,2) = c%ncr
       ind_comp(i,3) = c%nmaps
       ncr           = ncr + c%ncr
       i             = i+1
       c             => c%next()
    end do

  end subroutine initialize_signal_mod

  subroutine dump_components(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(i4b) :: i, unit
    class(comm_comp), pointer :: c

    unit = getlun()
    
    open(unit, file=trim(filename))
    c => compList
    do while (associated(c))
       write(unit,*) '# Component = ', trim(c%label)
       call c%dumpSED(unit)
       write(unit,*)
       c => c%next()
    end do
    close(unit)
    
  end subroutine dump_components

  subroutine sample_amps_by_CG(cpar, handle)
    implicit none

    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: stat, i
    real(dp)     :: Nscale = 1.d-4
    real(dp),           allocatable, dimension(:) :: rhs, x
    class(precondDiff), pointer                   :: P

    allocate(x(ncr))
    call cr_computeRHS(handle, rhs)
    P => precondDiff(cpar%comm_chain, Nscale)
    call solve_cr_eqn_by_CG(cpar, cr_matmulA, cr_invM, x, rhs, stat, P)
    call cr_x2amp(x)
    deallocate(rhs,x,P)

  end subroutine sample_amps_by_CG

  subroutine add_to_complist(c)
    implicit none
    class(comm_comp), pointer :: c
    
    if (.not. associated(compList)) then
       compList => c
    else
       call compList%add(c)
    end if
  end subroutine add_to_complist

  subroutine initialize_ptsrc_comp(cpar, id)
    implicit none
    class(comm_params), intent(in) :: cpar
    integer(i4b),       intent(in) :: id

    integer(i4b)        :: unit, i, npar, nmaps
    real(dp)            :: glon, glat, nu_ref
    logical(lgt)        :: pol
    character(len=1024) :: line
    character(len=128)  :: id_ptsrc
    class(comm_comp),       pointer :: c
    class(comm_ptsrc_comp), pointer :: P_ref
    real(dp), allocatable, dimension(:)   :: amp
    real(dp), allocatable, dimension(:,:) :: beta

    unit = getlun()

    nmaps = 1; if (cpar%cs_polarization(id)) nmaps = 3
    select case (trim(cpar%cs_type(id)))
    case ("radio")
       npar = 2
    case ("fir")
       npar = 2
    case ("sz")
       npar = 0
    end select
    allocate(amp(nmaps), beta(npar,nmaps))
    
    ! Initialize point sources based on catalog information
    nullify(P_ref)
    open(unit,file=trim(cpar%datadir) // '/' // trim(cpar%cs_catalog(id)),recl=1024)
    do while (.true.)
       read(unit,'(a)',end=1) line
       line = trim(line)
       if (line(1:1) == '#' .or. trim(line) == '') cycle
       write(*,*) trim(line)
       read(line,*) glon, glat, amp, beta, id_ptsrc
       c => comm_ptsrc_comp(cpar, id, glon, glat, amp, beta, id_ptsrc, P_ref)
       call add_to_complist(c)
       select type (c)
       class is (comm_ptsrc_comp)
          P_ref => c
       end select
    end do 
1   close(unit)

    call mpi_finalize(i)
    stop

  end subroutine initialize_ptsrc_comp
  
end module comm_signal_mod
