module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_template_comp_mod
  use comm_cr_mod
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
       select case (trim(cpar%cs_type(i)))
       case ("cmb")
          c => comm_cmb_comp(cpar, i)
       case ("power_law")
          c => comm_powlaw_comp(cpar, i)
       case default
          call report_error("Unknown component type: "//trim(cpar%cs_type(i)))
       end select

       ! Add object to list
       if (i == 1) then
          compList => c
       else
          call compList%setNext(c)
       end if
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

    call mpi_finalize(i)
    stop

    
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
    real(dp)     :: Nscale = 1.d0
    real(dp), allocatable, dimension(:) :: rhs, x

    allocate(x(ncr))
    call cr_computeRHS(handle, rhs)
    call solve_cr_eqn_by_CG(cpar, cr_matmulA, cr_invM, x, rhs, stat, Nscale)
    call cr_x2amp(x)
    deallocate(rhs,x)

  end subroutine sample_amps_by_CG

  
end module comm_signal_mod
