module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_template_comp_mod
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

end module comm_signal_mod
