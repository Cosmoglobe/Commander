module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_spindust_comp_mod
  use comm_MBB_comp_mod
  use comm_freefree_comp_mod
  use comm_line_comp_mod
  use comm_md_comp_mod
  use comm_template_comp_mod
  use comm_ptsrc_comp_mod
  use comm_cr_mod
  use comm_cr_utils
  use comm_hdf_mod
  use comm_data_mod
  implicit none

contains

  subroutine initialize_signal_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    integer(i4b) :: i, n
    class(comm_comp), pointer :: c
    
    ncomp = 0
    do i = 1, cpar%cs_ncomp_tot
       if (.not. cpar%cs_include(i)) cycle
       ncomp = ncomp + 1
       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') '  Initializing component ', i, ' : ', trim(cpar%cs_label(i))
       call update_status(status, "init_"//trim(cpar%cs_label(i)))

       ! Initialize object
       select case (trim(cpar%cs_class(i)))
       case ("diffuse")
          ! Diffuse components
          select case (trim(cpar%cs_type(i)))
          case ("cmb")
             c => comm_cmb_comp(cpar, ncomp, i)
          case ("power_law")
             c => comm_powlaw_comp(cpar, ncomp, i)
          case ("spindust")
             c => comm_spindust_comp(cpar, ncomp, i)
          case ("MBB")
             c => comm_MBB_comp(cpar, ncomp, i)
          case ("freefree")
             c => comm_freefree_comp(cpar, ncomp, i)
          case ("line")
             c => comm_line_comp(cpar, ncomp, i)
          case ("md")
             c => initialize_md_comps(cpar, ncomp, i, n)
             ncomp = ncomp + n - 1
          case default
             call report_error("Unknown component type: "//trim(cpar%cs_type(i)))
          end select
          call add_to_complist(c)
       case ("ptsrc")
          c => comm_ptsrc_comp(cpar, ncomp, i)
          call add_to_complist(c)
          !call c%dumpFITS('test', cpar%outdir)
          !call mpi_finalize(i)
          !stop
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

    allocate(x(ncr))
    call cr_computeRHS(handle, rhs)
    call update_status(status, "init_precond1")
    call initDiffPrecond(cpar%comm_chain)
    call initPtsrcPrecond(cpar%comm_chain)
    call update_status(status, "init_precond2")
    call solve_cr_eqn_by_CG(cpar, x, rhs, stat)
    call cr_x2amp(x)
    deallocate(rhs,x)

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

  subroutine initialize_from_chain(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    integer(i4b)              :: i
    character(len=4)          :: ctext
    character(len=6)          :: itext
    character(len=512)        :: chainfile, hdfpath
    class(comm_comp), pointer :: c
    type(hdf_file) :: file

    if (cpar%init_samp <= 0) return

    ! Open HDF file
    call int2string(cpar%mychain,   ctext)
    call int2string(cpar%init_samp, itext)
    chainfile = trim(adjustl(cpar%outdir)) // '/' // trim(adjustl(cpar%chain_prefix)) // &
         & '_c' // trim(adjustl(ctext)) // '.h5'
    call open_hdf_file(chainfile, file, 'r')
    
    ! Initialize instrumental parameters
    do i = 1, numband
       call read_hdf(file, trim(adjustl(itext))//'/gain/'//trim(adjustl(data(i)%label)), &
            & data(i)%gain)
       call read_hdf(file, trim(adjustl(itext))//'/bandpass/'//trim(adjustl(data(i)%label)), &
            & data(i)%bp%delta)
    end do
    
    ! Initialize component parameters
    c   => compList
    do while (associated(c))
       call c%initHDF(cpar, file, trim(adjustl(itext))//'/')
       c => c%next()
    end do

    ! Close HDF file
    call close_hdf_file(file)
    
  end subroutine initialize_from_chain
  
end module comm_signal_mod
