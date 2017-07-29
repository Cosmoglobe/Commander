module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_physdust_comp_mod
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
  use comm_chisq_mod
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
             call update_status(status, "init_done")
          case ("physdust")
             c => comm_physdust_comp(cpar, ncomp, i)
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
          c => initialize_template_comps(cpar, ncomp, i, n)
          !write(*,*) cpar%myid, ncomp, n
          ncomp = ncomp + n - 1
          call add_to_complist(c)
       case default
          call report_error("Unknown component class: "//trim(cpar%cs_class(i)))
       end select
    end do

    if (ncomp == 0) call report_error("Error: No signal components included in parameter file.")

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

!    call mpi_finalize(i)
!    stop
       
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

  subroutine sample_amps_by_CG(cpar, samp_group, handle)
    implicit none

    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: samp_group
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: stat, i
    real(dp)     :: Nscale = 1.d-4
    class(comm_comp), pointer :: c
    real(dp),           allocatable, dimension(:) :: rhs, x, mask

    allocate(x(ncr), mask(ncr))

    ! Set up component mask for current sample group
    c => compList
    do while (associated(c))
       call c%CG_mask(samp_group, mask)
       c => c%next()
    end do
    
    ! Solve the linear system
    call cr_computeRHS(cpar%operation, handle, mask, samp_group, rhs)
    call update_status(status, "init_precond1")
    call initPrecond(cpar%comm_chain)
    call update_status(status, "init_precond2")
    call solve_cr_eqn_by_CG(cpar, samp_group, x, rhs, stat)
    call cr_x2amp(samp_group, x)
    call update_status(status, "cr_end")
    deallocate(rhs,x)

  end subroutine sample_amps_by_CG

  subroutine initPrecond(comm)
    implicit none
    integer(i4b), intent(in) :: comm
    call initDiffPrecond(comm)
    call initPtsrcPrecond(comm)
    call initTemplatePrecond(comm)
  end subroutine initPrecond

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

    if (cpar%init_samp <= 0 .or. trim(cpar%init_chain_prefix) == 'none') return

    ! Open HDF file
    call int2string(cpar%mychain,   ctext)
    call int2string(cpar%init_samp, itext)
    if (trim(cpar%chain_prefix) == trim(cpar%init_chain_prefix)) then
       chainfile = trim(adjustl(cpar%outdir)) // '/' // trim(adjustl(cpar%chain_prefix)) // &
            & '_c' // trim(adjustl(ctext)) // '.h5'
    else
       chainfile = trim(adjustl(cpar%init_chain_prefix))
    end if
    call open_hdf_file(chainfile, file, 'r')
    
    ! Initialize instrumental parameters
    call update_status(status, "init_chain_inst")
    if (cpar%cs_init_inst_hdf) then
       do i = 1, numband
          if (cpar%ignore_gain_bp) then
             data(i)%gain     = 1.d0
             data(i)%bp%delta = 0.d0
          else
             call read_hdf(file, trim(adjustl(itext))//'/gain/'//trim(adjustl(data(i)%label)), &
                  & data(i)%gain)
             call read_hdf(file, trim(adjustl(itext))//'/bandpass/'//trim(adjustl(data(i)%label)), &
                  & data(i)%bp%delta)
          end if
       end do
    end if
    
    ! Initialize component parameters
    c   => compList
    do while (associated(c))
       if (.not. c%init_from_HDF) then
          c => c%next()
          cycle
       end if
       call update_status(status, "init_chain_"//trim(c%label))
       call c%initHDF(cpar, file, trim(adjustl(itext))//'/')
       c => c%next()
    end do

    ! Close HDF file
    call close_hdf_file(file)
    
  end subroutine initialize_from_chain


  subroutine sample_partialsky_tempamps(cpar, handle)
    implicit none

    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: i, ierr
    logical(lgt) :: skip
    real(dp)     :: vals(2), vals2(2), mu, sigma, amp, mu_p, sigma_p
    class(comm_map),           pointer :: res, invN_T
    class(comm_comp),          pointer :: c
    class(comm_template_comp), pointer :: pt
    
    ! Compute predicted signal for this band
    c => compList
    do while (associated(c))
       skip = .true.
       select type (c)
       class is (comm_template_comp)
          pt  => c
          if (associated(pt%mask)) skip = .false.
       end select
       if (skip) then
          c => c%next()
          cycle
       end if

       ! Get residual map
       res      => compute_residual(pt%band, exclude_comps=[pt%label]) 
       invN_T => comm_map(pt%T)
       call data(pt%band)%N%invN(invN_T)

       ! Compute mean and variance
       vals(1) = sum(invN_T%map * res%map  * pt%mask%map)
       vals(2) = sum(invN_T%map * pt%T%map * pt%mask%map)
       call mpi_reduce(vals, vals2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, res%info%comm, ierr)       

       if (res%info%myid == 0) then
          ! Compute mean and RMS from likelihood term
          mu    = vals2(1) / vals2(2)
          sigma = sqrt(1.d0/vals2(2))

          ! Add prior
          mu_p    = pt%P(1)
          sigma_p = pt%P(2)
          mu      = (mu*sigma_p**2 + mu_p * sigma**2) / (sigma_p**2 + sigma**2)
          sigma   = sqrt(sigma**2 * sigma_p**2 / (sigma**2 + sigma_p**2))

          if (trim(cpar%operation) == 'sample') then
             amp = mu + sigma * rand_gauss(handle)
          else
             amp = mu
          end if

          ! Update template amplitude in main structure
          pt%x       = amp  ! Amplitude
          pt%P_cg(1) = amp  ! CG prior mean
          
       end if

       deallocate(res, invN_T)

       c => c%next()
    end do

  end subroutine sample_partialsky_tempamps
  
end module comm_signal_mod
