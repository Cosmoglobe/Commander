module comm_signal_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_cmb_comp_mod
  use comm_powlaw_comp_mod
  use comm_physdust_comp_mod
  use comm_spindust_comp_mod
  use comm_spindust2_comp_mod
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
          case ("spindust2")
             c => comm_spindust2_comp(cpar, ncomp, i)
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

  subroutine sample_amps_by_CG(cpar, samp_group, handle, handle_noise)
    implicit none

    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: samp_group
    type(planck_rng),  intent(inout) :: handle, handle_noise

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
    call cr_computeRHS(cpar%operation, cpar%resamp_CMB, handle, handle_noise, mask, samp_group, rhs)
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

  subroutine initialize_from_chain(cpar, handle, init_samp, init_from_output)
    implicit none
    type(comm_params), intent(in)           :: cpar
    type(planck_rng),  intent(inout)        :: handle
    integer(i4b),      intent(in), optional :: init_samp
    logical(lgt),      intent(in), optional :: init_from_output

    integer(i4b)              :: i, j, ext(2), initsamp
    character(len=4)          :: ctext
    character(len=6)          :: itext
    character(len=512)        :: chainfile, hdfpath
    class(comm_comp), pointer :: c
    type(hdf_file) :: file
    class(comm_N),      pointer :: N
    class(comm_map),    pointer :: rms
    real(dp), allocatable, dimension(:,:) :: bp_delta, regnoise

    if ((cpar%init_samp <= 0 .or. trim(cpar%init_chain_prefix) == 'none') .and. &
         & .not. present(init_from_output)) return

    ! Open HDF file
    call int2string(cpar%mychain,   ctext)
    initsamp = cpar%init_samp; if (present(init_samp)) initsamp = init_samp
    call int2string(initsamp, itext)
    if (present(init_from_output)) then
       chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
            & '_c' // trim(adjustl(ctext)) // '.h5'
    else
       chainfile = trim(adjustl(cpar%init_chain_prefix))
    end if
    !write(*,*) 'chainfile', trim(chainfile)
    call open_hdf_file(chainfile, file, 'r')
    
    if (cpar%resamp_CMB .and. present(init_from_output)) then
    ! Initialize CMB component parameters; only once before starting Gibbs
       c   => compList
       do while (associated(c))
          if (trim(c%type) /= 'cmb') then
             c => c%next()
             cycle
          end if
          call update_status(status, "init_chain_"//trim(c%label))
          if (cpar%myid == 0) write(*,*) ' Initializing from chain = ', trim(c%label)
          call c%initHDF(cpar, file, trim(adjustl(itext))//'/')
          c => c%next()
       end do
       call close_hdf_file(file)
       return
    end if

    ! Initialize instrumental parameters
    call update_status(status, "init_chain_inst")
    if (cpar%cs_init_inst_hdf) then
       do i = 1, numband
          if (cpar%ignore_gain_bp) then
             data(i)%gain     = 1.d0
             data(i)%bp(0)%p%delta = 0.d0
          else
             call read_hdf(file, trim(adjustl(itext))//'/gain/'//trim(adjustl(data(i)%label)), &
                  & data(i)%gain)
  
             call get_size_hdf(file, trim(adjustl(itext))//'/bandpass/'//&
                  & trim(adjustl(data(i)%label)), ext)
             if (data(i)%ndet > ext(1)-1) then
                write(*,*) 'Error -- init HDF file does not contain enough bandpass information'
                stop
             end if
             allocate(bp_delta(0:ext(1)-1,ext(2)))
             call read_hdf(file, trim(adjustl(itext))//'/bandpass/'//trim(adjustl(data(i)%label)), &
                  & bp_delta)
             do j = 0, data(i)%ndet
                data(i)%bp(j)%p%delta = bp_delta(j,:)
             end do
             deallocate(bp_delta)
          end if
       end do
    end if
    
    ! Initialize component parameters
    c   => compList
    do while (associated(c))
       if (.not. c%init_from_HDF .or. (present(init_samp) .and. trim(c%type) == 'cmb')) then
          c => c%next()
          cycle
       end if
       call update_status(status, "init_chain_"//trim(c%label))
       if (cpar%myid == 0) write(*,*) ' Initializing from chain = ', trim(c%label)
       call c%initHDF(cpar, file, trim(adjustl(itext))//'/')
       c => c%next()
    end do

    ! Initialize TOD parameters
    if (cpar%enable_TOD_analysis) then
       do i = 1, numband  
          if (trim(data(i)%tod_type) == 'none') cycle
          if (.not. data(i)%tod%init_from_HDF)     cycle
          if (cpar%myid == 0) write(*,*) ' Initializing TOD par from chain = ', trim(data(i)%tod%freq)
          N => data(i)%N
          rms => comm_map(data(i)%info)
          select type (N)
          class is (comm_N_rms)
             call data(i)%tod%initHDF(file, initsamp, data(i)%map, rms)
          end select

          ! Update rms and data maps
          allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
          if (associated(data(i)%procmask)) then
             call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
          else
             call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, map=rms)
          end if
          if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
          data(i)%map%map = data(i)%map%map + regnoise         ! Add regularization noise
          data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask
          deallocate(regnoise)
          call rms%dealloc
       end do
    else if (cpar%resamp_CMB) then
       do i = 1, numband  
          cycle

          if (trim(data(i)%tod_type) == 'none') cycle
          !if (.not. data(i)%tod%init_from_HDF)  cycle
          if (cpar%myid == 0) write(*,*) ' Initializing map and rms from chain = ', trim(data(i)%label)

          hdfpath =  trim(adjustl(itext))//'/tod/'//trim(adjustl(data(i)%label))//'/'
          rms     => comm_map(data(i)%info)
          N       => data(i)%N
          call data(i)%map%readMapFromHDF(file, trim(adjustl(hdfpath))//'map')
          select type (N)
          class is (comm_N_rms)
             call rms%readMapFromHDF(file, trim(adjustl(hdfpath))//'rms')
          end select


          ! Update rms and data maps
          allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
          if (associated(data(i)%procmask)) then
             call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
          else
             call data(i)%N%update_N(data(i)%info, handle, data(i)%mask, regnoise, map=rms)
          end if
          if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
          data(i)%map%map = data(i)%map%map + regnoise         ! Add regularization noise
          data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask

!!$          call data(i)%map%writeFITS('map.fits')
!!$          call rms%writeFITS('rms.fits')
!!$          call mpi_finalize(j)
!!$          stop

          deallocate(regnoise)
          call rms%dealloc

       end do

    end if


    ! Close HDF file
    call close_hdf_file(file)
    
  end subroutine initialize_from_chain


  subroutine sample_powspec(handle, ok)
    implicit none
    type(planck_rng),  intent(inout) :: handle
    logical(lgt),      intent(out)   :: ok

    class(comm_comp), pointer :: c

    ok = .true.
    c  => compList
    do while (associated(c))
       if (trim(c%type) == 'md') then
          c => c%next()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          call c%Cl%sampleCls(c%x, handle, ok)
       end select
       c => c%next()
    end do

  end subroutine sample_powspec


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


  subroutine synchronize_bp_delta(initHDF)
    implicit none
    logical(lgt), intent(in) :: initHDF

    integer(i4b) :: i, j, ndet
    real(dp)     :: mu

    do i = 1, numband
       if (trim(data(i)%tod_type) == 'none') cycle
       ndet = data(i)%ndet
       do j = 1, ndet
          data(i)%bp(j)%p%delta = data(i)%tod%bp_delta(j,:)
       end do
       if (initHDF) then
          data(i)%bp(0)%p%delta = data(i)%tod%bp_delta(0,:)
       else
          do j = 1, size(data(i)%bp(0)%p%delta)
             mu = mean(data(i)%tod%bp_delta(1:data(i)%ndet,j))
             data(i)%tod%bp_delta(1:ndet,j) = data(i)%tod%bp_delta(1:ndet,j) - &
                  & mu + data(i)%bp(0)%p%delta(j)
          end do
       end if
    end do

  end subroutine synchronize_bp_delta
  

  subroutine sample_joint_alm_Cl(handle)
    implicit none
    type(planck_rng), intent(inout) :: handle

    integer(i4b) :: i, j, n, bin, pos, ierr
    logical(lgt) :: posdef, accept
    real(dp)     :: chisq_old, chisq_prop, U
    real(dp), allocatable, dimension(:)   :: Dl_old, Dl_prop, eta
    real(dp), allocatable, dimension(:,:) :: alm_prop, alm_old
    class(comm_comp), pointer :: c
    class(comm_map),  pointer :: res
    
    ! Initialize residual maps
    chisq_old = 0.d0
    do i = 1, numband
       res             => compute_residual(i)
!!$       write(*,*) sum(abs(res%map))
!!$       call mpi_finalize(ierr)
!!$       stop

       data(i)%res%map =  res%map
       chisq_old       =  chisq_old + data(i)%chisq()
!!$       write(*,*) chisq_old
!!$       call mpi_finalize(ierr)
!!$       stop
       call res%dealloc()
       nullify(res)
    end do

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)

          ! Add CMB 
          do i = 1, numband
             data(i)%c_old     => comm_map(data(i)%info)
             data(i)%c_prop    => comm_map(data(i)%info)
             data(i)%c_old%map =  c%getBand(i)
          end do

          allocate(alm_prop(0:c%x%info%nalm-1,0:c%x%info%nmaps))
          allocate(alm_old(0:c%x%info%nalm-1,0:c%x%info%nmaps))
          do bin = 1, c%Cl%nbin2
             if (trim(c%Cl%type) /= 'binned') cycle
             if (trim(c%Cl%bins2(bin)%stat) /= 'M') cycle

             n = c%Cl%bins2(bin)%ntot
if (c%x%info%myid ==0) write(*,*) bin, n
             pos = 0
             allocate(Dl_old(n), Dl_prop(n), eta(n))
             call c%Cl%set_Dl_bin(c%Cl%bins2(bin), Dl_old, pos, .false.)
             if (c%x%info%myid ==0) then   
                do i = 1, n
                   eta(i) = rand_gauss(handle)
                end do
                Dl_prop = Dl_old + matmul(c%Cl%bins2(bin)%M_prop, eta)
                posdef = c%Cl%check_posdef(bin, Dl_prop)
             end if
             call mpi_bcast(posdef, 1, MPI_LOGICAL, 0, c%x%info%comm, ierr)
             if (c%x%info%myid ==0) write(*,*) 'posdef', bin, c%Cl%bins2(bin)%lmin, posdef
             if (posdef) then
                call mpi_bcast(Dl_prop, n, MPI_DOUBLE_PRECISION, 0, c%x%info%comm, ierr)

                ! Compute rescaled alms
                alm_old  = c%x%alm
                alm_prop = c%x%alm
                pos      = 0
                call c%Cl%sqrtInvS(alm=alm_prop, info=c%x%info)
                call c%Cl%set_Dl_bin(c%Cl%bins2(bin), Dl_prop, pos, .true.)
                call c%Cl%updateS
                call c%Cl%sqrtS(alm=alm_prop, info=c%x%info)
                c%x%alm = alm_prop

                ! Compute proposal chisquare
                chisq_prop = 0.d0
                do i = 1, numband
                   data(i)%c_prop%map =  c%getBand(i)
                   data(i)%res%map    =  data(i)%res%map + data(i)%c_old%map - data(i)%c_prop%map
                   chisq_prop         =  chisq_prop      + data(i)%chisq()
                end do

                ! Apply Metropolis rule, and update data structures
                if (c%x%info%myid == 0) then
                   accept = rand_uni(handle) < exp(-0.5d0*(chisq_prop-chisq_old))
                   write(*,*) 'Old   = ', bin, real(Dl_old,sp)
                   write(*,*) 'Prop  = ', bin, real(Dl_prop,sp)
                   write(*,*) 'chisq = ', chisq_old, chisq_prop, accept
                end if
!!$                call mpi_finalize(ierr)
!!$                stop

                call mpi_bcast(accept, 1, MPI_LOGICAL, 0, c%x%info%comm, ierr)
                if (accept) then
                   chisq_old = chisq_prop
                   do i = 1, numband
                      data(i)%c_old%map = data(i)%c_prop%map 
                   end do
                else
                   c%x%alm = alm_old
                   pos     = 0
                   call c%Cl%set_Dl_bin(c%Cl%bins2(bin), Dl_old, pos, .true.)
                   call c%Cl%updateS
                   do i = 1, numband
                      data(i)%res%map   = data(i)%res%map - data(i)%c_old%map + data(i)%c_prop%map
                   end do
                end if

             end if
             deallocate(Dl_old, Dl_prop, eta)
          end do
          deallocate(alm_prop, alm_old)

          do i = 1, numband
             call data(i)%c_old%dealloc()
             call data(i)%c_prop%dealloc() 
          end do

       end select
       c => c%next()
    end do

  end subroutine sample_joint_alm_Cl


end module comm_signal_mod
