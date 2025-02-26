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
module comm_signal_mod
  use comm_chisq_mod
  use comm_cr_mod
  use comm_cmb_comp_mod
  use comm_cmb_relquad_comp_mod
  use comm_mbb_comp_mod
  use comm_mbbtab_comp_mod
  use comm_powlaw_comp_mod
  use comm_curvature_comp_mod
  use comm_spindust_comp_mod
  use comm_spindust2_comp_mod
  use comm_md_comp_mod
  use comm_line_comp_mod
  use comm_freefree_comp_mod
  use comm_freefreeEM_comp_mod
  use comm_exp_comp_mod
  use comm_physdust_comp_mod
  use comm_pah_comp_mod
  use comm_powlaw_break_comp_mod
  use comm_ame_lognormal_mod
  implicit none

  integer(i4b), parameter, private :: MAXSAMPGROUP     = 100


contains

  subroutine initialize_signal_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    integer(i4b) :: i, n
    class(comm_comp), pointer :: c => null()
    class(comm_comp), pointer :: c_two => null()
    logical(lgt) :: prior_exists

    ncomp = 0
    do i = 1, cpar%cs_ncomp_tot
       if (.not. cpar%cs_include(i)) cycle
       ncomp = ncomp + 1
       if (cpar%myid == 0 .and. cpar%verbosity > 0) &
            & write(*,fmt='(a,i5,a,a)') ' |  Initializing component ', i, ' : ', trim(cpar%cs_label(i))
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
          case ("exponential")
             c => comm_exp_comp(cpar, ncomp, i)
          case ("power_law_break")
             c => comm_powlaw_break_comp(cpar, ncomp, i)
          case ("curvature") 
             c => comm_curvature_comp(cpar, ncomp, i)
          case ("physdust")
             !if (cpar%myid == 0) write(*,*) 'dd beg', i, trim(cpar%cs_type(i))
             c => comm_physdust_comp(cpar, ncomp, i)
             !if (cpar%myid == 0) write(*,*) 'dd end', i, trim(cpar%cs_type(i))
          case ("spindust")
             c => comm_spindust_comp(cpar, ncomp, i)
          case ("spindust2")
             c => comm_spindust2_comp(cpar, ncomp, i)
          case ("lognormal")
             c => comm_ame_lognormal_comp(cpar, ncomp, i)
          case ("MBB")
             c => comm_MBB_comp(cpar, ncomp, i)
          case ("MBBtab")
             c => comm_MBBtab_comp(cpar, ncomp, i)
          case ("freefree")
             c => comm_freefree_comp(cpar, ncomp, i)
          case ("freefreeEM")
             c => comm_freefreeEM_comp(cpar, ncomp, i)
          case ("line")
             c => comm_line_comp(cpar, ncomp, i)
          case ("md")
             c => initialize_md_comps(cpar, ncomp, i, n)
             ncomp = ncomp + n - 1
          case ("pah")
             c => comm_pah_comp(cpar, ncomp, i)
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
          select case (trim(cpar%cs_type(i)))
          case ("template")
             c => initialize_template_comps(cpar, ncomp, i, n)
          !write(*,*) cpar%myid, ncomp, n
             ncomp = ncomp + n - 1
          case ("cmb_relquad")
             c => comm_cmb_relquad_comp(cpar, ncomp, i)
          end select
          call add_to_complist(c)
       case default
          call report_error("Unknown component class: "//trim(cpar%cs_class(i)))
       end select
       !if (cpar%myid == 0 .and. cpar%verbosity > 0) &
       !     & write(*,fmt='(a,i5,a,a)') ' |  finished  component ', i, ' : ', trim(cpar%cs_label(i))
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
       c             => c%nextComp()
    end do

       
    ! go through compList and check if any diffuse component is using a band monopole
    ! as the zero-level prior

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%mono_prior_type) == 'bandmono') then
             prior_exists=.false.
             c_two => compList
             do while (associated(c_two))
                select type (c_two)
                class is (comm_md_comp)
                   if (trim(c%mono_prior_band)==trim(c_two%label)) then
                      if (c_two%mono_from_prior) then
                         !Error, band already prior for another component
                         call report_error("Component '"//trim(c%label)//"'. Band monopole '"//trim(c_two%label)//"' already in use as a zero-level prior of another component")
                      end if
                      c_two%mono_from_prior = .true.
                      prior_exists = .true.
                   end if
                end select
                c_two => c_two%nextComp()
             end do

             if (.not. prior_exists) then
                !Error, could not find band monopole used as prior
                call report_error("Could not find band monopole '"//trim(c%mono_prior_band)//"' for zero-level prior of component "//trim(c%label))

             end if

          end if
       end select
       c => c%nextComp()
    end do

  end subroutine initialize_signal_mod

  subroutine dump_components(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(i4b) :: i, unit
    class(comm_comp), pointer :: c => null()

    unit = getlun()
    
    open(unit, file=trim(filename))
    c => compList
    do while (associated(c))
       write(unit,*) '# Component = ', trim(c%label)
       call c%dumpSED(unit)
       write(unit,*)
       c => c%nextComp()
    end do
    close(unit)
    
  end subroutine dump_components

  subroutine sample_amps_by_CG(cpar, samp_group, handle, handle_noise)
    implicit none

    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: samp_group
    type(planck_rng),  intent(inout) :: handle, handle_noise

    integer(i4b) :: stat, i, j, l, m, nactive
    real(dp)     :: Nscale = 1.d-4
    class(comm_comp), pointer :: c => null()
    character(len=32) :: cr_active_bands(100)
    real(dp),           allocatable, dimension(:) :: rhs, x, mask
    class(comm_map),     pointer :: res  => null()


    allocate(x(ncr), mask(ncr))

    ! Set up component mask for current sample group
    c => compList
    do while (associated(c))
       call c%CG_mask(samp_group, mask)
       c => c%nextComp()
    end do

    ! Activate requested frequency bands
    if(trim(cpar%cg_samp_group_bands(samp_group)) == 'all') then
      data%cr_active = .true.
    else
      call get_tokens(cpar%cg_samp_group_bands(samp_group), ',', cr_active_bands, nactive)
      data%cr_active = .false.
      do i = 1, nactive
         do j = 1, numband
            if (trim(cr_active_bands(i)) == trim(data(j)%label)) then
               data(j)%cr_active = .true.
               exit
            end if
         end do
      end do
    end if
   
    ! Sample point source amplitudes
    c => compList 
    do while (associated(c))
       select type (c)
       class is (comm_ptsrc_comp)
          if(c%precomputed_amps .and. c%active_samp_group(samp_group)) then
             ! Initialize residual maps
             do l = 1, numband
                res             => compute_residual(l)
                data(l)%res%map =  res%map
                call res%dealloc(); deallocate(res)
                nullify(res)
             end do
             ! Perform sampling
             call c%samplePtsrcAmp(cpar, handle, samp_group)
             return
          end if
       end select
       c => c%nextComp()
    end do
    
    ! If mono-/dipole are sampled, check if they are priors for a component zero-level
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_md_comp)
          if (c%active_samp_group(samp_group)) then
             if (c%mono_from_prior) then
                do i = 0, c%x%info%nalm-1
                   call c%x%info%i2lm(i,l,m)
                   if (l == 0) then ! save the monopole value
                      c%mono_alm = c%x%alm(i,1)
                   end if
                end do
             end if
          end if
       end select
       c => c%nextComp()
    end do

   
    ! Solve the linear system
    call cr_computeRHS(cpar%operation, cpar%resamp_CMB, cpar%only_pol,&
         & handle, handle_noise, mask, samp_group, rhs)
    call update_status(status, "init_precond1")
    call initPrecond(cpar%comm_chain)
    call update_status(status, "init_precond2")
    call solve_cr_eqn_by_CG(cpar, samp_group, x, rhs, stat)
    call cr_x2amp(samp_group, x)
    call update_status(status, "cr_end")
    deallocate(rhs,x)

    ! Apply monopole priors for active diffuse components
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%active_samp_group(samp_group)) call c%applyMonoDipolePrior(handle)
       end select
       c => c%nextComp()
    end do

    ! If mono-/dipole components is a zero-level prior, revert back to pre-sampling value if it has been sampled in the current CG group
    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_md_comp)
          if (c%active_samp_group(samp_group)) then
             if (c%mono_from_prior) then
                do i = 0, c%x%info%nalm-1
                   call c%x%info%i2lm(i,l,m)
                   if (l == 0) then ! monopole

                      write(*,fmt='(a)') " |  Band monopole of '"//&
                           & trim(c%label)//"' used as zero-level prior"
                      write(*,fmt='(a,f14.3)') " |     Revert back to pre-CG value: ",&
                           & c%mono_alm/sqrt(4.d0*pi)
                      write(*,fmt='(a,f14.3,a)') " |     (Sampled value in CG: ",&
                           & c%x%alm(i,1)/sqrt(4.d0*pi)," )"

                      c%x%alm(i,1) = c%mono_alm  ! revert to pre-CG search value 
                      !monopole in alm_uKRJ = mu_in_alm_uK_RJ + rms_in_alm_uKRJ * rand_gauss
                      !c%x%alm(i,1) = c%mu%alm(i,1) + sqrt(c%Cl%Dl(0,1))*rand_gauss(handle) 

                   end if
                end do
             end if
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine sample_amps_by_CG


  subroutine sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups)
    !
    !
    !  Convenience function for performing amplitude sampling over
    !  all sampling groups
    !
    !
    implicit none

    type(comm_params), intent(in)            :: cpar
    type(planck_rng),  intent(inout)         :: handle, handle_noise
    character(len=512), intent(in), optional :: cg_groups


    integer(i4b)                          :: samp_group, i, n
    integer(i4b), dimension(MAXSAMPGROUP) :: group_inds
    character(len=3) :: toks(MAXSAMPGROUP)


    if (present(cg_groups)) then
      group_inds = 0
      call get_tokens(cg_groups, ',', toks, n)
      do i = 1, n
        read(toks(i), *) group_inds(i)
      end do
    else
      group_inds = [(samp_group, samp_group=1, MAXSAMPGROUP)]
    end if


    call timer%start(TOT_AMPSAMP)
    do samp_group = 1, cpar%cg_num_user_samp_groups
       if (findloc(group_inds, samp_group, dim=1) == 0) cycle
       if (cpar%myid_chain == 0) then
          write(*,fmt='(a,i4,a,i4,a,i4,a,a)') ' |  Chain = ', cpar%mychain, &
          & ' -- CG sample group = ', samp_group, ' of ', cpar%cg_num_user_samp_groups, ': ', &
          & trim(cpar%cg_samp_group(samp_group))
       end if
       call sample_amps_by_CG(cpar, samp_group, handle, handle_noise)

       if (trim(cpar%cmb_dipole_prior_mask) /= 'none') call apply_cmb_dipole_prior(cpar, handle)

    end do
    call timer%stop(TOT_AMPSAMP)

  end subroutine sample_all_amps_by_CG

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
       call compList%addComp(c)
    end if
  end subroutine add_to_complist

  subroutine initialize_from_chain(cpar, handle, init_samp, init_from_output, first_call)
    implicit none
    type(comm_params), intent(in)           :: cpar
    type(planck_rng),  intent(inout)        :: handle
    integer(i4b),      intent(in), optional :: init_samp
    logical(lgt),      intent(in), optional :: init_from_output, first_call

    integer(i4b)              :: i, j, ext(2), initsamp, initsamp2
    character(len=4)          :: ctext
    character(len=6)          :: itext, itext2
    character(len=512)        :: chainfile, hdfpath
    class(comm_comp), pointer :: c => null()
    type(hdf_file) :: file, file2
    class(comm_N),      pointer :: N => null() 
    class(comm_map),    pointer :: rms => null()
    real(dp), allocatable, dimension(:,:) :: bp_delta, regnoise

    if (trim(cpar%init_chain_prefix) == 'none' .and. &
         & .not. present(init_from_output)) return
    
    ! Open default HDF file
    call int2string(cpar%mychain,   ctext)
    if (present(init_from_output)) then
       chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
            & '_c' // trim(adjustl(ctext)) // '.h5'
       initsamp = init_samp
    else
       call get_chainfile_and_samp(cpar%init_chain_prefix, chainfile, initsamp)
       if (present(init_samp)) initsamp = init_samp
    end if
    call int2string(initsamp, itext)
    call open_hdf_file(chainfile, file, 'r')
    
    !TODO: I get a crash here when the init file is missing or doesn't have the
    !required sample number, but it's some sort of MPI crash with no good
    !explanation or description so we should check for that somehow and write a
    !better error

    if (cpar%resamp_CMB .and. present(init_from_output)) then
       ! Initialize CMB component parameters; only once before starting Gibbs
       c   => compList
       do while (associated(c))
          if (.not. c%output) then
             c => c%nextComp()
             cycle
          end if
          call update_status(status, "init_chain_"//trim(c%label))
          if (cpar%myid == 0) write(*,*) '|  Initializing from chain = ', trim(c%label)
          call c%initHDFComp(cpar, file, trim(adjustl(itext))//'/')
          c => c%nextComp()
       end do
       call close_hdf_file(file)
       return
    end if
    
    ! Initialize instrumental parameters
    call update_status(status, "init_chain_inst")
    if (trim(cpar%cs_init_inst_hdf) /= 'none' .or. present(init_from_output)) then
       do i = 1, numband
          if (cpar%ignore_gain_bp) then
             data(i)%gain          = 1.d0
             data(i)%bp(0)%p%delta = 0.d0
          else
             if (trim(cpar%cs_init_inst_hdf) == 'default' .or. present(init_from_output)) then
                call read_hdf(file, trim(adjustl(itext))//'/gain/'//trim(adjustl(data(i)%label)), &
                  & data(i)%gain)
  
                call get_size_hdf(file, trim(adjustl(itext))//'/bandpass/'//&
                     & trim(adjustl(data(i)%label)), ext)
                if (data(i)%ndet > ext(1)-1) then
                   write(*,*) 'Error -- init HDF file ', trim(chainfile), ' does not contain enough bandpass information'
                   stop
                end if
                allocate(bp_delta(0:ext(1)-1,ext(2)))
                call read_hdf(file, trim(adjustl(itext))//'/bandpass/'//trim(adjustl(data(i)%label)), &
                     & bp_delta)
                do j = 0, data(i)%ndet
                   data(i)%bp(j)%p%delta = bp_delta(j,:)
                end do
                deallocate(bp_delta)
             else
                call get_chainfile_and_samp(cpar%cs_init_inst_hdf, &
                     & chainfile, initsamp2)
                call int2string(initsamp2, itext2)
                call open_hdf_file(chainfile, file2, 'r')
                call read_hdf(file2, trim(adjustl(itext2))//'/gain/'//trim(adjustl(data(i)%label)), &
                     & data(i)%gain)
  
                call get_size_hdf(file2, trim(adjustl(itext2))//'/bandpass/'//&
                     & trim(adjustl(data(i)%label)), ext)
                if (data(i)%ndet > ext(1)-1) then
                   write(*,*) 'Error -- init HDF file ', trim(chainfile), ' does not contain enough bandpass information'
                   stop
                end if
                allocate(bp_delta(0:ext(1)-1,ext(2)))
                call read_hdf(file2, trim(adjustl(itext2))//'/bandpass/'//trim(adjustl(data(i)%label)), &
                     & bp_delta)
                do j = 0, data(i)%ndet
                   data(i)%bp(j)%p%delta = bp_delta(j,:)
                end do
                deallocate(bp_delta)
                call close_hdf_file(file2)
             end if
          end if
       end do
    end if

    
    ! Initialize component parameters
    c   => compList
    do while (associated(c))
       if (.not. present(init_from_output) .and. &
            & (trim(c%init_from_HDF) == 'none' .or. &
            & (trim(c%type) == 'cmb' .and. .not. present(first_call)))) then
          c => c%nextComp()
          cycle
       end if
       call update_status(status, "init_chain_"//trim(c%label))
       if (cpar%myid == 0) write(*,*) '|  Initializing from chain = ', trim(c%label)
       if (trim(c%init_from_HDF) == 'default' .or. trim(c%init_from_HDF) == 'none') then
          call c%initHDFComp(cpar, file, trim(adjustl(itext))//'/')
       else
          call get_chainfile_and_samp(c%init_from_HDF, &
                     & chainfile, initsamp2)
          call int2string(initsamp2, itext2)
          call open_hdf_file(chainfile, file2, 'r')
          call c%initHDFComp(cpar, file2, trim(adjustl(itext2))//'/')
          call close_hdf_file(file2)
       end if
       c => c%nextComp()
    end do


    ! Initialize TOD parameters
    if (cpar%enable_TOD_analysis) then
       do i = 1, numband  
          if (trim(data(i)%tod_type) == 'none') cycle
          if (trim(data(i)%tod%init_from_HDF) == 'none' .and. .not. present(init_from_output))     cycle
          if (cpar%myid == 0) write(*,*) '|  Initializing TOD par from chain = ', trim(data(i)%tod%freq)
          N => data(i)%N
          rms => comm_map(data(i)%rmsinfo)
          select type (N)
          class is (comm_N_rms)
             if (trim(data(i)%tod%init_from_HDF) == 'default' .or. present(init_from_output)) then
                call data(i)%tod%initHDF(file, initsamp, data(i)%map, rms)
             else
                call get_chainfile_and_samp(data(i)%tod%init_from_HDF, &
                     & chainfile, initsamp2)
                call open_hdf_file(chainfile, file2, 'r')
                call data(i)%tod%initHDF(file2, initsamp2, data(i)%map, rms)
                call close_hdf_file(file2)
             end if
          class is (comm_N_rms_qucov)
             if (trim(data(i)%tod%init_from_HDF) == 'default' .or. present(init_from_output)) then
                call data(i)%tod%initHDF(file, initsamp, data(i)%map, rms)
             else
                call get_chainfile_and_samp(data(i)%tod%init_from_HDF, &
                     & chainfile, initsamp2)
                call open_hdf_file(chainfile, file2, 'r')
                call data(i)%tod%initHDF(file2, initsamp2, data(i)%map, rms)
                call close_hdf_file(file2)
             end if
          class default 
             write(*,*) 'Noise type is not covered'
          end select

          ! Update rms and data maps; add regularization noise if needed, no longer already included in the sample on disk
          allocate(regnoise(0:data(i)%info%np-1,data(i)%info%nmaps))
          if (associated(data(i)%procmask)) then
             call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
          else
             call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, map=rms)
          end if
          if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
          data(i)%map0%map = data(i)%map%map
          data(i)%map%map = data(i)%map%map + regnoise
          data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask
          call rms%dealloc
          deallocate(regnoise)
       end do

    else if (cpar%resamp_CMB) then
       do i = 1, numband  
          if (trim(data(i)%tod_type) == 'none') cycle
          !if (.not. data(i)%tod%init_from_HDF)  cycle
          if (cpar%myid == 0) write(*,*) '|  Initializing map and rms from chain = ', trim(data(i)%label), trim(data(i)%tod_type)

          hdfpath =  trim(adjustl(itext))//'/tod/'//trim(adjustl(data(i)%label))//'/'
          rms     => comm_map(data(i)%rmsinfo)
          N       => data(i)%N
          !write(*,*) trim(file%filename), trim(adjustl(hdfpath))//'map'
          call data(i)%map%readMapFromHDF(file, trim(adjustl(hdfpath))//'map')
          select type (N)
          class is (comm_N_rms)
             call rms%readMapFromHDF(file, trim(adjustl(hdfpath))//'rms')
          class is (comm_N_rms_qucov)
             call rms%readMapFromHDF(file, trim(adjustl(hdfpath))//'rms')
          class default
             write(*,*) 'resamp_CMB noise class not defined'
          end select

          ! Update rms and data maps; add regularization noise if needed, no longer already included in the sample on disk
          allocate(regnoise(0:data(i)%rmsinfo%np-1,data(i)%rmsinfo%nmaps))
          if (associated(data(i)%procmask)) then
             call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, procmask=data(i)%procmask, map=rms)
          else
             call data(i)%N%update_N(data(i)%rmsinfo, handle, data(i)%mask, regnoise, map=rms)
          end if
          if (cpar%only_pol) data(i)%map%map(:,1) = 0.d0
          data(i)%map%map = data(i)%map%map + regnoise
          data(i)%map%map = data(i)%map%map * data(i)%mask%map ! Apply mask
          call rms%dealloc
          deallocate(regnoise)
       end do
    end if


    ! Close HDF file
    call close_hdf_file(file)
    
  end subroutine initialize_from_chain


  subroutine sample_powspec(handle, ok)
    implicit none
    type(planck_rng),  intent(inout) :: handle
    logical(lgt),      intent(out)   :: ok

    class(comm_comp), pointer :: c => null()

    ok = .true.
    c  => compList
    do while (associated(c))
       if (trim(c%type) == 'md') then
          c => c%nextComp()
          cycle
       end if
       select type (c)
       class is (comm_diffuse_comp)
          call c%Cl%sampleCls(c%x, handle, ok)
       end select
       c => c%nextComp()
    end do

  end subroutine sample_powspec


  subroutine sample_partialsky_tempamps(cpar, handle)
    implicit none

    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: i, ierr
    logical(lgt) :: skip
    real(dp)     :: vals(2), vals2(2), mu, sigma, amp, mu_p, sigma_p
    class(comm_map),           pointer :: res => null(), invN_T => null()
    class(comm_comp),          pointer :: c => null()
    class(comm_template_comp), pointer :: pt => null()
    
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
          c => c%nextComp()
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

       c => c%nextComp()
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
    class(comm_comp), pointer :: c => null()
    class(comm_map),  pointer :: res => null()
    
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
       call res%dealloc(); deallocate(res)
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
             !if (c%x%info%myid ==0) write(*,*) bin, n
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
             !if (c%x%info%myid ==0) write(*,*) 'posdef', bin, c%Cl%bins2(bin)%lmin, posdef
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
                   if (chisq_prop < chisq_old) then
                      accept = .true.
                   else
                      accept = rand_uni(handle) < exp(-0.5d0*(chisq_prop-chisq_old))
                   end if
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
             call data(i)%c_old%dealloc(); deallocate(data(i)%c_old)
             call data(i)%c_prop%dealloc(); deallocate(data(i)%c_prop)
          end do

       end select
       c => c%nextComp()
    end do

  end subroutine sample_joint_alm_Cl


end module comm_signal_mod
