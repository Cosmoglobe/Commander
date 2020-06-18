module comm_nonlin_mod
  use comm_param_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_gain_mod
  use comm_line_comp_mod
  use comm_diffuse_comp_mod
  use comm_signal_mod
  implicit none

contains

!!$  subroutine sample_mono_dipole_with_mask(cpar, iter, handle)
!!$    implicit none
!!$    type(comm_params),  intent(in)    :: cpar
!!$    integer(i4b),       intent(in)    :: iter
!!$    type(planck_rng),   intent(inout) :: handle    
!!$
!!$    integer(i4b) :: i
!!$    class(comm_map),     pointer :: res
!!$    class(comm_comp),    pointer :: c
!!$    real(dp),          allocatable, dimension(:,:) :: m
!!$
!!$    ! Find monopole and dipole component
!!$    c => compList
!!$    do while (associated(c))
!!$       if (trim(c%label) /= 'md') then
!!$          c => c%next()
!!$          cycle
!!$       else
!!$          exit
!!$       end if
!!$    end do
!!$
!!$    ! Estimate monopoles and dipoles for each frequency
!!$    do i = 1, numband
!!$       ! Compute residual
!!$       res     => compute_residual(i)
!!$       m       = c%getBand(i)
!!$       res%map = res%map + m
!!$
!!$       call res%dealloc()
!!$       nullify(res)
!!$       deallocate(m)
!!$    end do
!!$    
!!$
!!$    nullify(c)
!!$
!!$
!!$    ! Sample spectral parameters for each signal component
!!$    allocate(status_fit(numband))
!!$    c => compList
!!$    do while (associated(c))
!!$       if (c%npar == 0) then
!!$          c => c%next()
!!$          cycle
!!$       end if
!!$       if (all(c%p_gauss(2,:) == 0.d0)) then
!!$          c => c%next()
!!$          cycle
!!$       end if
!!$       
!!$       do j = 1, c%npar
!!$
!!$          if (c%p_gauss(2,j) == 0.d0) cycle
!!$
!!$
!!$  end subroutine sample_mono_dipole_with_mask

  subroutine sample_nonlin_params(cpar, iter, handle, handle_noise)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle, handle_noise    

    integer(i4b) :: i, j, p
    real(dp)     :: t1, t2
    logical(lgt) :: samp_cg
    class(comm_comp),    pointer :: c    => null()

    call wall_time(t1)
                    
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%next()
          cycle
       end if
                    
       do j = 1, c%npar
          if (c%p_gauss(2,j) == 0.d0) cycle
          select type (c)
          class is (comm_diffuse_comp)
             !lmax_ind_pol is the lmax of poltype index p, for spec. ind. j 
             if (any(c%lmax_ind_pol(1:c%poltype(j),j) >= 0)) &
                  & call sample_specind_alm(cpar, iter, handle, c%id, j)
             if (any(c%lmax_ind_pol(1:c%poltype(j),j) < 0)) then
                call sample_specind_local(cpar, iter, handle, c%id, j)

                !check if any poltype has been sampled with ridge/marginal lnL
                samp_cg = .false.
                do p = 1,c%poltype(j)
                   if (p > c%nmaps) cycle
                   if (c%lmax_ind_pol(p,j) < 0 .and. &
                        & trim(c%pol_lnLtype(p,j)) /= 'chisq') samp_cg = .true.
                end do

                if (samp_cg) then !need to resample amplitude
                   !call sample amplitude for the component specific cg_sample group
                   if (cpar%myid == cpar%root) then
                      write(*,*) 'Sampling component amplitude of ',trim(c%label),' after spectral index sampling of ', &
                           & trim(c%indlabel(j))
                   end if
                   call sample_amps_by_CG(cpar, c%cg_unique_sampgroup, handle, handle_noise)
                end if
                !if/when 3x3 cov matrices are implemented, this CG-search needs to go inside local sampler routine (after every poltype index has been sampled)
             end if !any local sampling

          class is (comm_line_comp) !these codes should (maybe) not need to change
             call sample_specind_local(cpar, iter, handle, c%id, j)

          class is (comm_ptsrc_comp)
             call sample_specind_local(cpar, iter, handle, c%id, j)
          
          end select

       end do
       
       !go to next component
       c => c%next()
           
    end do

    ! Sample calibration factors
    do i = 1, numband
       if (.not. data(i)%sample_gain) cycle
       call sample_gain(cpar%operation, i, cpar%outdir, cpar%mychain, iter, mod(iter,cpar%resamp_hard_gain_prior_nth_iter)==0, handle)
    end do

    ! Update mixing matrices if gains have been sampled
    if (any(data%sample_gain)) then
       c => compList
       do while (associated(c))
          call c%updateMixmat
          c => c%next()
       end do
    end if

    call wall_time(t2)
    if (cpar%myid_chain == 0) write(*,*) 'CPU time specind = ', real(t2-t1,sp)
    
  end subroutine sample_nonlin_params


  subroutine sample_specind_alm(cpar, iter, handle, comp_id, par_id)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, r, q, p, pl, np, nlm, l_, m_, idx, delta, burnin, cholesky_calc
    integer(i4b) :: nsamp, out_every, check_every, num_accepted, smooth_scale, id_native, ierr, ind, nalm_tot_reg
    integer(i4b) :: p_min, p_max, nalm_tot, pix, region
    real(dp)     :: t1, t2, ts, dalm, thresh, steplen
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior, alms_mean, alms_var, chisq_jeffreys
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext
    character(len=3) :: tag
    character(len=9) :: ar_tag
    character(len=120) :: outmessage
    character(len=512) :: filename

    logical :: accepted, exist, doexit, optimize, apply_prior
    class(comm_mapinfo), pointer :: info => null()
    class(comm_mapinfo), pointer :: info_theta => null()
    class(comm_map),     pointer :: theta => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_map),     pointer :: theta_smooth => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_comp),    pointer :: c    => null()
    type(map_ptr),     allocatable, dimension(:) :: df

    real(dp),          allocatable, dimension(:,:,:)  :: alms, regs, buffer3
    real(dp),          allocatable, dimension(:,:)    :: m
    real(dp),          allocatable, dimension(:)      :: buffer, rgs, chisq, theta_pixreg_prop
    integer(c_int),    allocatable, dimension(:)      :: maxit


    ! Sample spectral parameter (parid) for the given signal component
    allocate(status_fit(numband))
    c => compList
    do while (c%id /= comp_id)
       c => c%next()
    end do

    select type (c)
    class is (comm_diffuse_comp)
       
       j = par_id !quick fix to only sample spec. ind. parameter par_id
             
       ! Set up smoothed data
       if (cpar%myid_chain == 0) write(*,*) '   Sampling '//trim(c%label)//' '//trim(c%indlabel(j))
       call update_status(status, "spec_alm start "//trim(c%label)//' '//trim(c%indlabel(j)))

       if (c%apply_jeffreys) then
          allocate(df(numband))
          do k = 1, numband
             df(k)%p => comm_map(data(k)%info)
          end do
       end if

       call wall_time(t1)

       info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
            & c%x%info%lmax, c%x%info%nmaps, c%x%info%pol)

       info_theta  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
            & 2*c%x%info%nside, 1, .false.)
       theta => comm_map(info_theta)

       ! Params
       out_every = 10
       check_every = 25 !100
       nsamp = cpar%almsamp_nsamp !2000
       burnin = cpar%almsamp_burnin ! Gibbs iter burnin. Tunes steplen.
       cholesky_calc = 1 ! Which gibbs iter to calculate cholesky, then corrlen.
       optimize = cpar%almsamp_optimize
       apply_prior = cpar%almsamp_apply_prior
       thresh = FLOAT(check_every)*0.8d0 !40.d0 ! 40.d0

       if (info%myid == 0 .and. c%L_read(j)) then
          write(*,*) "Sampling with cholesky matrix"
       end if

       if (info%myid == 0 .and. maxval(c%corrlen(j,:)) > 0) nsamp = maxval(c%corrlen(j,:))
       call mpi_bcast(nsamp, 1, MPI_INTEGER, 0, c%comm, ierr)

       ! Static variables
       doexit = .false.

       allocate(chisq(0:nsamp))
       chisq = 0.d0
       allocate(alms(0:nsamp, 0:c%nalm_tot-1,info%nmaps))             
       if (cpar%almsamp_pixreg) allocate(regs(0:nsamp, 0:MAXVAL(c%npixreg(:,j)),info%nmaps)) ! Region values            
       allocate(maxit(info%nmaps)) ! maximum iteration 
       maxit = 0

       ! Open log files
       if (info%myid == 0) open(69, file=trim(cpar%outdir)//'/nonlin-samples_'//trim(c%label)//'_'//trim(c%indlabel(j))//'.dat', status = 'unknown', access = 'append', recl=10000)
       if (info%myid == 0) open(66, file=trim(cpar%outdir)//'/region-samples_'//trim(c%label)//'_'//trim(c%indlabel(j))//'.dat', status = 'unknown', access = 'append', recl=10000)

       ! Save initial alm        
       alms = 0.d0
       regs = 0.d0
       ! Gather alms from threads to alms array with correct indices
       do pl = 1, c%theta(j)%p%info%nmaps
          call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, 0, pl, pl)
          allocate(buffer(c%nalm_tot))
          call mpi_allreduce(alms(0,:,pl), buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          alms(0,:,pl) = buffer
          if (cpar%almsamp_pixreg) regs(0,:,pl) = c%theta_pixreg(:,pl,j)
          deallocate(buffer)
       end do

       ! uniform fix
       !if (c%lmax_ind > 0 .and. alms(0,1:,:) == 0.d0) then
       !   alms(0,1:,:) = alms(0,1:,:) + 1e-6
       !   c%theta(j)%p%alm = c%theta(j)%p%alm + 1e-6
       !end if
       do pl = 1, c%theta(j)%p%info%nmaps
          ! if sample only pol, skip T
          if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 

          ! HKE -- disabling T for now
          if (pl==1) cycle 

          ! p already calculated if larger than poltype 
          if (pl > c%poltype(j)) cycle

          ! p to be sampled with a local sampler 
          if (c%lmax_ind_pol(pl,j) < 0) cycle

          if (cpar%almsamp_pixreg) then
             allocate(theta_pixreg_prop(0:c%npixreg(pl,j))) 
             allocate(rgs(0:c%npixreg(pl,j))) ! Allocate random vector
          else 
             allocate(rgs(0:c%nalm_tot-1)) ! Allocate random vector
          end if

          ! Get sampling tag
          if (c%poltype(j) == 1) then
             tag = "TQU"
          else if (c%poltype(j) == 2) then
             if (pl == 1) then
                tag = "T"
             else
                tag = "QU"
             end if
          else if (c%poltype(j) == 3) then
             if (pl == 1) then
                tag = "T"
             else if (pl == 2) then
                tag = "Q"
             else if (pl == 3) then
                tag = "U"
             end if
          end if

          !   if (c%theta(j)%p%info%nalm > 0) c%theta(j)%p%alm = (-4.d0 + 0.02d0*p)*sqrt(4.d0*pi)
          if (allocated(c%indmask)) then
             call compute_chisq(c%comm, chisq_fullsky=chisq(0), mask=c%indmask, lowres_eval=.true.)
          else
             call compute_chisq(c%comm, chisq_fullsky=chisq(0), lowres_eval=.true.)
          end if

          ! Use chisq from last iteration
          if (optimize .and. iter > 1 .and. chisq(0)>c%chisq_min(j,pl)) chisq(0) = c%chisq_min(j,pl)

          if (c%apply_jeffreys) then
             call c%updateMixmat(df=df, par=j)
             call compute_jeffreys_prior(c, df, pl, j, chisq_jeffreys)
             chisq(0) = chisq(0) + chisq_jeffreys
          end if

          call wall_time(t1)
          if (info%myid == 0) then 
             if (apply_prior) then
                if (.not. cpar%almsamp_pixreg) then
                   ! Add prior 
                   chisq_prior = ((alms(0,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                   if (c%nalm_tot > 1) then
                      do p = 1, c%nalm_tot-1
                         chisq_prior = chisq_prior + (alms(0,p,pl)/c%sigma_priors(p,j))**2
                      end do
                   end if
                else
                   ! Apply a prior per region
                   chisq_prior = 0.d0
                   do p = 1, c%npixreg(pl,j)
                      !write(*,*) "theta", c%theta_pixreg(p,pl,j), p, c%p_gauss(1,j)
                      chisq_prior = chisq_prior + (((c%theta_pixreg(p,pl,j) - c%p_gauss(1,j))/c%p_gauss(2,j))**2)
                   end do
                end if

                ! Output init sample
                write(*,fmt='(a, i6, a, f12.2, a, f6.2, a, 3f7.2)') "# sample: ", 0, " - chisq: " , chisq(0), " prior: ", chisq_prior,  " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
                if (cpar%almsamp_pixreg) write(*,fmt='(a,*(f7.3))') " regs:", c%theta_pixreg(1:,pl,j)
                chisq(0) = chisq(0) + chisq_prior
             else 
                write(*,fmt='(a, i6, a, f12.2, a, 3f7.2)') "# sample: ", 0, " - chisq: " , chisq(0),  " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
             end if

          end if

          ! Sample new alms (Account for poltype)
          num_accepted = 0

          nalm_tot = (c%lmax_ind_pol(pl,j) + 1)**2
          do i = 1, nsamp                   

             ! Gather alms from threads to alms array with correct indices
             call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)

             ! Send all alms to 0 (Dont allreduce because only root will do calculation)
             allocate(buffer(c%nalm_tot))
             call mpi_reduce(alms(i,:,pl), buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)
             alms(i,:,pl) = buffer
             deallocate(buffer)

             if (.not. cpar%almsamp_pixreg) then
                ! Propose new alms
                if (info%myid == 0) then
                   rgs = 0.d0
                   !cahnge nalm_tot
                   do p = 0, nalm_tot-1
                      rgs(p) = c%steplen(pl,j)*rand_gauss(handle)     
                   end do
                   alms(i,:,pl) = alms(i-1,:,pl) + matmul(c%L(:,:,pl,j), rgs)
                end if
             else
                ! --------- region sampling start
                !c%theta_pixreg(c%npixreg(pl,j),pl,j) = 0.d0 ! Just remove the last one for safe measure
                if (info%myid == 0) then
                   rgs = 0.d0
                   do p = 1, c%npixreg(pl,j)
                      rgs(p) = c%steplen(pl,j)*rand_gauss(handle)     
                   end do

                   ! Propose new pixel regions
                   theta_pixreg_prop = c%theta_pixreg(:,pl,j) + matmul(c%L(:c%npixreg(pl,j), :c%npixreg(pl,j), pl, j), rgs)  !0.05d0*rgs
                end if

                call mpi_bcast(theta_pixreg_prop, c%npixreg(pl,j)+1, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)

                ! Loop over pixels in region
                do pix = 0, theta%info%np-1
                   ! Else, use calculated change
                   theta%map(pix,1) = theta_pixreg_prop(c%ind_pixreg_arr(pix,pl,j))
                end do

                ! Smooth after regions are set, if smoothing scale > 0 and beam FWHM > 0.0
                smooth_scale = c%smooth_scale(j)
                if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
                   if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                      call smooth_map(info_theta, .false., &
                           & c%B_pp_fr(j)%p%b_l*0.d0+1.d0, theta, &  
                           & c%B_pp_fr(j)%p%b_l, theta_smooth)
                   else
                      theta_smooth => comm_map(info_theta)
                      theta_smooth%map=theta%map
                   end if
                else
                   theta_smooth => comm_map(info_theta)
                   theta_smooth%map=theta%map
                end if

                
                call theta_smooth%YtW_scalar
                call mpi_allreduce(theta_smooth%info%nalm, nalm_tot_reg, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
                allocate(buffer3(0:1, 0:nalm_tot_reg-1, 2)) ! Denne er nalm for 1! trenger alle

                call gather_alms(theta_smooth%alm, buffer3, theta_smooth%info%nalm, theta_smooth%info%lm, 0, 1, 1)
                call mpi_allreduce(MPI_IN_PLACE, buffer3, nalm_tot_reg, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                alms(i,:,pl) = buffer3(0,:c%nalm_tot,1)
                deallocate(buffer3)

                call theta_smooth%dealloc()
                ! ------- region sampling end
             end if

             ! Broadcast proposed alms from root
             allocate(buffer(c%nalm_tot))
             buffer = alms(i,:,pl)
             call mpi_bcast(buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)                   
             alms(i,:,pl) = buffer
             deallocate(buffer)

             ! Save to correct poltypes
             if (c%poltype(j) == 1) then      ! {T+E+B}
                do q = 1, c%theta(j)%p%info%nmaps
                   alms(i,:,q) = alms(i,:,pl) ! Save to all maps
                   call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)
                end do
             else if (c%poltype(j) == 2) then ! {T,E+B}
                if (pl == 1) then
                   call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                else
                   do q = 2, c%theta(j)%p%info%nmaps
                      alms(i,:,q) = alms(i,:,pl)
                      call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)                            
                   end do
                end if
             else if (c%poltype(j) == 3) then ! {T,E,B}
                call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
             end if

             ! Update mixing matrix with new alms
             if (c%apply_jeffreys) then
                call c%updateMixmat(df=df, par=j)
                call compute_jeffreys_prior(c, df, pl, j, chisq_jeffreys)
             else
                call c%updateMixmat
             end if

             ! Calculate proposed chisq
             if (allocated(c%indmask)) then
                call compute_chisq(c%comm, chisq_fullsky=chisq(i), mask=c%indmask, lowres_eval=.true.)
             else
                call compute_chisq(c%comm, chisq_fullsky=chisq(i), lowres_eval=.true.)
             end if

             ! Accept/reject test
             ! Reset accepted bool
             accepted = .false.
             if (info%myid == 0) then

                ! Adding prior
                if (c%apply_jeffreys) chisq(i) = chisq(i) + chisq_jeffreys

                if (apply_prior) then
                   if (.not. cpar%almsamp_pixreg) then
                      ! Apply priors per alm
                      chisq_prior = ((alms(i,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                      if (nalm_tot > 1) then
                         do p = 1, nalm_tot-1
                            !write(*,*) "alms ", p, alms(i,p,pl), c%sigma_priors(p,j)
                            chisq_prior = chisq_prior + (alms(i,p,pl)/c%sigma_priors(p,j))**2
                         end do
                      end if
                   else
                      ! Apply a prior per region
                      chisq_prior = 0.d0
                      do p = 1, c%npixreg(pl,j)
                         chisq_prior = chisq_prior + (((theta_pixreg_prop(p) - c%p_gauss(1,j))/c%p_gauss(2,j))**2)
                      end do
                      !write(*,*) "prior ", chisq_prior
                   end if

                   chisq(i) = chisq(i) + chisq_prior
                end if

                diff = chisq(i-1)-chisq(i)
                if ( chisq(i) > chisq(i-1)) then             
                   ! Small chance of accepting this too
                   ! Avoid getting stuck in local mminimum
                   if ( .not. (iter > burnin .and. trim(c%operation) == 'optimize')) then
                      accepted = (rand_uni(handle) < exp(0.5d0*diff))
                   end if
                else
                   accepted = .true.
                end if

                ! Count accepted and assign chisq values
                if (accepted) then
                   num_accepted = num_accepted + 1
                   ar_tag = achar(27)//'[92m'
                   if (cpar%almsamp_pixreg) c%theta_pixreg(:,pl,j) = theta_pixreg_prop
                else
                   chisq(i) = chisq(i-1)
                   ar_tag = achar(27)//'[91m'
                end if
                if (cpar%almsamp_pixreg) regs(i,:,pl) = c%theta_pixreg(:,pl,j)
                write(outmessage,fmt='(a, i6, a, f12.2, a, f8.2, a, f7.2, a, f7.4)') tag, i, " - chisq: " , chisq(i)-chisq_prior, " ", chisq_prior, " diff: ", diff, " - a00: ", alms(i,0,pl)/sqrt(4.d0*PI)
                write(*,*) adjustl(trim(ar_tag)//trim(outmessage)//trim(achar(27)//'[0m'))
                if (cpar%almsamp_pixreg) then
                   write(outmessage, fmt='(a, *(f7.3))') "regs:", theta_pixreg_prop(1:) ! Max space
                   write(*,*) adjustl(trim(ar_tag)//trim(outmessage)//trim(achar(27)//'[0m'))
                end if
             end if

             ! Broadcast result of accept/reject test
             call mpi_bcast(accepted, 1, MPI_LOGICAL, 0, c%comm, ierr)
             if (cpar%almsamp_pixreg) call mpi_bcast(c%theta_pixreg(:,pl,j), c%npixreg(pl,j)+1, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)

             if (.not. accepted) then
                ! If rejected, restore old values and send to 
                if (c%poltype(j) == 1) then      ! {T+E+B}
                   do q = 1, c%theta(j)%p%info%nmaps
                      alms(i,:,q) = alms(i-1,:,pl) ! Save to all mapsb
                      call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, q, q)
                   end do
                else if (c%poltype(j) == 2) then ! {T,E+B}
                   if (pl == 1) then
                      alms(i,:,pl) = alms(i-1,:,pl) 
                      call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, pl, 1)
                   else
                      do q = 2, c%theta(j)%p%info%nmaps
                         alms(i,:,q) = alms(i-1,:,pl)
                         call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, q, q)                            
                      end do
                   end if
                else if (c%poltype(j) == 3) then ! {T,E,B}
                   alms(i,:,pl) = alms(i-1,:,pl)
                   call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, pl, pl)
                end if

                call c%updateMixmat                                   
             end if


             if (info%myid == 0) then 
                ! Output log to file
                write(69, *) iter, tag, i, chisq(i), alms(i,:,pl)
                write(66, *) iter, tag, i, chisq(i), c%theta_pixreg(:, pl, j)

                ! Write to screen every out_every'th
                if (mod(i,out_every) == 0) then
                   diff = chisq(i-out_every) - chisq(i) ! Output diff
                   write(*,fmt='(a, i6, a, f12.2, a, f8.2, a, f7.2, a, f7.4)') " "//tag, i, " - chisq: " , chisq(i)-chisq_prior, " ", chisq_prior, " diff: ", diff, " - a00: ", alms(i,0,pl)/sqrt(4.d0*PI)
                   if (cpar%almsamp_pixreg) write(*,fmt='(a,*(f7.3))') " regs:", c%theta_pixreg(1:,pl,j)
                end if
                ! Adjust learning rate every check_every'th
                if (mod(i, check_every) == 0) then
                   diff = chisq(i-check_every)-chisq(i)                  

                   ! Accept rate
                   accept_rate = num_accepted/FLOAT(check_every)
                   num_accepted = 0

                   ! Write to screen
                   call wall_time(t2)
                   ts = (t2-t1)/DFLOAT(check_every) ! Average time per sample
                   write(*, fmt='(a, i6, a, i4, a, f8.2, a, f5.3, a, f5.2)') " "//tag, i, " - diff last ", check_every, " ", diff, " - accept rate: ", accept_rate, " - time/sample: ", ts
                   call wall_time(t1)

                   ! Adjust steplen in tuning iteration
                   if (iter <= burnin) then !( .not. c%L_read(j) .and. iter == 1) then ! Only adjust if tuning

                      if (accept_rate < 0.2) then                 
                         c%steplen(pl,j) = c%steplen(pl,j)*0.5d0
                         write(*,fmt='(a,f10.5)') "Reducing steplen -> ", c%steplen(pl,j)
                      else if (accept_rate > 0.45 .and. accept_rate < 0.55) then
                         c%steplen(pl,j) = c%steplen(pl,j)*2.0d0
                         write(*,fmt='(a,f10.5)') "Equilibrium - Increasing steplen -> ", c%steplen(pl,j)
                      else if (accept_rate > 0.8) then
                         c%steplen(pl,j) = c%steplen(pl,j)*2.d0
                         write(*,fmt='(a,f10.5)') "Increasing steplen -> ", c%steplen(pl,j)
                      end if
                   end if

                   ! Exit if threshold in tuning stage (First 2 iterations if not initialized on L)
                   if (c%corrlen(j,pl) == 0 .and. diff < thresh .and. accept_rate > 0.2 .and. i>=500) then
                      doexit = .true.
                      write(*,*) "Chisq threshold and accept rate reached for tuning iteration", thresh
                   end if
                end if
             end if

             call mpi_bcast(doexit, 1, MPI_LOGICAL, 0, c%comm, ierr)
             if (doexit .or. i == nsamp) then
                if (optimize) c%chisq_min(j,pl) = chisq(i) ! Stop increase in chisq
                if (info%myid == 0 .and. i == nsamp) write(*,*) "nsamp samples reached", nsamp

                ! Save max iteration for this signal
                if (c%poltype(j) == 1) then
                   maxit(:) = i 
                else if (c%poltype(j) == 2) then
                   if (pl==1) then
                      maxit(pl) = i
                   else
                      maxit(pl:) = i
                   end if
                else 
                   maxit(pl) = i
                end if
                doexit = .false. 
                exit
             end if

          end do ! End samples
          deallocate(rgs)
          if (cpar%almsamp_pixreg) deallocate(theta_pixreg_prop) 
       end do ! End pl



       if (info%myid == 0) close(58)

       ! Calculate correlation length and cholesky matrix 
       ! (Only if first iteration and not initialized from previous)
       if (info%myid == 0 .and. maxval(c%corrlen(j,:)) == 0) then
          if (c%L_read(j)  .and. iter >= burnin) then
             write(*,*) "Computing correlation function"

             do pl = 1, c%theta(j)%p%info%nmaps
                ! Skip signals with poltype tag
                if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                if (pl > c%poltype(j)) cycle
                if (c%lmax_ind_pol(pl,j) < 0) cycle

                if (cpar%almsamp_pixreg) then
                   call compute_corrlen(regs(:,1:,pl), c%npixreg(pl,j), maxit(pl), c%corrlen(j,pl))
                else
                   call compute_corrlen(alms(:,:,pl), nalm_tot, maxit(pl), c%corrlen(j,pl))
                end if

                write(*,*) "Correlation length (< 0.1): ", c%corrlen(j,pl) 
             end do

             ! If both corrlen and L have been calulated then output
             filename = trim(cpar%outdir)//'/init_alm_'//trim(c%label)//'_'//trim(c%indlabel(j))//'.dat'
             write(*,*) "Writing tuning parameters to file: ", trim(filename)
             open(58, file=filename, recl=10000)
             write(58,*) c%corrlen(j,:)
             write(58,*) c%L(:,:,:,j)
             close(58)

          else if (iter == cholesky_calc) then
             ! If L does not exist yet, calculate
             write(*,*) 'Calculating cholesky matrix'

             do pl = 1, c%theta(j)%p%info%nmaps
                if (maxit(pl) == 0) cycle ! Cycle if not sampled
                if (cpar%almsamp_pixreg) then
                   call compute_covariance_matrix(regs(INT(maxit(pl)/2):maxit(pl),:,pl), c%L(:c%npixreg(pl,j), :c%npixreg(pl,j), pl, j), .true.)
                else
                   call compute_covariance_matrix(alms(INT(maxit(pl)/2):maxit(pl),0:c%nalm_tot-1,pl), c%L(0:c%nalm_tot-1,0:c%nalm_tot-1,pl,j), .true.)
                end if
             end do
             c%steplen(:,j) = 1.d0
             c%L_read(j) = .true. ! L now exists!
          end if
       end if

       if (info%myid == 0) close(69)   
       if (info%myid == 0) close(66)   

       deallocate(alms, regs, chisq, maxit)
       call theta%dealloc()

       ! Clean up
       if (c%apply_jeffreys) then
          do k = 1, numband
             call df(k)%p%dealloc()
          end do
          deallocate(df)
       end if
    end select
    
    deallocate(status_fit)

  end subroutine sample_specind_alm

  subroutine sample_specind_local(cpar, iter, handle, comp_id, par_id)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx, p_ind, p_min, p_max
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind
    real(dp)     :: t1, t2, ts, dalm, fwhm_prior, temp_theta
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext
    logical :: accepted, exist, doexit, skip
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_comp),    pointer :: c    => null()
    class(comm_map),     pointer :: theta_single_pol => null() ! Spectral parameter of one poltype index (too be smoothed)
    real(dp),          allocatable, dimension(:,:,:)   :: alms
    real(dp),          allocatable, dimension(:,:) :: m
    real(dp),          allocatable, dimension(:) :: buffer, rgs, chisq

    integer(c_int),    allocatable, dimension(:,:) :: lm
    integer(i4b), dimension(MPI_STATUS_SIZE) :: mpistat

    call wall_time(t1)
    
    ! Initialize residual maps
    do i = 1, numband
       res             => compute_residual(i)
       data(i)%res%map =  res%map
       call res%dealloc()
       nullify(res)
    end do

    c           => compList     ! Extremely ugly hack...
    do while (c%id /= comp_id)
       c => c%next()
    end do

    ! Sample poltype index poltype_id of spec. ind. par_id of component comp_id
    allocate(status_fit(numband))

    ! Add current component back into residual
    if (trim(c%class) /= 'ptsrc') then
       do i = 1, numband
          allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
          m               = c%getBand(i)
          data(i)%res%map = data(i)%res%map + m
          deallocate(m)
       end do
    end if


    ! Set up smoothed data
    select type (c)
    class is (comm_line_comp)
       if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(par_id))

       ! Sample spectral parameters
       call c%sampleSpecInd(cpar, handle, par_id, iter)

    class is (comm_ptsrc_comp)
       if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label)

       ! Sample spectral parameters
       call c%sampleSpecInd(cpar, handle, par_id, iter)

    class is (comm_diffuse_comp)
       if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(par_id))
       call update_status(status, "nonlin start " // trim(c%label)// ' ' // trim(c%indlabel(par_id)))

       ! Set up type of smoothing scale
       id_native    = 0

       ! Compute smoothed residuals
       nullify(info)
       status_amp   = 0
       status_fit   = 0
       smooth_scale = c%smooth_scale(par_id)
       do i = 1, numband
          if (cpar%num_smooth_scales == 0) then
             status_fit(i)   = 1    ! Native
          else
             if (.not. associated(data(i)%N_smooth(smooth_scale)%p) .or. &
                  & data(i)%bp(0)%p%nu_c < c%nu_min_ind(par_id) .or. &
                  & data(i)%bp(0)%p%nu_c > c%nu_max_ind(par_id)) then
                status_fit(i) = 0
             else
                if (.not. associated(data(i)%B_smooth(smooth_scale)%p)) then
                   status_fit(i)   = 1 ! Native
                else
                   status_fit(i)   = 2 ! Smooth
                end if
             end if
          end if
                
          if (status_fit(i) == 0) then
             ! Channel is not included in fit
             nullify(res_smooth(i)%p)
             nullify(rms_smooth(i)%p)
          else if (status_fit(i) == 1) then
             ! Fit is done in native resolution
             id_native          = i
             info               => data(i)%res%info
             res_smooth(i)%p    => data(i)%res
             tmp                => data(i)%N
             select type (tmp)
                class is (comm_N_rms)
                rms_smooth(i)%p    => tmp
             end select
          else if (status_fit(i) == 2) then
             ! Fit is done with downgraded data
             info  => comm_mapinfo(data(i)%res%info%comm, cpar%nside_smooth(smooth_scale), cpar%lmax_smooth(smooth_scale), &
                  & data(i)%res%info%nmaps, data(i)%res%info%pol)
             call smooth_map(info, .false., data(i)%B(0)%p%b_l, data(i)%res, &
                  & data(i)%B_smooth(smooth_scale)%p%b_l, res_smooth(i)%p)
             rms_smooth(i)%p => data(i)%N_smooth(smooth_scale)%p
          end if

       end do
       status_amp = maxval(status_fit)

       ! Compute smoothed amplitude map
       if (.not. associated(info) .or. status_amp == 0) then
          write(*,*) 'Error: No bands contribute to index fit!'
          call mpi_finalize(i)
          stop
       end if
       if (status_amp == 1) then
          ! Smooth to the beam of the last native channel
          info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, c%x%info%lmax, &
               & c%x%info%nmaps, c%x%info%pol)
          call smooth_map(info, .true., data(id_native)%B(0)%p%b_l*0.d0+1.d0, c%x, &  
               & data(id_native)%B(0)%p%b_l, c%x_smooth)
       else if (status_amp == 2) then
          ! Smooth to the common FWHM
          info  => comm_mapinfo(c%x%info%comm, cpar%nside_smooth(smooth_scale), cpar%lmax_smooth(smooth_scale), &
               & c%x%info%nmaps, c%x%info%pol)
          call smooth_map(info, .true., &
               & data(1)%B_smooth(smooth_scale)%p%b_l*0.d0+1.d0, c%x, &  
               & data(1)%B_smooth(smooth_scale)%p%b_l,           c%x_smooth)
       end if

       ! Compute smoothed spectral index maps
       allocate(c%theta_smooth(c%npar))
       do k = 1, c%npar
          if (k == par_id) cycle
          if (status_amp == 1) then ! Native resolution
             info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
                  & c%x%info%lmax, c%x%info%nmaps, c%x%info%pol)
             call smooth_map(info, .false., &
                  & data(id_native)%B(0)%p%b_l*0.d0+1.d0, c%theta(k)%p, &  
                  & data(id_native)%B(0)%p%b_l,           c%theta_smooth(k)%p)
          else if (status_amp == 2) then ! Common FWHM resolution
             info  => comm_mapinfo(c%theta(k)%p%info%comm, cpar%nside_smooth(smooth_scale), &
                  & cpar%lmax_smooth(smooth_scale), c%theta(k)%p%info%nmaps, c%theta(k)%p%info%pol)
             call smooth_map(info, .false., &
                  & data(1)%B_smooth(smooth_scale)%p%b_l*0.d0+1.d0, c%theta(k)%p, &  
                  & data(1)%B_smooth(smooth_scale)%p%b_l,           c%theta_smooth(k)%p)
          end if
       end do

       ! Sample spectral parameters for diffuse component
       call sampleDiffuseSpecInd_nonlin(cpar, handle, c%id, par_id, iter)

    end select
          
          
    ! Clean up temporary data structures
    select type (c)
    class is (comm_line_comp)
    class is (comm_diffuse_comp)
       
       if (associated(c%x_smooth)) then
          call c%x_smooth%dealloc()
          nullify(c%x_smooth)
       end if
       do k =1, c%npar
          if (k == par_id) cycle
          if (allocated(c%theta_smooth)) then
             if (associated(c%theta_smooth(k)%p)) then
                call c%theta_smooth(k)%p%dealloc()
             end if
          end if
       end do
       if (allocated(c%theta_smooth)) deallocate(c%theta_smooth)
       do i = 1, numband
          if (.not. associated(rms_smooth(i)%p)) cycle
          if (status_fit(i) == 2) then
             call res_smooth(i)%p%dealloc()
          end if
          nullify(res_smooth(i)%p)
       end do

       smooth_scale = c%smooth_scale(par_id)
       if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
          if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
             ! Smooth index map with a postprocessing beam
             !deallocate(c%theta_smooth)             
             !We rewrite this to smooth the parameter to the same nside as the component, 
             !so that the evaluation/sampling can be done at an arbitrary nside (given by smoothing scale),
             !then be smoothed with a post processing beam at the components full resolution
             info  => comm_mapinfo(c%theta(par_id)%p%info%comm, c%theta(par_id)%p%info%nside, &
                  & 2*c%theta(par_id)%p%info%nside, 1, c%theta(par_id)%p%info%pol) !only want 1 map

             !spec. ind. map with 1 map (will be smoothed like zero spin map using the existing code)
             theta_single_pol => comm_map(info)
             
             do p_ind = 1,c%poltype(par_id)
                if (c%lmax_ind_pol(p_ind,par_id) >= 0) cycle !alm sampled
                if (c%poltype(par_id) == 1) then
                   p_min = 1; p_max = c%nmaps
                   if (cpar%only_pol) p_min = 2
                else if (c%poltype(par_id) == 2) then
                   if (p_ind == 1) then
                      p_min = 1; p_max = 1
                   else
                      p_min = 2; p_max = c%nmaps
                   end if
                else if (c%poltype(par_id) == 3) then
                   p_min = p_ind
                   p_max = p_ind
                end if

                if (p_min<=p_max) then !just a guaranty that we dont smooth for nothing
                   !copy poltype id theta map to the single map
                   theta_single_pol%map(:,1)=c%theta(par_id)%p%map(:,p_min)
                   allocate(c%theta_smooth(c%npar))
                   !smooth single map as intensity (i.e. zero spin)
                   call smooth_map(info, .false., &
                        & c%B_pp_fr(par_id)%p%b_l*0.d0+1.d0, theta_single_pol, &  
                        & c%B_pp_fr(par_id)%p%b_l, c%theta_smooth(par_id)%p)

                   do p = p_min,p_max
                      ! assign smoothed theta map to relevant polarizations
                      c%theta(par_id)%p%map(:,p) = c%theta_smooth(par_id)%p%map(:,1)
                   end do
                   call c%theta_smooth(par_id)%p%dealloc()
                   deallocate(c%theta_smooth)
                end if
             end do
             call theta_single_pol%dealloc()
             theta_single_pol => null()

          end if
       end if

       call update_status(status, "nonlin stop " // trim(c%label)// ' ' // trim(c%indlabel(par_id)))

    end select

    call c%updateMixmat 


    !This is a nice cleanup, but if new residuals are drawn for each time sample_specind_local is called,
    !then this is a step that is not necessary to perform. Leaving it in for now. 
    ! Subtract updated component from residual
    if (trim(c%class) /= 'ptsrc') then
       do i = 1, numband
          allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
          m               = c%getBand(i)
          data(i)%res%map = data(i)%res%map - m
          deallocate(m)
       end do
    end if
   
  end subroutine sample_specind_local
  

  subroutine gather_alms(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(in)    :: alm
    integer(c_int), dimension(1:,0:), intent(in) :: lm
    real(dp), dimension(0:,0:,1:), intent(inout) :: alms
    integer(i4b),                intent(in)    :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind

    do k = 0, nalm-1
       ! Gather all alms
       l = lm(1,k)
       m = lm(2,k)
       ind = l**2 + l + m
       alms(i,ind,pl_tar) = alm(k,pl)
    end do

  end subroutine gather_alms

  subroutine distribute_alms(alm, alms, nalm, lm, i, pl, pl_tar)
    implicit none

    real(dp), dimension(0:,1:),    intent(inout)    :: alm
    integer(c_int), dimension(1:,0:), intent(in)   :: lm
    real(dp), dimension(0:,0:,1:),  intent(in)       :: alms
    integer(i4b),                intent(in)       :: nalm, i, pl, pl_tar
    integer(i4b) :: k, l, m, ind
    
    do k = 0, nalm-1
       ! Distribute alms
       l = lm(1,k)
       m = lm(2,k)
       ind = l**2 + l + m
       alm(k,pl_tar) = alms(i,ind,pl)
    end do

  end subroutine distribute_alms

  subroutine compute_corrlen(x, n, maxit, corrlen)
    implicit none

    real(dp), dimension(:,:),    intent(in)    :: x
    integer(i4b),                  intent(in)    :: n
    integer(i4b),                  intent(in)    :: maxit
    integer(i4b),                  intent(out)   :: corrlen

    real(dp),          allocatable, dimension(:) :: C_
    integer(c_int),    allocatable, dimension(:) :: N_      
    real(dp)     :: x_mean, x_var
    integer(i4b) :: pl, p, q, k, corrlen_init, delta

    ! Calculate Correlation length
    delta = 100
    corrlen_init = 1
    allocate(C_(delta))
    allocate(N_(delta))
    
    !open(58, file='correlation_function.dat', recl=10000)

    corrlen = corrlen_init
           
    ! Calculate correlation function per parameter
    do p = 1, n
       x_mean = mean(x(1:maxit,p))
       x_var = variance(x(1:maxit,p))
     
       
       N_ = 0
       C_ = 0.d0
       do q = 1, maxit
          do k = 1, delta
             if (q+k > maxit) cycle
             C_(k) = C_(k) + (x(q,p)-x_mean)*(x(q+k,p)-x_mean)
             N_(k) = N_(k) + 1 ! Less samples every q
          end do
       end do
       
       where (N_>0) C_ = C_/N_
       if ( x_var > 0 ) C_ = C_/x_var

      ! write(58,*) p, C_ ! Write to file

       ! Find correlation length
       do k = corrlen_init, delta
          if (C_(k) > 0.1) then
             if (k > corrlen) corrlen = k
          else
             exit
          end if
       end do
    end do
    deallocate(C_, N_)
    close(58)
  end subroutine compute_corrlen


  !Here comes all subroutines for sampling diffuse components locally
  ! Sample spectral parameters
  subroutine sampleDiffuseSpecInd_nonlin(cpar, handle, comp_id, par_id, iter)
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: par_id, comp_id, iter

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer_lnL, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, accept_rate
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample
    class(comm_mapinfo),            pointer :: info => null()

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    first_call          = .false.


    id = par_id
    c           => compList     ! Extremely ugly hack...
    do while (comp_id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
       info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
            & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)


       npar      = c_lnL%npar
       np        = c_lnL%x_smooth%info%np
       nmaps     = c_lnL%x_smooth%info%nmaps

       theta_min = c_lnL%p_uni(1,id)
       theta_max = c_lnL%p_uni(2,id)

       !needs rewriting to only sample polarizations covered by polt_id
       if (trim(cpar%operation) == 'optimize') then
          call c_lnL%sampleSpecInd(cpar, handle, par_id, iter) !use code in diffuse_comp_mod

       else if (trim(cpar%operation) == 'sample') then

          allocate(buffer_lnL(0:c_lnL%theta(id)%p%info%np-1,c_lnL%theta(id)%p%info%nmaps))
          buffer_lnL=max(min(c_lnL%theta(id)%p%map,theta_max),theta_min) 
          

          do p = 1,c_lnL%poltype(id)
             if (c_lnL%lmax_ind_pol(p,id) >= 0) cycle !this set of polarizations are not to be local sampled (is checked before this point)
             if (c_lnL%poltype(id) > 1 .and. cpar%only_pol .and. p == 1) cycle !only polarization (poltype > 1)
             if (p > c_lnL%nmaps) cycle ! poltype > number of maps



             call wall_time(t1)
             if (c_lnL%pol_pixreg_type(p,id) > 0) then
                if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) 'Sampling poltype index', p, &
                     & 'of ', c_lnL%poltype(id) !Needed?
                call sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
                call wall_time(t2)
                if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) 'poltype:',c_lnL%poltype(id),' pol:', &
                     & p,'CPU time specind = ', real(t2-t1,sp)
             else
                write(*,*) 'Undefined spectral index sample region'
                write(*,*) 'Component:',trim(c_lnL%label),'ind:',trim(c_lnL%indlabel(id))
                stop
             end if



          end do

          !after sampling is done we assign the spectral index its new value(s)
          c_lnL%theta(id)%p%map = buffer_lnL
          
          if (info%myid == 0 .and. cpar%verbosity > 2) write(*,*) 'Updating Mixing matrix'
          ! Update mixing matrix
          call c_lnL%updateMixmat

          ! Ask for CG preconditioner update
          if (c_lnL%cg_unique_sampgroup > 0) recompute_diffuse_precond = .true.

          ! deallocate

          deallocate(buffer_lnL)

       end if !operation

    end select

  end subroutine sampleDiffuseSpecInd_nonlin

  subroutine sampleDiffuseSpecIndSinglePix_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    real(dp),               dimension(0:,:), intent(inout)        :: buffer_lnL
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, par_id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL
    real(dp)     :: accept_rate, accept_scale, lnL_sum, sum_proplen, sum_accept
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, burn_in, count_accept
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample, use_det, burned_in
    real(dp),     allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),     allocatable, dimension(:) :: theta_corr_arr
    integer(i4b), allocatable, dimension(:) :: band_i, pol_j
    class(comm_mapinfo),            pointer :: info => null()

    c           => compList     ! Extremely ugly hack...
    do while (comp_id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select
    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)
                
                
    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.

    id = par_id  !hack to not rewrite the entire code from diffuse_comp_mod

    n_prop_limit = 40
    n_corr_limit = 40
    count_accept = 0
    sum_accept   = 0.d0
    sum_proplen   = 0.d0

    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (cpar%only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np        = c_lnL%x_smooth%info%np
    nmaps     = c_lnL%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)

    !set up which bands and polarizations to include
    !this is to exclude the need of checking this for each time the mixing matrix needs to be computed per pixel
    allocate(band_i(1000),pol_j(1000))
    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the component, i.e. frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < c_lnL%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > c_lnL%nu_max_ind(id)) cycle

       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set all thata to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (info%myid == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (info%myid == 0 .and. .true.) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
    end if

    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))
    do pix = 0,np-1 !loop over the number of pixels the processor is handling
       
       if (c_lnL%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) then     ! if pixel is masked out
          buffer_lnL(pix,p_min:p_max) = c_lnL%p_gauss(1,id) ! set spec. ind to prior
          c_lnL%pol_nprop(id)%p%map(pix,p)=0.d0 !to mark that the pixel is masked out, we set nprop to zero
          cycle
       end if
       
       if (c_lnL%pol_sample_nprop(p,id) .or. c_lnL%pol_sample_proplen(p,id)) then
          pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
          allocate(theta_corr_arr(pixreg_nprop))
          n_spec_prop = 0
          n_corr_prop = 0 
          !c_lnL%pol_sample_nprop(p,id) = boolean array of size (poltype,npar)
          !c_lnL%pol_sample_proplen(p,id) = boolean array of size (poltype,npar)
       else
          pixreg_nprop = IDINT(c_lnL%pol_nprop(id)%p%map(pix,p)) ! comm_map is real(dp), use IDINT to get the integer value
       end if

       n_accept = 0
       first_sample = .true.

       do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)
          
          if (first_sample) then
             old_theta = buffer_lnL(pix,p)

             !get the values of the remaining spec inds of the component
             do i = 1, npar
                if (i == id) cycle
                all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p)
             end do
          end if

          !draw new spec ind
          new_theta = old_theta + c_lnL%pol_proplen(id)%p%map(pix,p)*rand_gauss(handle)

          !init log-like for new sample
          if (first_sample) lnL_old = 0.d0
          lnL_new = 0.d0
          
          if (new_theta > theta_max .or. new_theta < theta_min) then !new proposal is outside limits
             lnL_new = -1.d30 
             ! skip the true caclulation of lnL, we reject the sample ~100%
          else
             ! do a clever way of computing the lnL for the different evaluation types
             if (trim(c_lnL%pol_lnLtype(p,id))=='chisq') then
                !write the chisq sampling here

                do k = 1,band_count !run over all active bands
                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind.
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                   end if
                   all_thetas(id)=new_theta
                   mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   !compute chisq for pixel 'pix' of the associated band
                   if (first_sample) then
                      res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old*c_lnL%x_smooth%map(pix,pol_j(k))
                      lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                   end if
                   res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new*c_lnL%x_smooth%map(pix,pol_j(k))
                   lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                end do

             else if ((trim(c_lnL%pol_lnLtype(p,id))=='ridge') .or. &
                  & (trim(c_lnL%pol_lnLtype(p,id))=='marginal')) then
                if (trim(c_lnL%pol_lnLtype(p,id))=='ridge') then
                   use_det = .false.
                else
                   use_det = .true.
                end if

                !! build mixing matrix, invN and data for pixel
                do k = 1,band_count !run over all active bands
                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind.
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                   end if
                   all_thetas(id)=new_theta
                   mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                   invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                end do

                !compute the marginal/ridge log-likelihood for the pixel
                if (first_sample) then
                   lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end if
                lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                     & use_det, band_count)

             else 
                write(*,*) 'invalid polarized lnL sampler type'
                stop

             end if !chisq/marginal/ridge

             !add spec. ind. priors
             if (first_sample) then
                all_thetas(id)=old_theta
                do l = 1, npar
                   if (c_lnL%p_gauss(2,l) > 0.d0) then
                      lnL_old = lnL_old - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                   end if
                   all_thetas(id)=new_theta
                end do
             end if
             do l = 1, npar
                if (c_lnL%p_gauss(2,l) > 0.d0) then
                   lnL_new = lnL_new - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                end if
             end do
             !first sample done
             first_sample = .false.

          end if !outside spec limits

          !accept/reject new spec ind
          delta_lnL = lnL_new-lnL_old

          if (delta_lnL < 0.d0) then
             !if delta_lnL is more negative than -25.d0, limit to -25.d0
             if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

             !draw random uniform number
             a = rand_uni(handle) !draw uniform number from 0 to 1
             if (exp(delta_lnL) > a) then
                !accept
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
             end if
          else
             !accept new sample, higher likelihood
             old_theta = new_theta
             lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
             n_accept = n_accept + 1
          end if
          
          if (j > burn_in) burned_in = .true.
          ! evaluate proposal_length, then correlation length. 
          !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
          if (c_lnL%pol_sample_proplen(p,id)) then
             n_spec_prop = n_spec_prop + 1
             accept_rate = n_accept*1.d0/n_spec_prop
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_spec_prop = 0
                n_accept = 0
             else if (n_spec_prop > n_prop_limit) then
                if (accept_rate > 0.7d0) then
                   accept_scale = 1.5d0
                else if (accept_rate < 0.3d0) then
                   accept_scale = 0.5d0
                else 
                   accept_scale = -1.d0
                end if

                if (accept_scale > 0.d0) then
                   c_lnL%pol_proplen(id)%p%map(pix,p) = c_lnL%pol_proplen(id)%p%map(pix,p)* &
                        & accept_scale !set new proposal length
                   n_accept = 0 !reset with new prop length
                   n_spec_prop = 0 !reset with new prop length
                else
                   sum_accept = sum_accept + accept_rate
                   count_accept = count_accept + 1
                   sum_proplen = sum_proplen + c_lnL%pol_proplen(id)%p%map(pix,p)
                   if (.not. c_lnL%pol_sample_nprop(p,id)) then
                      exit 
                      !ugly, but we only end up here in the first gibbs sample if we are only fitting proposal
                      !lengths and not correlation length
                   end if
                end if
             end if

          else if (c_lnL%pol_sample_nprop(p,id)) then
             ! if we are keeping track of spec inds for corr. length, then save current spec ind
             n_corr_prop = n_corr_prop + 1
             theta_corr_arr(n_corr_prop) = old_theta
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_corr_prop = 0
             else if (n_corr_prop > n_corr_limit) then
                corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                if (corr_len > 0) then ! .or. n_corr_prop > 4*c_lnL%nprop_uni(2,id)) then
                   !set nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                   c_lnL%pol_nprop(id)%p%map(pix,p) = 1.d0*min(max(2*corr_len,c_lnL%nprop_uni(1,id)),c_lnL%nprop_uni(2,id))
                   exit 
                   !ugly, but we only end up here in the first gibbs sample if we are fitting correlation length

                end if
             end if
          else
             accept_rate = n_accept*1.d0/j
             sum_accept = sum_accept + accept_rate
             count_accept = count_accept + 1
             sum_proplen = sum_proplen + c_lnL%pol_proplen(id)%p%map(pix,p)
          end if
          
       end do !pixreg_nprop


       if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)
       
       !assign last accepted theta value (old_theta) to the pixels of the current pixel region, 
       !in the polarizations given by poltype and p
       buffer_lnL(pix,p_min:p_max) = old_theta
       
    end do !pix = 1,np


    call mpi_allreduce(MPI_IN_PLACE, sum_accept, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, count_accept, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, sum_proplen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    if (info%myid==0) then
       write(*,fmt='(a,f6.4)') '   avg. accept rate     = ', sum_accept/count_accept
       write(*,fmt='(a,e14.4)') '   avg. proposal length = ', sum_proplen/count_accept
    end if
    if (iter >= 2) then !only stop sampling after 2nd iteration, in case first iter is far off.
       c_lnL%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals and
       c_lnL%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done first gibbs iteration)
    end if
    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)


  end subroutine sampleDiffuseSpecIndSinglePix_nonlin


  subroutine sampleDiffuseSpecIndFullsky_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    real(dp),               dimension(0:,:), intent(inout)        :: buffer_lnL
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, par_id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: buff1_r(1), buff2_r(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL
    real(dp)     :: accept_rate, accept_scale, lnL_sum, init_theta
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, pix_count, buff1_i(1), buff2_i(1), burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in
    real(dp),     allocatable, dimension(:)   :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),     allocatable, dimension(:)   :: theta_corr_arr
    integer(i4b), allocatable, dimension(:)   :: band_i, pol_j
    class(comm_mapinfo),            pointer   :: info => null()

    c           => compList     ! Extremely ugly hack...
    do while (comp_id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select
    info  => comm_mapinfo(c_lnL%x%info%comm, c_lnL%x%info%nside, &
         & c_lnL%x%info%lmax, c_lnL%x%info%nmaps, c_lnL%x%info%pol)
                

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.


    id = par_id !hack to not rewrite too much from diffuse_comp_mod
    
    n_prop_limit = 20
    n_corr_limit = 20
    out_every    = 10

    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (cpar%only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np        = c_lnL%x_smooth%info%np
    nmaps     = c_lnL%x_smooth%info%nmaps

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)

    !set up which bands and polarizations to include

    allocate(band_i(3*numband),pol_j(3*numband))

    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the smoothed component, and if band frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < c_lnL%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > c_lnL%nu_max_ind(id)) cycle
                   
       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set theta to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (info%myid == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (info%myid == 0 .and. .true.) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
    end if

    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))

    n_spec_prop = 0
    n_accept = 0
    n_corr_prop = 0 
    if (c_lnL%pol_sample_nprop(p,id) .or. c_lnL%pol_sample_proplen(p,id)) then
       pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
       allocate(theta_corr_arr(pixreg_nprop))
       !c_lnL%pol_sample_nprop(j,p,id) = boolean array of size (n_pixreg,poltype)
       !c_lnL%pol_sample_proplen(j,p,id) = boolean array of size (n_pixreg,poltype)
    else
       pixreg_nprop = IDINT(c_lnL%pol_nprop(id)%p%map(0,p)) ! comm_map is real(dp), use IDINT() to convert to integer
    end if

    first_sample  = .true.
    loop_exit     = .false.
    
    call wall_time(t1)
    do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)
       if (info%myid == 0) then
          if (first_sample) then
             !taking the avg fullsky value
             old_theta = sum(buffer_lnL(:,p_min)) 
             pix_count = np
             call mpi_allreduce(MPI_IN_PLACE, old_theta, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, pix_count, 1, MPI_INTEGER, & 
                  & MPI_SUM, info%comm, ierr)
             old_theta = old_theta/pix_count

             init_theta = old_theta
             write(*,fmt='(a, f10.5)') "  initial sepc. ind. value: ", init_theta
          end if

          !draw new spec ind, 
          new_theta = old_theta + c_lnL%pol_proplen(id)%p%map(0,p)*rand_gauss(handle)

       end if


       !broadcast new_theta, and new old_theta, and all_thetas to the other processors
       if (first_sample) call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, &
            & 0, c_lnL%comm, ierr)
       call mpi_bcast(new_theta, 1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)

       if (info%myid==0 .and. mod(j,out_every)==0) then
          write(*,fmt='(a, i6, a, f10.5, a, f10.5)') "  proposal: ", j," -- Current ind: ", &
               & old_theta, " -- proposed ind: ", new_theta
       end if


       !init log-like for new sample
       if (first_sample) lnL_old = 0.d0
       lnL_new = 0.d0       

       if (new_theta > theta_max .or. new_theta < theta_min) then !new proposal is outside limits
          lnL_new = -1.d30 
          ! skip the true caclulation of lnL, we reject the sample ~100%
          if (info%myid==0) write(*,fmt='(a, f10.5, a, f10.5)') "    Proposed ind outside limits.  min: ", &
               & theta_min," -- max: ", theta_max
       else
          pix_count = 0
          all_thetas(id)=new_theta
          !lnL type should split here
          if (trim(c_lnL%pol_lnLtype(p,id))=='chisq') then
             !loop bands and then pixels

             do pix = 0,np-1 !loop over pixels covered by the processor
                if (c_lnL%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                end do

                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                   ! both for old and new spec. ind. and calc. log likelihood (chisq)
                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta

                      res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old* &
                           & c_lnL%x_smooth%map(pix,pol_j(k))
                      lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2

                   end if
                   mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale

                   res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                        & c_lnL%x_smooth%map(pix,pol_j(k))
                   lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                end do
             end do

          else if ((trim(c_lnL%pol_lnLtype(p,id))=='ridge') .or. &
               & (trim(c_lnL%pol_lnLtype(p,id))=='marginal')) then
             if (trim(c_lnL%pol_lnLtype(p,id))=='ridge') then
                use_det = .false.
             else
                use_det = .true.
             end if

             do pix = 0,np-1 !More precise, we need to loop over pixels covered by the processor
                
                if (c_lnL%pol_ind_mask(id)%p%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                !get the values of the remaining spec inds of the component for the given pixel
                do i = 1, npar
                   if (i == id) cycle
                   all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                end do

                !build mixing matrix
                do k = 1,band_count !run over all active bands
                   pix_count = pix_count + 1

                   if (first_sample) then
                      all_thetas(id)=old_theta
                      mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      all_thetas(id)=new_theta
                   end if
                   mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                        & data(band_i(k))%gain * c_lnL%cg_scale
                   pix_count = pix_count + 1
                   data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                   invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                end do

                !compute the marginal/ridge log-like for the pixel
                if (first_sample) then
                   lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end if
                lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                     & use_det, band_count)
             end do
          else 
             if (info%myid==0) write(*,*) 'invalid polarized lnL sampler type'
             stop

          end if !chisq/marginal/ridge


          !Sum lnL_new (and lnL_old) among all processors
          if (first_sample) then
             !there is a problem with the mpi_allreduce call. Commander3 crashes
             call mpi_allreduce(lnL_old, buff2_r, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info%comm, ierr)
             lnL_old = buff2_r(1)

             call mpi_allreduce(pix_count, k, 1, MPI_INTEGER, & 
                  & MPI_SUM, info%comm, ierr)
             pix_count = k

             if (pix_count == 0) then! there are no data that is valid
                old_theta = c_lnL%p_gauss(1,id) !set theta to prior
                !set nprop to 0, write error message, then exit sampling
                c_lnL%pol_nprop(id)%p%map(:,p) = 0.d0
                if (info%myid==0) write(*,*) "    No valid pixels to sample spectral index from"
                exit
             end if

          end if

          call mpi_allreduce(lnL_new, lnl_sum, 1, MPI_DOUBLE_PRECISION, & 
               & MPI_SUM, info%comm, ierr)
          lnL_new = lnL_sum

          !add spec. ind. priors
          if (info%myid == 0) then
             if (first_sample) then
                all_thetas(id)=old_theta
                do l = 1, npar
                   if (c_lnL%p_gauss(2,l) > 0.d0) then
                      lnL_old = lnL_old - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                   end if
                   all_thetas(id)=new_theta
                end do
             end if
             do l = 1, npar
                if (c_lnL%p_gauss(2,l) > 0.d0) then
                   lnL_new = lnL_new - 0.5d0 * (all_thetas(l)-c_lnL%p_gauss(1,l))**2 / c_lnL%p_gauss(2,l)**2 
                end if
             end do
          end if

          !first sample done
          first_sample = .false.

          if (pix_count==0) then
             !set theta to prior
             old_theta = c_lnL%p_gauss(1,id)
             !set number of proposals to zero, i.e. spectral index will not be sampled in later iterations 
             !as there are no unmasked data
             c_lnL%pol_nprop(id)%p%map(:,p)=0.d0
             exit   !exit sampling
          end if

       end if !new_theta outside spec limits

       if (info%myid == 0) then
          !accept/reject new spec ind
          delta_lnL = lnL_new-lnL_old

          if (mod(j,out_every)==0) write(*,fmt='(a, e14.5)') "    lnL_new - lnL_old = ", delta_lnL


          if (delta_lnL < 0.d0) then
             !if delta_lnL is more negative than -25.d0, limit to -25.d0
             if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

             !draw random uniform number
             a = rand_uni(handle) !draw uniform number from 0 to 1
             if (exp(delta_lnL) > a) then
                !accept
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
             end if
          else
             !accept new sample, higher likelihood
             old_theta = new_theta
             lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
             n_accept = n_accept + 1
          end if

          if (j > burn_in) burned_in = .true.
          ! evaluate proposal_length, then correlation length. 
          !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
          if (c_lnL%pol_sample_proplen(p,id)) then
             n_spec_prop = n_spec_prop + 1
             accept_rate = n_accept*1.d0/n_spec_prop
             if (mod(n_spec_prop,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
             if (.not. burned_in) then 
                !giving som time to let the first steps "burn in" if one starts far away from max lnL
                n_spec_prop = 0
                n_accept = 0
             else if (n_spec_prop > n_prop_limit) then
                if (accept_rate > 0.7d0) then
                   accept_scale = 1.5d0
                else if (accept_rate < 0.3d0) then
                   accept_scale = 0.5d0
                else 
                   accept_scale = -1.d0
                end if

                if (accept_scale > 0.d0) then
                   !prop_len only ever used by root processor
                   c_lnL%pol_proplen(id)%p%map(0,p) = c_lnL%pol_proplen(id)%p%map(0,p)*accept_scale !set new proposal length
                   n_accept = 0 !reset with new prop length
                   n_spec_prop = 0 !reset with new prop length
                   write(*,fmt='(a, e14.5)') "     New prop. len. = ", c_lnL%pol_proplen(id)%p%map(0,p)
                else
                   c_lnL%pol_sample_proplen(p,id)=.false. !making sure we dont go back into prop. len. sampling
                   !only necessary to update for myid==0
                   if (.not. c_lnL%pol_sample_nprop(p,id)) then
                      loop_exit=.true. !we have found the ideal proposal length, not looking for correlation length
                   end if
                end if
             end if

          else if (c_lnL%pol_sample_nprop(p,id)) then
             ! if we are keeping track of spec inds for corr. length, then save current spec ind
             n_corr_prop = n_corr_prop + 1
             theta_corr_arr(n_corr_prop) = old_theta
             if (.not. burned_in) then
                n_corr_prop = 0
             else if (n_corr_prop > n_corr_limit) then
                corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                if (corr_len > 0) then ! .or. n_corr_prop > 4*c_lnL%nprop_uni(2,id)) then
                   !set pixreg_nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                   c_lnL%pol_nprop(id)%p%map(:,p) = 1.d0*min(max(2*corr_len,c_lnL%nprop_uni(1,id)),c_lnL%nprop_uni(2,id))
                   loop_exit=.true. !we have found the correlation length
                   c_lnL%pol_sample_nprop(p,id)=.false. !making sure we dont go back into corr. len. sampling
                end if
             end if
          else
             accept_rate = n_accept*1.d0/j
             if (mod(j,out_every)==0) write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate
          end if
          
       end if

       call wall_time(t2)

       if (info%myid == 0 .and. mod(j,out_every)==0) write(*,*) '   Sample:',j,'  Wall time per sample:',real((t2-t1)/j,sp)

       call mpi_bcast(loop_exit, 1, MPI_LOGICAL, 0, c_lnL%comm, ierr)

       if (loop_exit) exit
       
    end do !nprop

    !bcast proposal length
    call mpi_bcast(c_lnL%pol_proplen(id)%p%map(0,p), 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    c_lnL%pol_proplen(id)%p%map(:,p) = c_lnL%pol_proplen(id)%p%map(0,p)

    !bcast number of proposals
    call mpi_bcast(c_lnL%pol_nprop(id)%p%map(0,p), 1, MPI_DOUBLE_PRECISION, 0, &
         & c_lnL%comm, ierr)
    c_lnL%pol_nprop(id)%p%map(:,p) = c_lnL%pol_nprop(id)%p%map(0,p)

    !bcast last valid theta to all procs to update theta map
    call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, 0, c_lnL%comm, ierr)

    buffer_lnL(:,p_min:p_max) = old_theta

    if (info%myid==0) then
       write(*,fmt='(a, f10.5)') "    Final spec. ind. value: ", old_theta
       write(*,fmt='(a, i5, a, e14.5)') "    Number of proposals: ",IDINT(c_lnL%pol_nprop(id)%p%map(0,p)), &
            & "  --  Proposal length: ", c_lnL%pol_proplen(id)%p%map(0,p)
       write(*,fmt='(a, e14.5)') "    Difference in spec. ind., new - old: ", old_theta-init_theta
       write(*,*) ''
    end if

    if (iter >= 2) then !only stop sampling after 2nd iteration, in case first iter is far off.
       c_lnL%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals) and
       c_lnL%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done in the first 
                                               !gibbs iteration anyways)
    end if

    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)
    if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)

  end subroutine sampleDiffuseSpecIndFullsky_nonlin


  subroutine sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    real(dp),               dimension(0:,:), intent(inout)        :: buffer_lnL
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, par_id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, m, n, q, pr, max_pr, pix, ierr, ind(1), counter, n_ok, id
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps, nsamp
    real(dp)     :: a, b, a_tot, b_tot, s, t0, t1, t2, t3, t4, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: buff_init, theta_min, theta_max, running_accept, running_dlnL
    logical(lgt) :: ok, outside
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    !Following is for the local sampler
    real(dp)     :: mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, lnL_prior, lnL_init
    real(dp)     :: accept_rate, accept_scale, lnL_sum, proplen, chisq_jeffreys, avg_dlnL, lnL_total_init
    real(dp)     :: old_theta, new_theta
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop, band_count, pix_count, buff1_i(1), buff2_i(1), burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    integer(i4b) :: npixreg, smooth_scale, arr_ind, np_lr, np_fr, myid_pix, unit
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in, sampled_nprop, sampled_proplen
    character(len=512) :: filename, postfix, fmt_pix, npixreg_txt
    character(len=6) :: itext
    character(len=4) :: ctext
    character(len=2) :: pind_txt
    real(dp),      allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_old_arr, mixing_new_arr
    real(dp),      allocatable, dimension(:) :: theta_corr_arr, old_thetas, new_thetas, init_thetas, sum_theta
    real(dp),      allocatable, dimension(:) :: old_theta_smooth, new_theta_smooth, dlnL_arr
    real(dp),      allocatable, dimension(:,:) :: theta_MC_arr
    integer(i4b),  allocatable, dimension(:) :: band_i, pol_j, accept_arr
    class(comm_mapinfo),             pointer :: info_fr => null() !full resolution
    class(comm_mapinfo),             pointer :: info_fr_single => null() !full resolution, nmaps=1
    class(comm_mapinfo),             pointer :: info_lr => null() !low resolution
    class(comm_mapinfo),             pointer :: info_lr_single => null() !low resolution, nmaps=1
    class(comm_mapinfo),             pointer :: info_data => null() !data band map info
    class(comm_map),                 pointer :: theta_single_fr => null() ! Spectral parameter of polt_id index
    class(comm_map),                 pointer :: theta_single_lr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_lr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_lr_hole => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_fr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: mask_lr => null() ! lowres mask
    class(comm_map),                 pointer :: res_map => null() ! lowres  residual map
    class(comm_map),                 pointer :: temp_res => null() 
    type(map_ptr), allocatable, dimension(:) :: lr_chisq
    type(map_ptr), allocatable, dimension(:) :: df

    c           => compList     ! Extremely ugly hack...
    do while (comp_id /= c%id)
       c => c%next()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select

    id = par_id !hack to not rewrite too much from diffuse_comp_mod

    info_fr  => comm_mapinfo(c_lnL%theta(id)%p%info%comm, c_lnL%theta(id)%p%info%nside, &
         & 2*c_lnL%theta(id)%p%info%nside, c_lnL%theta(id)%p%info%nmaps, c_lnL%theta(id)%p%info%pol)

    info_fr_single  => comm_mapinfo(c_lnL%theta(id)%p%info%comm, c_lnL%theta(id)%p%info%nside, &
         & 2*c_lnL%theta(id)%p%info%nside, 1, c_lnL%theta(id)%p%info%pol)

    info_lr  => comm_mapinfo(c_lnL%x_smooth%info%comm, c_lnL%x_smooth%info%nside, &
         & c_lnL%x_smooth%info%lmax, c_lnL%x_smooth%info%nmaps, c_lnL%x_smooth%info%pol)

    info_lr_single  => comm_mapinfo(c_lnL%x_smooth%info%comm, c_lnL%x_smooth%info%nside, &
         & c_lnL%x_smooth%info%lmax, 1, c_lnL%x_smooth%info%pol)
                
    theta_single_fr => comm_map(info_fr_single)
    theta_fr => comm_map(info_fr_single)

    myid_pix = info_fr%myid ! using this proc id and info_fr%comm in all mpi calls to ensure that values are properly dispersed between processors 

    !ud_grade mask
    mask_lr => comm_map(info_lr)
    call c_lnL%pol_ind_mask(id)%p%udgrade(mask_lr)

    !init lowres residual map
    res_map => comm_map(info_lr)


    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 20
    burned_in           = .false.
    first_call          = .false.

    
    n_prop_limit = 40
    n_corr_limit = 100
    out_every    = 40
   
    !     buffer_lnL (theta, limited by min/max)

    if (c_lnL%poltype(id) == 1) then
       p_min = 1; p_max = c_lnL%nmaps
       if (cpar%only_pol) p_min = 2
    else if (c_lnL%poltype(id) == 2) then
       if (p == 1) then
          p_min = 1; p_max = 1
       else
          p_min = 2; p_max = c_lnL%nmaps
       end if
    else if (c_lnL%poltype(id) == 3) then
       p_min = p
       p_max = p
    else
       write(*,*) 'Unsupported polarization type'
       stop
    end if

    npar      = c_lnL%npar
    np_fr     = info_fr%np
    np_lr     = info_lr%np
    nmaps     = c_lnL%x_smooth%info%nmaps
    npixreg   = c_lnL%npixreg(p,id)

    theta_min = c_lnL%p_uni(1,id)
    theta_max = c_lnL%p_uni(2,id)


    !set up which bands and polarizations to include
    allocate(band_i(3*numband),pol_j(3*numband))

    band_count=0
    do k = 1,numband !run over all active bands
       !check if the band is associated with the smoothed component, and if band frequencies are within 
       !freq. limits for the component
       if (.not. associated(rms_smooth(k)%p)) cycle
       if (data(k)%bp(0)%p%nu_c < c_lnL%nu_min_ind(id) .or. data(k)%bp(0)%p%nu_c > c_lnL%nu_max_ind(id)) cycle
                   
       do l = p_min,p_max !loop all maps of band k associated with p (given poltype)
          if (l > data(k)%info%nmaps) cycle !no map in band k for polarization l
          band_count = band_count + 1
          band_i(band_count)=k
          pol_j(band_count)=l
       end do
    end do

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set theta to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (myid_pix == 0)  write(*,*) 'no data bands available for sampling of spec ind'
       return
    else
       if (cpar%verbosity>3 .and. myid_pix == 0) then
             write(*,fmt='(a)') '  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') '   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do
       end if
       if (trim(c_lnL%pol_lnLtype(p,id))=='chisq' .and. .true.) then !debug chisq (RMS scaling) for smoothing scale
          allocate(lr_chisq(band_count))
          do k = 1,band_count
             lr_chisq(k)%p => comm_map(info_lr_single)
             lr_chisq(k)%p%map=0.d0
          end do
       end if
    end if

    ! This is used for marginal/ridge sampling
    allocate(all_thetas(npar))
    allocate(mixing_old_arr(band_count),mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))

    !This is used for all
    allocate(old_thetas(0:npixreg),new_thetas(0:npixreg),init_thetas(0:npixreg))
    allocate(old_theta_smooth(0:np_lr-1), new_theta_smooth(0:np_lr-1))
    allocate(accept_arr(n_prop_limit), dlnL_arr(n_prop_limit))

    !This is used for fullres chisq
    if (c_lnL%apply_jeffreys) then
       allocate(df(numband))
       do k = 1, numband
          df(k)%p => comm_map(data(k)%info)
       end do
    end if

    sampled_nprop = .false.
    sampled_proplen = .false.
    if (c_lnL%pol_sample_nprop(p,id)) sampled_nprop = .true.
    if (c_lnL%pol_sample_proplen(p,id)) sampled_proplen = .true.

    old_thetas = c_lnL%theta_pixreg(:,p,id)
    old_thetas = min(max(old_thetas,theta_min),theta_max)
    new_thetas = old_thetas
    ! that the root processor operates on
    init_thetas = old_thetas
    if (cpar%verbosity>2 .and. info_fr%myid == 0 .and. npixreg > 0) write(*,fmt='(a, f10.5)') "  initial (avg) spec. ind. value: ", &
         & sum(init_thetas(1:npixreg))/npixreg

    if (myid_pix==0) then
       allocate(theta_MC_arr(10000,npixreg))
       theta_MC_arr = 0.d0
    end if

    do pr = 1,npixreg

       call wall_time(t0)
       !debug
       !if (myid_pix==0) then
       !   write(*,*) myid_pix,info_lr%myid,'init',init_thetas
       !   write(*,*) myid_pix,info_lr%myid,'old',old_thetas
       !end if
       !if (myid_pix==1) then
       !   write(*,*) myid_pix,info_lr%myid,'old',old_thetas
       !end if

       nsamp=0
       n_spec_prop = 0
       n_accept = 0
       n_corr_prop = 0 
       avg_dlnL = 0.d0
       accept_rate=0.d0
       accept_arr(:)=0
       dlnL_arr(:)=0.d0
       running_dlnL=1.d3
       running_accept=0.d0

       theta_fr%map = old_thetas(0) !prior value
       theta_single_fr%map = 0.d0

       if (sampled_nprop) c_lnL%pol_sample_nprop(p,id) = .true.
       if (sampled_proplen) c_lnL%pol_sample_proplen(p,id) = .true.

       if (c_lnL%npix_pixreg(pr,p,id) == 0) cycle !no pixels in pixel region

       pix_count = 0
       do pix = 0,np_fr-1
          if (pr == c_lnL%ind_pixreg_arr(pix,p,id)) then
             if (c_lnL%pol_ind_mask(id)%p%map(pix,p) > 0.5d0) pix_count = pix_count + 1
          end if
       end do

       call mpi_allreduce(MPI_IN_PLACE, pix_count, 1, MPI_INTEGER, & 
            & MPI_SUM, info_fr%comm, ierr)
       if (pix_count == 0) cycle !all pixels in pixreg is masked out
       
       if (c_lnL%pol_sample_nprop(p,id) .or. c_lnL%pol_sample_proplen(p,id)) then
          pixreg_nprop = 10000 !should be enough to find correlation length (and proposal length if prompted) 
          allocate(theta_corr_arr(pixreg_nprop))
          !c_lnL%pol_sample_nprop(j,p,id) = boolean array of size (n_pixreg,poltype)
          !c_lnL%pol_sample_proplen(j,p,id) = boolean array of size (n_pixreg,poltype)
       else
          pixreg_nprop = c_lnL%nprop_pixreg(pr,p,id)
       end if

       first_sample  = .true.
       loop_exit     = .false.

       !assign values to pixel regions that are not in sampled pixreg, those in pix.reg. are set to 1.d0 in single map
       !write(*,*) info_fr%myid,c_lnL%pol_ind_mask(id)%p%info%myid, np_fr, c_lnL%pol_ind_mask(id)%p%info%np
       do pix = 0,np_fr-1
          if (pr /= c_lnL%ind_pixreg_arr(pix,p,id)) then
             theta_fr%map(pix,1) = old_thetas(c_lnL%ind_pixreg_arr(pix,p,id))
          else
             theta_fr%map(pix,1) = 0.d0
             theta_single_fr%map(pix,1) = 1.d0
          end if
       end do

       !then we smooth to lower resolution
       smooth_scale = c_lnL%smooth_scale(id)
       if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
          if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then !smooth to correct resolution
             call smooth_map(info_lr_single, .false., &
                  & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, theta_fr, &  
                  & data(1)%B_postproc(smooth_scale)%p%b_l, theta_lr_hole)
             call smooth_map(info_lr_single, .false., &
                  & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, theta_single_fr, &  
                  & data(1)%B_postproc(smooth_scale)%p%b_l, theta_single_lr)
          else !no postproc smoothing, ud_grade to correct resolution
             theta_single_lr => comm_map(info_lr_single)
             theta_lr_hole => comm_map(info_lr_single)
             call theta_single_fr%udgrade(theta_single_lr)
             call theta_fr%udgrade(theta_lr_hole)
          end if
       else !native resolution, fr = lr
          theta_single_lr => comm_map(info_lr_single)
          theta_lr_hole => comm_map(info_lr_single)
          theta_single_lr%map(:,1) = theta_single_fr%map(:,1)
          theta_lr_hole%map(:,1) = theta_fr%map(:,1)
       end if

       call wall_time(t1)
       lnl_init=0.d0

       do j = 1,pixreg_nprop !propose new sample n times (for faster mixing in total)

          nsamp = nsamp + 1
          if (myid_pix == 0) then

             if (first_sample) then
                old_theta = old_thetas(pr)
                proplen=c_lnL%proplen_pixreg(pr,p,id)
                if (cpar%verbosity>2) then                
                   write(*,fmt='(a, i3, a, f10.5)') "  initial pix.reg. ",pr,"  -- spec. ind. value: ", init_thetas(pr)
                   write(*,fmt='(a, e14.5)') "  initial proposal length: ",proplen
                end if
                new_thetas = old_thetas
             end if

             new_theta = old_theta + proplen*rand_gauss(handle)
          end if

          !broadcast new_theta (and old_theta)
          if (first_sample) then
             call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, &
               & 0, info_fr%comm, ierr)
             old_thetas(pr)=old_theta
          end if

          call mpi_bcast(new_theta, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
          new_thetas(pr) = new_theta

          if (cpar%verbosity>3 .and. myid_pix==0 .and. mod(j,out_every)==0) then
             write(*,fmt='(a, i6, a, i3, a, f10.5, a, f10.5)') "  proposal: ", j," -- Pixreg ", pr, &
                  & " -- Current ind: ", old_theta, " -- proposed ind: ", new_theta
          end if


          !init log-like for new sample
          if (first_sample) lnL_old = 0.d0
          lnL_new = 0.d0       

          outside = .false.

          if (new_theta > theta_max .or. new_theta < theta_min) outside = .true.

          if (outside) then
             lnL_new = -1.d30 
             ! skip the true caclulation of lnL, we reject the sample ~100%
             if (myid_pix==0 .and. cpar%verbosity > 2) then
                write(*,fmt='(a, f10.5, a, f10.5)') "    Proposed ind outside limits.  min: ", &
                     & theta_min," -- max: ", theta_max
                write(*,fmt='(a, f10.5)') "    Proposed ind ",new_theta
             end if
          else

             if (first_sample) then
                !set up the old theta map
                old_theta_smooth(:)=theta_lr_hole%map(:,1)+ old_theta*theta_single_lr%map(:,1)
             end if
             !set up the new theta map
             new_theta_smooth(:)=theta_lr_hole%map(:,1)+ new_theta*theta_single_lr%map(:,1)

             !lnL type should split here
             if (trim(c_lnL%pol_lnLtype(p,id))=='chisq') then

                do k = 1,band_count !run over all active bands

                   if (data(band_i(k))%N%type == "rms") then !normal chisq

                      do pix = 0,np_lr-1 !loop over pixels covered by the processor (on the lowres (smooth scale) map)
                         if (mask_lr%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out, go to next pixel
                         all_thetas(id)=new_theta_smooth(pix)
                         !get the values of the remaining spec inds of the component for the given pixel
                         do i = 1, npar
                            if (i == id) cycle
                            all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                         end do

                         ! get conversion factor from amplitude to data (i.e. mixing matrix element)           
                         ! both for old and new spec. ind. and calc. log likelihood (chisq)
                         if (first_sample) then
                            all_thetas(id)=old_theta_smooth(pix)
                            mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                                 & data(band_i(k))%gain * c_lnL%cg_scale
                            all_thetas(id)=new_theta_smooth(pix)

                            res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old* &
                                 & c_lnL%x_smooth%map(pix,pol_j(k))
                            lnL_old = lnL_old -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2

                         end if

                         mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                              & data(band_i(k))%gain * c_lnL%cg_scale

                         res_lnL = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                              & c_lnL%x_smooth%map(pix,pol_j(k))

                         if (allocated(lr_chisq)) then !for debugging
                            lr_chisq(k)%p%map(pix,1) = (res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                         end if

                         lnL_new = lnL_new -0.5d0*(res_lnL*rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k)))**2
                      end do

                   else if (data(band_i(k))%N%type == "QUcov") then !QU_cov rms, udgrade if necessary
                      if (first_sample) then
                         !build residual map
                         res_map%map = 0.d0
                         do pix = 0,np_lr-1 !loop over pixels covered by the processor
                            if (mask_lr%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out, go to next pixel

                            !get the values of the remaining spec inds of the component for the given pixel
                            do i = 1, npar
                               if (i == id) cycle
                               all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                            end do

                            all_thetas(id)=old_theta_smooth(pix)

                            mixing_old = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                                 & data(band_i(k))%gain * c_lnL%cg_scale

                            res_map%map(pix,pol_j(k)) = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_old* &
                              & c_lnL%x_smooth%map(pix,pol_j(k))
                         end do

                         if (data(band_i(k))%info%nside /= info_lr%nside) then
                            info_data  => comm_mapinfo(data(band_i(k))%info%comm, data(band_i(k))%info%nside, &
                                 & 0, data(band_i(k))%info%nmaps, data(band_i(k))%info%nmaps==3)
                            temp_res => comm_map(info_data)
                            call res_map%udgrade(temp_res)
                         else
                            temp_res => comm_map(info_lr)
                            temp_res%map = res_map%map
                         end if

                         call data(i)%N%sqrtInvN(temp_res)

                         temp_res%map = temp_res%map*temp_res%map
                        
                         lnL_old = lnL_old -0.5d0*sum(temp_res%map(:,pol_j(k)))
                         call temp_res%dealloc()
                      end if
                      
                      !build residual map
                      res_map%map = 0.d0
                      do pix = 0,np_lr-1 !loop over pixels covered by the processor
                         if (mask_lr%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out, go to next pixel

                         do i = 1, npar
                            if (i == id) cycle
                            all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                         end do

                         all_thetas(id)=new_theta_smooth(pix)

                         mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                              & data(band_i(k))%gain * c_lnL%cg_scale

                         res_map%map(pix,pol_j(k)) = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                              & c_lnL%x_smooth%map(pix,pol_j(k))
                      end do

                      if (data(band_i(k))%info%nside /= info_lr%nside) then
                         info_data  => comm_mapinfo(data(band_i(k))%info%comm, data(band_i(k))%info%nside, &
                              & 0, data(band_i(k))%info%nmaps, data(band_i(k))%info%nmaps==3)
                         temp_res => comm_map(info_data)
                         call res_map%udgrade(temp_res)
                      else
                         temp_res => comm_map(info_lr)
                         temp_res%map = res_map%map
                      end if

                      call data(i)%N%sqrtInvN(temp_res)

                      temp_res%map = temp_res%map*temp_res%map

                      lnL_new = lnL_new -0.5d0*sum(temp_res%map(:,pol_j(k)))
                      call temp_res%dealloc()

                   end if
                end do

             else if ((trim(c_lnL%pol_lnLtype(p,id))=='ridge') .or. &
                  & (trim(c_lnL%pol_lnLtype(p,id))=='marginal')) then
                if (trim(c_lnL%pol_lnLtype(p,id))=='ridge') then
                   use_det = .false.
                else
                   use_det = .true.
                end if

                do pix = 0,np_lr-1 !More precise, we need to loop over pixels covered by the processor

                   if (mask_lr%map(pix,p) < 0.5d0) cycle     ! if pixel is masked out

                   all_thetas(id)=new_theta_smooth(pix)
                   !get the values of the remaining spec inds of the component for the given pixel
                   do i = 1, npar
                      if (i == id) cycle
                      all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                   end do

                   !build mixing matrix
                   do k = 1,band_count !run over all active bands
                      pix_count = pix_count + 1

                      if (first_sample) then
                         all_thetas(id)=old_theta_smooth(pix)
                         mixing_old_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                              & data(band_i(k))%gain * c_lnL%cg_scale
                         all_thetas(id)=new_theta_smooth(pix)
                      end if
                      mixing_new_arr(k) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale
                      pix_count = pix_count + 1
                      data_arr(k)=res_smooth(band_i(k))%p%map(pix,pol_j(k))
                      if (data(band_i(k))%N%type == "rms") then !normal diagonal noise
                         invN_arr(k)=rms_smooth(band_i(k))%p%siN%map(pix,pol_j(k))**2 !assumed diagonal and uncorrelated 
                      else
                         invN_arr(k)=0.d0 ! Neglect the band (should just get the diagonal element of the cov matrix)
                      end if
                   end do

                   !compute the marginal/ridge log-like for the pixel
                   if (first_sample) then
                      lnL_old = lnL_old + comp_lnL_marginal_diagonal(mixing_old_arr, invN_arr, data_arr, &
                           & use_det, band_count)
                   end if
                   lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                        & use_det, band_count)
                end do
             else 
                write(*,*) 'invalid polarized lnL sampler type'
                stop

             end if !chisq/marginal/ridge


             !Sum lnL_new (and lnL_old) among all processors
             if (first_sample) then
                call mpi_allreduce(MPI_IN_PLACE, lnL_old, 1, MPI_DOUBLE_PRECISION, & 
                     & MPI_SUM, info_fr%comm, ierr)

                if (.not. ( (trim(c_lnL%pol_lnLtype(p,id))=='chisq') .or. (trim(c_lnL%pol_lnLtype(p,id))=='chisq_lowres'))) then
                end if
             end if

             call mpi_allreduce(MPI_IN_PLACE, lnL_new, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info_fr%comm, ierr)

             !add spec. ind. prior
             if (c_lnL%p_gauss(2,id) > 0.d0) then
                !Find prior "chisq" and add it to lnL
                if (first_sample) then
                   lnL_prior = (old_theta-c_lnL%p_gauss(1,id))**2
                   lnL_prior = -0.5d0 * lnl_prior/c_lnL%p_gauss(2,id)**2
                   lnL_old = lnL_old + lnL_prior
                   lnl_init = lnL_old
                end if
                lnL_prior = (new_thetas(pr)-c_lnL%p_gauss(1,id))**2
                lnL_prior = -0.5d0 * lnl_prior/c_lnL%p_gauss(2,id)**2
                lnL_new = lnL_new + lnL_prior
             end if

             !first sample done
             if (first_sample .and. myid_pix == 0 .and. cpar%verbosity > 2) write(*,fmt='(a, e14.5)') "    lnL_init = ", lnL_init
             first_sample = .false.

          end if !new_theta outside spec limits

          if (myid_pix == 0) then
             !accept/reject new spec ind
             delta_lnL = lnL_new-lnL_old

             if (cpar%verbosity>3 .and. mod(j,out_every)==0) then
                write(*,fmt='(a, e14.5)') "    lnL_new = ", lnL_new
                write(*,fmt='(a, e14.5)') "    lnL_old = ", lnL_old
                write(*,fmt='(a, e14.5)') "    lnL_new - lnL_old = ", delta_lnL
             end if

             avg_dlnL = (avg_dlnL*(n_spec_prop) + abs(delta_lnL))/(n_spec_prop+1) 

             n_spec_prop = n_spec_prop + 1
             arr_ind = mod(n_spec_prop,n_prop_limit)
             if (arr_ind==0) arr_ind=n_prop_limit
             dlnL_arr(arr_ind) = abs(delta_lnL)
             running_dlnL = sum(dlnL_arr(1:min(n_spec_prop,n_prop_limit)))/max(min(n_spec_prop,n_prop_limit),1)

             if (delta_lnL < 0.d0) then
                !if delta_lnL is more negative than -25.d0, limit to -25.d0
                if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

                !draw random uniform number
                a = rand_uni(handle) !draw uniform number from 0 to 1
                if (exp(delta_lnL) > a) then
                   !accept
                   old_theta = new_theta
                   lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                   n_accept = n_accept + 1
                   accept_arr(arr_ind) = 1 !accept
                else
                   accept_arr(arr_ind) = 0 !reject
                end if
             else
                !accept new sample, higher likelihood
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
                accept_arr(arr_ind) = 1 !accept
             end if
             
             theta_MC_arr(j,pr) = old_theta

             running_accept = (1.d0*sum(accept_arr(1:min(n_spec_prop,n_prop_limit))))/max(min(n_spec_prop,n_prop_limit),1)

             if (j > burn_in) burned_in = .true.
             ! evaluate proposal_length, then correlation length. 
             !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
             accept_rate = n_accept*1.d0/n_spec_prop
             if (c_lnL%pol_sample_proplen(p,id)) then
                if (cpar%verbosity>3 .and. mod(n_spec_prop,out_every)==0) then
                   write(*,fmt='(a, f6.4)') "   accept rate = ", running_accept
                   write(*,fmt='(a, e14.5)') "   avg. abs. delta_lnL = ", running_dlnL
                   write(*,fmt='(a, e14.5)') "    lnL_new = ", lnL_new
                   write(*,fmt='(a, e14.5)') "    lnL_old = ", lnL_old
                   write(*,fmt='(a, e14.5)') "    lnL_new - lnL_old = ", delta_lnL
                end if
                
                if (.not. burned_in) then 
                   !giving som time to let the first steps "burn in" if one starts far away from max lnL
                   n_spec_prop = 0
                   n_accept = 0
                   avg_dlnL = 0.d0
                else if (n_spec_prop > n_prop_limit) then
                   if (running_accept > 0.7d0 .and. running_dlnL < 10.d0) then
                      accept_scale = 1.5d0
                   else if (running_accept < 0.3d0) then
                      accept_scale = 0.5d0
                   else 
                      accept_scale = -1.d0
                   end if

                   if (accept_scale > 0.d0) then
                      !prop_len only ever used by root processor
                      proplen = proplen*accept_scale !set new proposal length
                      n_accept = 0 !reset with new prop length
                      n_spec_prop = 0 !reset with new prop length
                      avg_dlnL = 0.d0
                      if (cpar%verbosity>3) then
                         write(*,fmt='(a, f6.4)')  "      accept rate =         ", running_accept
                         write(*,fmt='(a, e14.5)') "      avg. abs. delta_lnL = ", running_dlnL
                         write(*,fmt='(a, e14.5)') "      New prop. len. =      ", proplen
                      end if
                   else
                      !add additional requirement that avg. absolute delta chisq less than 10
                      ! else continue running (have not gotten close enough to Max lnL)
                      if (running_dlnL < 10.d0) then
                         c_lnL%proplen_pixreg(pr,p,id) = proplen
                         c_lnL%pol_sample_proplen(p,id)=.false. !making sure we dont go back into prop. len. sampling
                         !only necessary to update for myid==0
                         if (.not. c_lnL%pol_sample_nprop(p,id)) then
                            loop_exit=.true. !we have found the ideal proposal length, not looking for correlation length
                         end if
                      end if
                   end if
                end if

             else if (c_lnL%pol_sample_nprop(p,id)) then
                ! if we are keeping track of spec inds for corr. length, then save current spec ind
                n_corr_prop = n_corr_prop + 1
                theta_corr_arr(n_corr_prop) = old_theta
                if (.not. burned_in) then
                   n_corr_prop = 0
                else if (n_corr_prop > n_corr_limit) then
                   corr_len = calc_corr_len(theta_corr_arr(1:n_corr_prop),n_corr_prop) !returns -1 if no corr.len. is found
                   if (corr_len > 0) then ! .or. n_corr_prop > 4*c_lnL%nprop_uni(2,id)) then
                      !set pixreg_nprop (for pix region) to x*corr_len, x=2, but inside limits from parameter file
                      c_lnL%nprop_pixreg(pr,p,id) = 1.d0*min(max(2*corr_len,c_lnL%nprop_uni(1,id)),c_lnL%nprop_uni(2,id))
                      loop_exit=.true. !we have found the correlation length
                      c_lnL%pol_sample_nprop(p,id)=.false. !making sure we dont go back into corr. len. sampling
                   end if
                end if
             else
                if (cpar%verbosity>2 .and. mod(j,out_every)==0) then
                   write(*,fmt='(a, f6.4)') "   accept rate = ", accept_rate                   
                   if (cpar%verbosity>3) write(*,fmt='(a, f6.4)') "   avg. abs. delta_lnL = ", avg_dlnL
                end if
             end if

          end if !myid_pix == 0

          call wall_time(t2)

          if (cpar%verbosity>3 .and. myid_pix == 0 .and. mod(j,out_every)==0) then
             write(*,*) '   Sample:',j,'  Wall time per sample:',real((t2-t1)/j,sp)
          end if

          call mpi_bcast(loop_exit, 1, MPI_LOGICAL, 0, info_fr%comm, ierr)

          if (loop_exit) exit

       end do !nprop

       if (pr == 1) lnl_total_init = lnl_init

       call wall_time(t2)

       !save old theta for further pixelregions
       call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
       old_thetas(pr)=old_theta

       if (cpar%verbosity>2) then !pixreg pr
          if (myid_pix==0) then
             write(*,fmt='(a, i5)')    "    Pixel region: ",pr
             write(*,fmt='(a, f10.5)') "      Final spec. ind. value:                   ", old_thetas(pr)
             write(*,fmt='(a, e14.5)') "      Difference in spec. ind., new - old:      ", &
                  & (old_thetas(pr)-init_thetas(pr))
             write(*,*) '      Samples:',nsamp
             if (nsamp > 0) write(*,*) '        Wall time per sample (sec):                 ',real((t2-t1)/nsamp,sp)
             write(*,*) '        Initialization wall time pixel region (sec):',real((t1-t0),sp)
             write(*,*) '        Total wall time pixel region (sec):         ',real((t2-t0),sp)

             if (sampled_nprop) write(*,fmt='(a, i5)') "      Number of proposals after tuning: ",c_lnL%nprop_pixreg(pr,p,id)
             if (sampled_proplen) write(*,fmt='(a, e14.5)') "      Proposal length after tuning: ",c_lnL%proplen_pixreg(pr,p,id)
             write(*,fmt='(a, e14.5)') "      New Log-Likelihood:                       ", lnl_old
             write(*,fmt='(a, e14.5)') "      Difference in Log-Likelihood (new - old): ", lnl_old-lnl_init
             write(*,*) ''
          end if
       end if

       if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)
       call theta_single_lr%dealloc()
       call theta_lr_hole%dealloc()

    end do !pr = 1,max_pr

    !bcast proposal length
    call mpi_bcast(c_lnL%proplen_pixreg(1:npixreg,p,id), npixreg, MPI_DOUBLE_PRECISION, 0, &
         & info_fr%comm, ierr)
    do pix = 0,np_fr-1
       c_lnL%pol_proplen(id)%p%map(pix,p) = c_lnL%proplen_pixreg(c_lnL%ind_pixreg_arr(pix,p,id),p,id)
    end do
    
    !bcast number of proposals
    call mpi_bcast(c_lnL%nprop_pixreg(1:npixreg,p,id), npixreg, MPI_INTEGER, 0, &
         & info_fr%comm, ierr)
    do pix = 0,np_fr-1
       c_lnL%pol_nprop(id)%p%map(pix,p) = 1.d0*c_lnL%nprop_pixreg(c_lnL%ind_pixreg_arr(pix,p,id),p,id)
    end do

    !bcast last valid theta to all procs to update theta map
    call mpi_bcast(old_thetas(0:npixreg), npixreg+1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)

    !assign thetas to pixel regions, thetas will be smoothed to same FWHM as in lnL eval when exiting sampler 
    do pix=0,np_fr-1
       theta_fr%map(pix,1)=old_thetas(c_lnL%ind_pixreg_arr(pix,p,id))
    end do

    do i = p_min,p_max
       buffer_lnL(:,i) = theta_fr%map(:,1)
    end do
    c_lnL%theta_pixreg(0:npixreg,p,id)=old_thetas

    if (cpar%verbosity>2 .and. myid_pix==0 .and. npixreg > 1) then
       write(*,*) "    Average values from pixel region sampling"
       write(*,fmt='(a, i5)') "    Number of proposals (after tuning):  ", &
            & sum(c_lnL%nprop_pixreg(1:npixreg,p,id))/npixreg
       write(*,fmt='(a, e14.5)') "  --  Proposal length (after tuning): ", &
            & sum(c_lnL%proplen_pixreg(1:npixreg,p,id))/npixreg
       write(*,fmt='(a, f10.5)') "    Final spec. ind. value:           ", sum(old_thetas(1:npixreg))/npixreg
       write(*,fmt='(a, e14.5)') "    Difference in spec. ind., new - old: ", &
            & sum(old_thetas(1:npixreg)-init_thetas(1:npixreg))/npixreg
       write(*,fmt='(a, e14.5)') "    New Log-Likelihood:               ", &
            & lnl_old
       write(*,fmt='(a, e14.5)') "    Difference in Log-Likelihood (new - old): ", &
            & lnl_old-lnl_total_init
       write(*,*) ''
    end if

    !debug
    if (allocated(lr_chisq)) then
       call int2string(cpar%mychain, ctext)
       call int2string(iter,         itext)

       do k = 1,band_count
          filename=trim(cpar%outdir)//'/'//'chisq_local_'//trim(data(band_i(k))%label)

          if (pol_j(k) == 1) then
             filename=trim(filename)//'_I'
          else if (pol_j(k) == 2) then
             filename=trim(filename)//'_Q'
          else if (pol_j(k) == 3) then
             filename=trim(filename)//'_U'             
          end if
          filename=trim(filename)//'_'//trim(postfix)//'.fits'

          call lr_chisq(k)%p%writeFITS(trim(filename))
          call lr_chisq(k)%p%dealloc()
       end do
       deallocate(lr_chisq)
    end if

    if (iter >= c_lnL%sample_first_niter) then !only stop sampling after n'th iteration, in case first iter is far off.
       c_lnL%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals) and
       c_lnL%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done in the first 
                                               !gibbs iteration anyways)
    else !reset sampling to true if we have sampled
       if (sampled_nprop) c_lnL%pol_sample_nprop(p,id) = .true. 
       if (sampled_proplen) c_lnL%pol_sample_proplen(p,id) = .true. 
    end if

    call int2string(iter,         itext)
    call int2string(p,         pind_txt)
    call int2string(cpar%mychain, ctext)
    postfix = 'c'//ctext//'_k'//itext//'_p'//pind_txt

    !print MC theta to file
    if (myid_pix==0) then
       unit = getlun()
       filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//trim(c_lnL%indlabel(id))//&
            & '_theta_MC_'//trim(postfix)//'.dat'
       open(unit,file=trim(filename))
       !read(npixreg_txt,*) npixreg
       !fmt_pix=trim(npixreg_txt)//'f12.6'
       !write(*,*) fmt_pix
       do i = 1,10000
          !write(unit,'(i8,'//trim(fmt_pix)//')') i,theta_MC_arr(i,:)
          write(unit,'(i8,f12.6)') i,theta_MC_arr(i,1)
       end do
       close(unit)
       deallocate(theta_MC_arr)
    end if

    if (.true.) then
       do i = 1,band_count
          filename=trim(cpar%outdir)//'/'//'reduced_data_band_'//trim(data(band_i(i))%label)//'_'// &
               & trim(c_lnl%label)//'_'//trim(c_lnL%indlabel(id))//'_'//trim(postfix)//'.fits'
          call res_smooth(band_i(i))%p%writeFITS(trim(filename))
       end do
    end if

    !deallocate
    deallocate(mixing_old_arr,mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)
    deallocate(old_thetas,new_thetas,init_thetas,old_theta_smooth, new_theta_smooth)
    if (allocated(df)) deallocate(df)
    deallocate(accept_arr,dlnL_arr)

    call theta_fr%dealloc()
    call theta_single_fr%dealloc()
    call mask_lr%dealloc()
    call res_map%dealloc()
    theta_fr => null()
    theta_single_fr => null()
    mask_lr => null()
    res_map => null()

  end subroutine sampleDiffuseSpecIndPixReg_nonlin

  function calc_corr_len(spec_arr,n_spec) 
    implicit none
    integer(i4b),               intent(in) :: n_spec
    real(dp),     dimension(:), intent(in) :: spec_arr
    integer(i4b)                           :: calc_corr_len

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y
    real(dp), dimension(:), allocatable :: xarr, yarr,covarr

    calc_corr_len = -1 !default if no corr length is found

    maxcorr = n_spec/2
    allocate(xarr(maxcorr),yarr(maxcorr),covarr(maxcorr))

    do i = 1,maxcorr
       ns=maxcorr-i
       xarr=spec_arr(1:ns)
       yarr=spec_arr(1+i:ns+i)
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       covarr(i) = sum((xarr-x_mean)*(yarr-y_mean))/(ns*sig_x*sig_y)
       if (covarr(i) < 0.d0) then
          calc_corr_len = i
          exit
       end if
    end do

    deallocate(xarr,yarr,covarr)

  end function calc_corr_len

  function comp_lnL_marginal_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_marginal_diagonal

    integer(i4b) :: i, j, mat_len
    real(dp)     :: MNd,MNM
    real(dp), dimension(:), allocatable :: MN

    if (present(arr_len)) then
       allocate(MN(arr_len))
       MN=mixing(1:arr_len)*invN_arr(1:arr_len)
       MNd=sum(MN*data(1:arr_len))
       MNM=sum(MN*mixing(1:arr_len))
    else
       allocate(MN(size(mixing)))
       MN=mixing*invN_arr
       MNd=sum(MN*data)
       MNM=sum(MN*mixing)
    end if

    comp_lnL_marginal_diagonal = 0.d0

    if (MNM /= 0.d0) then 
       MNM=1.d0/MNM !invert 1x1 matrix
    else
       comp_lnL_marginal_diagonal=1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    comp_lnL_marginal_diagonal = 0.5d0*MNd*MNM*MNd

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_marginal_diagonal = comp_lnL_marginal_diagonal - 0.5d0*log(MNM) 

    deallocate(MN)


  end function comp_lnL_marginal_diagonal



end module comm_nonlin_mod
