module comm_nonlin_mod
  use comm_param_mod
  use comm_data_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_gain_mod
  use comm_line_comp_mod
  use comm_diffuse_comp_mod
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

  subroutine sample_nonlin_params(cpar, iter, handle)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    

    integer(i4b) :: i
    real(dp)     :: t1, t2
    class(comm_comp),    pointer :: c    => null()

    call wall_time(t1)
    
    ! Sample spectral indices with standard chisquare sampler, conditional on amplitudes
    call sample_specind_alm(cpar, iter, handle)

    ! Sample spectral indices with local sampler (per-pixel, ptsrc, templates)
    call sample_specind_local(cpar, iter, handle)

    ! Sample calibration factors
    do i = 1, numband
       if (.not. data(i)%sample_gain) cycle
       call sample_gain(cpar%operation, i, cpar%outdir, cpar%mychain, iter, handle)
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


  subroutine sample_specind_alm(cpar, iter, handle)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, nalm_tot, ind
    real(dp)     :: t1, t2, ts, steplen, dalm, fwhm_prior, sigma_prior
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext, jtext
    logical :: accepted, exist, doexit
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_comp),    pointer :: c    => null()
    real(dp),          allocatable, dimension(:,:,:)   :: alms, alms_covmat, L
    real(dp),          allocatable, dimension(:,:) :: m
    real(dp),          allocatable, dimension(:) :: buffer, rgs, sigma_priors, chisq

    integer(c_int),    allocatable, dimension(:,:) :: lm
    integer(i4b), dimension(MPI_STATUS_SIZE) :: mpistat

    ! Sample spectral parameters for each signal component
    allocate(status_fit(numband))
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%next()
          cycle
       end if
       if (all(c%p_gauss(2,:) == 0.d0)) then
          c => c%next()
          cycle
       end if
       
       select type (c)
       class is (comm_diffuse_comp)
          
          do j = 1, c%npar
             
             if (c%p_gauss(2,j) == 0.d0 .or. c%lmax_ind >= 0) cycle
             
             ! Set up smoothed data
             if (cpar%myid_chain == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(j))
             call update_status(status, "spec_alm start " // trim(c%label)// ' ' // trim(c%indlabel(j)))
             
             ! TODO
             ! Reading cholesky for every sample.
             ! Steplength globvar
             ! Put in separate routine? Avoid circle import
             ! Burnin-runde
             
             call wall_time(t1)
             ! Params
             nalm_tot = (c%lmax_ind+1)**2
             steplen = 0.3d0 !c%p_gauss(2,j) ! Init learning rate for proposals
             fwhm_prior = 1200.d0
             sigma_prior = 0.01d0
             out_every = 30
             nsamp = 2000 ! Maxsamp
             num_accepted = 0
             doexit = .false.
             
             info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
                  & c%x%info%lmax, c%x%info%nmaps, c%x%info%pol)
             
             allocate(chisq(0:nsamp))
             allocate(alms(0:nsamp, 0:nalm_tot-1,info%nmaps))                         
             allocate(rgs(0:nalm_tot-1)) ! Allocate random vector
             allocate(L(0:nalm_tot-1, 0:nalm_tot-1, info%nmaps))                         
             if (info%myid == 0) then
                if (nalm_tot > 1) then
                   allocate(sigma_priors(1:nalm_tot-1)) !a_00 is given by different one
                   ! Saving and smoothing priors
                   idx = 1
                   do l_ = 1, c%lmax_ind ! Skip a_00 - m-major ordering (0,0)(1,0)(2,0)(1,1)(1,-1)(2,1)(2,-1)(2,2)(2,-2)
                      sigma_priors(idx) = sigma_prior*exp(-0.5d0*l_*(l_+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                      idx = idx + 1
                   end do
                   do l_ = 1, c%lmax_ind
                      do m_ = 1, l_
                         sigma_priors(idx) = sigma_prior*exp(-0.5d0*l_*(l_+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                         idx = idx + 1
                         sigma_priors(idx) = sigma_prior*exp(-0.5d0*l_*(l_+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                         idx = idx + 1
                      end do
                   end do
                end if
                open(69, file=trim(cpar%outdir)//'/nonlin-samples.dat', recl=10000)
                
                ! Read cholesky matrix. Only root needs this
                inquire(file=trim(cpar%datadir)//'/alm_cholesky.dat', exist=exist)
                if (exist .and. trim(c%label)=="synch" ) then ! If present cholesky file
                   write(*,*) "Reading cholesky matrix"
                   open(unit=11, file=trim(cpar%datadir)//'/alm_cholesky.dat', form='unformatted')
                   read(11) steplen
                   read(11) L
                   close(11)
                else
                   write(*,*) "No cholesky matrix found"
                   L(:,:,:) = 0.d0 ! Set diagonal to 0.001
                   L(0,0,:) = c%p_gauss(2,j) ! Set diagonal to 0.001
                   if (nalm_tot > 1) then
                      do p = 1, nalm_tot-1
                         !L(p,p,:) = 0.1*c%p_gauss(2,j) ! Set diagonal to 0.001
                         L(p,p,:) = 0.1d0*sigma_priors(p)
                      end do
                   end if
                end if
             end if
             
             ! Save initial alm        
             alms(:,:,:) = 0.d0
             ! Gather alms from threads to alms array with correct indices
             do pl = 1, c%theta(j)%p%info%nmaps
                call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, 0, pl, pl)
                allocate(buffer(nalm_tot))
                call mpi_allreduce(alms(0,:,pl), buffer, nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                alms(0,:,pl) = buffer
                deallocate(buffer)
             end do
             
             ! Calculate initial chisq
             if (allocated(c%indmask)) then
                call compute_chisq(c%comm, chisq_fullsky=chisq(0), mask=c%indmask)
             else
                call compute_chisq(c%comm, chisq_fullsky=chisq(0))
             end if
             
             call wall_time(t1)
             if (info%myid == 0) then 
                chisq_prior = 0.d0 
                do pl = 1, c%theta(j)%p%info%nmaps
                   ! if sample only pol, skip T
                   if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                   ! p already calculated if larger than poltype ( smart ;) )
                   if (pl > c%poltype(j)) cycle
                      
                   !chisq_prior = chisq_prior + ((alms(0,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                   if (nalm_tot > 1) then
                      do p = 1, nalm_tot-1
                         chisq_prior = chisq_prior + (alms(0,p,pl)/sigma_priors(p))**2
                      end do
                   end if
                   chisq(0) = chisq(0) + chisq_prior
                end do
                
                ! Output init sample
                write(*,fmt='(a, i6, a, f16.2, a, 3f7.2)') "# sample: ", 0, " - chisq: " , chisq(0), " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
             end if
             do i = 1, nsamp
                chisq_prior = 0.d0
                ! Sample new alms (Account for poltype)
                alms(i,:,:) = alms(i-1,:,:)
                do pl = 1, c%theta(j)%p%info%nmaps
                   
!!$                   if (mod(iter,2) == 1 .and. pl == 1) cycle
!!$                   if (mod(iter,2) == 0 .and. pl >  1) cycle
                   
                   ! if sample only pol, skip T
                   if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                   
                   ! p already calculated if larger than poltype ( smart ;) )
                   if (pl > c%poltype(j)) cycle
                   
                   ! Gather alms from threads to alms array with correct indices
                   call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                   
                   ! Send all alms to 0 (Dont allreduce because only root will do calculation)
                   allocate(buffer(nalm_tot))
                   call mpi_reduce(alms(i,:,pl), buffer, nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)
                   alms(i,:,pl) = buffer
                   deallocate(buffer)
                   
                   ! Propose new alms
                   if (info%myid == 0) then
                      do p = 0, nalm_tot-1
                         rgs(p) = rand_gauss(handle)                       
                      end do
                      alms(i,:,pl) = alms(i-1,:,pl) + steplen*matmul(L(:,:,pl), rgs)
!!$                      write(*,*) 'a', alms(i,:,pl)
!!$                      write(*,*) 'b', matmul(L(:,:,pl), rgs)
!!$                      write(*,*) 'c', rgs
                      
                         ! Adding prior
                      ! Currently applying same prior on all signals
                      !chisq_prior = chisq_prior + ((alms(i,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                      if (nalm_tot > 1) then
                         do p = 1, nalm_tot-1
                            chisq_prior = chisq_prior + (alms(i,p,pl)/sigma_priors(p))**2
                         end do
                      end if
                   end if
                   
                   ! Broadcast proposed alms from root
                   allocate(buffer(nalm_tot))
                   buffer = alms(i,:,pl)
                   call mpi_bcast(buffer, nalm_tot, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)                   
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
                         call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, 1)
                      else
                         do q = 2, c%theta(j)%p%info%nmaps
                            alms(i,:,q) = alms(i,:,pl)
                            call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)                            
                         end do
                      end if
                   else if (c%poltype(j) == 3) then ! {T,E,B}
                      call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                   end if
                end do
                
                ! Update mixing matrix with new alms
                call c%updateMixmat
                
                ! Calculate proposed chisq
                if (allocated(c%indmask)) then
                   call compute_chisq(c%comm, chisq_fullsky=chisq(i), mask=c%indmask)
                else
                   call compute_chisq(c%comm, chisq_fullsky=chisq(i))
                end if
                
                ! Accept/reject test
                ! Reset accepted bool
                accepted = .false.
                if (info%myid == 0) then
                   chisq(i) = chisq(i) + chisq_prior
                   !write(*,fmt='(i6,3f12.2)') i, chisq(i), chisq(i-1), alms(i,0,pl)/sqrt(4*pi)
                   if ( chisq(i) > chisq(i-1) ) then                 
                      ! Small chance of accepting this too
                      ! Avoid getting stuck in local mminimum
                      diff = chisq(i-1)-chisq(i)
                      accepted = (rand_uni(handle) < exp(0.5d0*diff))
                   else
                      accepted = .true.
                   end if
                   
                   ! Count accepted and assign chisq values
                   if (accepted) then
                      num_accepted = num_accepted + 1
                   else
                      chisq(i) = chisq(i-1)
                   end if
                end if
                
                ! Broadcast result of accept/reject test
                call mpi_bcast(accepted, 1, MPI_LOGICAL, 0, c%comm, ierr)
                
                if (.not. accepted) then
                   ! If rejected, restore old values and send to 
                   do pl = 1, c%theta(j)%p%info%nmaps
                      call distribute_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, info%lm, i-1, pl, pl)
                   end do
                   alms(i,:,:) = alms(i-1,:,:)
                   call c%updateMixmat                                   
                end if
                
                ! Write to screen
                if (mod(i,out_every) == 0 .and. info%myid == 0) then
                   call wall_time(t2)                   
                   diff = chisq(i-out_every) - chisq(i) ! Output diff
                   ts = (t2-t1)/DFLOAT(out_every) ! Average time per sample
                   write(*,fmt='(a,i6, a, f16.2, a, f10.2, a, f7.2, a, 3f7.2)') "- sample: ", i, " - chisq: " , chisq(i), " - diff: ", diff, " - time/sample: ", ts, " - a_00: ", alms(i,0,:)/sqrt(4.d0*PI)
                   call wall_time(t1)
                end if
                
                ! Adjust learning rate every 100th
                if (mod(i, 30) == 0 .and. info%myid == 0) then 
                   ! Accept rate
                   accept_rate = num_accepted/30.d0
                   num_accepted = 0
                   
                   diff = chisq(i-30)-chisq(i)
                   
                   ! Write to screen
                   write(*, fmt='(a, i6, a, f8.2, a, f5.3)') "# sample: ", i, " - diff last 30: ", diff, " - accept rate: ", accept_rate
                   
                   ! Adjust steplen
                   if (accept_rate < 0.2.and. iter == 1) then                 
                      steplen = steplen*0.5d0
                      if (info%myid == 0) write(*,fmt='(a,f10.5)') "Reducing steplen -> ", steplen
                   else if (accept_rate > 0.8 .and. iter == 1) then
                      steplen = steplen*2.d0
                      if (info%myid == 0) write(*,fmt='(a,f10.5)') "Increasing steplen -> ", steplen
                   end if
                   
                   ! Burnin
                   if (iter == 1 .and. diff < 20.d0 .and. i > 1000) doexit = .true.
                   if (iter  >= 1 .and. i>=100) doexit = .true.
                   !! Check corrlen seudocode
                   !x = alms(i-100:i-50,:,:)
                   !y = alms(i-50:i,:,:)
                   !xmean = mean(alms(i-100:i-50,:,:))
                   !ymean = mean(alms(i-50:i,:,:))
                   !do l = 1, 50
                   !   cov(l, :) = ((x(l,:,:) - xmean)*(y(l,:,:) - ymean))/50.d0
                   !   sigx(l, :) = sqrt((x(l,:,:) - xmean)**2)/50.d0
                   !   sigy(l, :) = sqrt((y(l,:,:) - ymean)**2)/50.d0
                   !   corr_coeff = cov(l,:)/(sigx(l,:)*sigy(l,:))
                   !   if (corr_coeff >= 0.d0) then
                   !      corrlen = l
                   !   end if
                   !end do
                      
                end if
               
                
                ! Output samples, chisq, alms to file
                if (info%myid == 0) write(69, *) i, chisq(i), alms(i,:,:)
                
                call mpi_bcast(doexit, 1, MPI_LOGICAL, 0, c%comm, ierr)
                if (doexit) exit
             end do
             if (info%myid == 0) close(58)
             !deallocate(sigma_prior)
             
             ! Calculate cholesky
             ! Output like nprocs*nalm:+k, Does not work. Need to be saved in order to be distributed correctly
             if (info%myid == 0 .and. iter == 1) then! At this point, id=0 should have all values
                allocate(alms_covmat(0:nalm_tot-1, 0:nalm_tot-1, c%theta(j)%p%info%nmaps))
                do p = 1, c%theta(j)%p%info%nmaps
                   call compute_covariance_matrix(alms(nsamp/2:nsamp,0:nalm_tot-1,p), alms_covmat(0:nalm_tot-1,0:nalm_tot-1,p), .true.)
                   !call get_eigenvalues(alms_covmat(:,:,p), rgs)
                   !write(*,*) 'W = ', rgs
                end do
                open(58, file=trim(cpar%outdir)//'/alm_cholesky.dat', form='unformatted')
                write(58) steplen
                write(58) alms_covmat
                close(58)
                close(69)          
                deallocate(alms_covmat)
                if (nalm_tot > 1) deallocate(sigma_priors)
             end if
             deallocate(alms, rgs, chisq, L)
          end do
       end select

       ! Loop to next component
       c => c%next()
    end do
    deallocate(status_fit)

  end subroutine sample_specind_alm

  subroutine sample_specind_local(cpar, iter, handle)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, nalm_tot, ind
    real(dp)     :: t1, t2, ts, steplen, dalm, fwhm_prior, sigma_prior
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext, jtext
    logical :: accepted, exist, doexit, skip
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_comp),    pointer :: c    => null()
    real(dp),          allocatable, dimension(:,:,:)   :: alms, alms_covmat, L
    real(dp),          allocatable, dimension(:,:) :: m
    real(dp),          allocatable, dimension(:) :: buffer, rgs, sigma_priors, chisq

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
    
    ! Sample spectral parameters for each signal component
    allocate(status_fit(numband))
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%next()
          cycle
       end if
       if (all(c%p_gauss(2,:) == 0.d0)) then
          c => c%next()
          cycle
       end if

       ! Only sample components with lmax < 0 with local sampler; others are done with the alm sampler
       skip = .false.
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_ind >= 0) skip = .true.
       end select
       if (skip) then
          c => c%next()
          cycle
       end if


       do j = 1, c%npar

          if (c%p_gauss(2,j) == 0.d0) cycle

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
             if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(j))
          class is (comm_diffuse_comp)
             if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(j))
             call update_status(status, "nonlin start " // trim(c%label)// ' ' // trim(c%indlabel(j)))

             ! Set up type of smoothing scale
             id_native    = 0

             ! Compute smoothed residuals
             nullify(info)
             status_amp   = 0
             status_fit   = 0
             smooth_scale = c%smooth_scale(j)
             do i = 1, numband
                if (cpar%num_smooth_scales == 0) then
                   status_fit(i)   = 1    ! Native
                else
                   if (.not. associated(data(i)%N_smooth(smooth_scale)%p) .or. &
                        & data(i)%bp(0)%p%nu_c < c%nu_min_ind(j) .or. &
                        & data(i)%bp(0)%p%nu_c > c%nu_max_ind(j)) then
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
                   info  => comm_mapinfo(data(i)%res%info%comm, cpar%nside_smooth(j), cpar%lmax_smooth(j), &
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
                info  => comm_mapinfo(c%x%info%comm, cpar%nside_smooth(j), cpar%lmax_smooth(j), &
                     & c%x%info%nmaps, c%x%info%pol)
                call smooth_map(info, .true., &
                     & data(1)%B_smooth(smooth_scale)%p%b_l*0.d0+1.d0, c%x, &  
                     & data(1)%B_smooth(smooth_scale)%p%b_l,           c%x_smooth)
             end if

             ! Compute smoothed spectral index maps
             allocate(c%theta_smooth(c%npar))
             do k = 1, c%npar
                if (k == j) cycle
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

          end select
          
          ! Sample spectral parameters
          call c%sampleSpecInd(handle, j)
          
          ! Clean up temporary data structures
          select type (c)
          class is (comm_line_comp)
          class is (comm_diffuse_comp)
             
             if (associated(c%x_smooth)) then
                call c%x_smooth%dealloc()
                nullify(c%x_smooth)
             end if
             do k =1, c%npar
                if (k == j) cycle
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

             smooth_scale = c%smooth_scale(j)
             if (cpar%num_smooth_scales > 0) then
                if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                   ! Smooth index map with a postprocessing beam
                   !deallocate(c%theta_smooth)
                   allocate(c%theta_smooth(c%npar))
                   info  => comm_mapinfo(c%theta(j)%p%info%comm, cpar%nside_smooth(smooth_scale), &
                        & cpar%lmax_smooth(smooth_scale), c%theta(j)%p%info%nmaps, c%theta(j)%p%info%pol)
                   call smooth_map(info, .false., &
                        & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, c%theta(j)%p, &  
                        & data(1)%B_postproc(smooth_scale)%p%b_l,           c%theta_smooth(j)%p)
                   c%theta(j)%p%map = c%theta_smooth(j)%p%map
                   call c%theta_smooth(j)%p%dealloc()
                   deallocate(c%theta_smooth)
                end if
             end if

             call update_status(status, "nonlin stop " // trim(c%label)// ' ' // trim(c%indlabel(j)))

          end select

          ! Subtract updated component from residual
          if (trim(c%class) /= 'ptsrc') then
             do i = 1, numband
                allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
                m               = c%getBand(i)
                data(i)%res%map = data(i)%res%map - m
                deallocate(m)
             end do
          end if

       end do

       ! Loop to next component
       c => c%next()
    end do
    deallocate(status_fit)
    
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


end module comm_nonlin_mod
