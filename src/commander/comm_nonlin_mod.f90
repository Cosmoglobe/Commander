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

    integer(i4b) :: i, j, k, r, q, p, pl, np, nlm, l_, m_, idx, delta, corrlen_init, burnin, cholesky_calc
    integer(i4b) :: nsamp, out_every, check_every, num_accepted, smooth_scale, id_native, ierr, ind
    real(dp)     :: t1, t2, ts, dalm, thresh, steplen
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior, alms_mean, alms_var, chisq_jeffreys
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext, jtext
    character(len=3) :: tag
    character(len=9) :: ar_tag
    character(len=512) :: filename

    logical :: accepted, exist, doexit
    class(comm_mapinfo), pointer :: info => null()
    class(comm_comp),    pointer :: c    => null()
    type(map_ptr),     allocatable, dimension(:) :: df

    real(dp),          allocatable, dimension(:,:,:)  :: alms
    real(dp),          allocatable, dimension(:,:)    :: m
    real(dp),          allocatable, dimension(:)      :: buffer, rgs, chisq, N, C_
    integer(c_int),    allocatable, dimension(:)      :: maxit

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
             !write(*,*) "L ", c%L(:,:,1,j)
             if (c%p_gauss(2,j) == 0.d0 .or. c%lmax_ind < 0) cycle
             
             ! Set up smoothed data
             if (cpar%myid_chain == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(j))
             call update_status(status, "spec_alm start " // trim(c%label)// ' ' // trim(c%indlabel(j)))
             
             if (c%apply_jeffreys) then
                allocate(df(numband))
                do k = 1, numband
                   df(k)%p => comm_map(data(k)%info)
                end do
             end if
             

             call wall_time(t1)

             info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
                  & c%x%info%lmax, c%x%info%nmaps, c%x%info%pol)
             
             ! Params
             write(jtext, fmt = '(I1)') j ! Create j string
             out_every = 10
             check_every = 100
             nsamp = 500
             burnin = 3 ! Gibbs iter burnin. Tunes steplen.
             cholesky_calc = 1 ! Which gibbs iter to calculate cholesky, then corrlen.

             thresh = FLOAT(check_every)*0.8d0 !40.d0 ! 40.d0

             corrlen_init = 1
             if (info%myid == 0 .and. c%L_read(j)) then
                write(*,*) "Sampling with cholesky matrix"
             end if

             if (info%myid == 0 .and. maxval(c%corrlen(j,:)) > 0) nsamp = maxval(c%corrlen(j,:))
             call mpi_bcast(nsamp, 1, MPI_INTEGER, 0, c%comm, ierr)

             ! Static variables
             doexit = .false.
             
             allocate(chisq(0:nsamp))
             allocate(alms(0:nsamp, 0:c%nalm_tot-1,info%nmaps))                         
             allocate(rgs(0:c%nalm_tot-1)) ! Allocate random vector
             allocate(maxit(info%nmaps)) ! maximum iteration 
             maxit = 0

             if (info%myid == 0) open(69, file=trim(cpar%outdir)//'/nonlin-samples_'//trim(c%label)//'_par'//trim(jtext)//'.dat', status = 'unknown', access = 'append', recl=10000)

             ! Save initial alm        
             alms = 0.d0
             ! Gather alms from threads to alms array with correct indices
             do pl = 1, c%theta(j)%p%info%nmaps
                call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, 0, pl, pl)
                allocate(buffer(c%nalm_tot))
                call mpi_allreduce(alms(0,:,pl), buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                alms(0,:,pl) = buffer
                deallocate(buffer)
             end do
             
             do pl = 1, c%theta(j)%p%info%nmaps
                ! if sample only pol, skip T
                if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                
                ! p already calculated if larger than poltype 
                if (pl > c%poltype(j)) cycle

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


                

                ! Calculate initial chisq
                !if (info%myid == 0) open(54,file='P_jeffreys.dat')
                !do p = 1, 100

                !   if (c%theta(j)%p%info%nalm > 0) c%theta(j)%p%alm = (-4.d0 + 0.02d0*p)*sqrt(4.d0*pi)

                   if (allocated(c%indmask)) then
                      call compute_chisq(c%comm, chisq_fullsky=chisq(0), mask=c%indmask)
                   else
                      call compute_chisq(c%comm, chisq_fullsky=chisq(0))
                   end if
                   if (c%apply_jeffreys) then
                      call c%updateMixmat(df=df, par=j)
                      call compute_jeffreys_prior(c, df, pl, j, chisq_jeffreys)
                   end if

                !   if (info%myid == 0) write(54,*) -4.d0 + 0.02d0*p, chisq_jeffreys, chisq(0)
                !   if (info%myid == 0) write(*,*) -4.d0 + 0.02d0*p, chisq_jeffreys, chisq(0)
                   if (c%apply_jeffreys) chisq(0) = chisq(0) + chisq_jeffreys

                !end do
                !if (info%myid == 0) close(54)
                !call mpi_finalize(p)
                !stop
             

                call wall_time(t1)
                if (info%myid == 0) then 
                   ! Add prior 
                   chisq_prior = ((alms(0,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                   if (c%nalm_tot > 1) then
                      do p = 1, c%nalm_tot-1
                         chisq_prior = chisq_prior + (alms(0,p,pl)/c%sigma_priors(p,j))**2
                      end do
                   end if
                   chisq(0) = chisq(0) + chisq_prior
                
                   ! Output init sample
                   write(*,fmt='(a, i6, a, f16.2, a, 3f7.2)') "# sample: ", 0, " - chisq: " , chisq(0), " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
                end if


                ! Sample new alms (Account for poltype)
                num_accepted = 0
                chisq(1:) = 0.d0
                do i = 1, nsamp                   
                   
                   ! Gather alms from threads to alms array with correct indices
                   call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                   
                   ! Send all alms to 0 (Dont allreduce because only root will do calculation)
                   allocate(buffer(c%nalm_tot))
                   call mpi_reduce(alms(i,:,pl), buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)
                   alms(i,:,pl) = buffer
                   deallocate(buffer)
                   
                   ! Propose new alms
                   if (info%myid == 0) then
                      do p = 0, c%nalm_tot-1
                         rgs(p) = c%steplen(pl,j)*rand_gauss(handle)     
                      end do                                               
                      alms(i,:,pl) = alms(i-1,:,pl) + matmul(c%L(:,:,pl,j), rgs)

                      !write(*,*) "steplen", c%steplen(pl,j), "cholesky", c%L(:,:,pl,j)
                      !write(*,*) "alms", alms(i-1,:,pl)

                      ! Adding prior
                      chisq_prior = ((alms(i,0,pl) - sqrt(4*PI)*c%p_gauss(1,j))/c%p_gauss(2,j))**2
                      if (c%nalm_tot > 1) then
                         do p = 1, c%nalm_tot-1
                            !write(*,*) "alms ", p, alms(i,p,pl), c%sigma_priors(p,j)
                            chisq_prior = chisq_prior + (alms(i,p,pl)/c%sigma_priors(p,j))**2
                         end do
                      end if
                      
                      ! Prior adjustments (Nessecary because of loop adjustment)
!                      if (c%poltype(j) == 1) chisq_prior = chisq_prior*c%theta(j)%p%info%nmaps ! IF sampling pol
!                      if (c%poltype(j) == 2 .and. pl == 2) chisq_prior = chisq_prior!*2.d0
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
                      call compute_chisq(c%comm, chisq_fullsky=chisq(i), mask=c%indmask)
                   else
                      call compute_chisq(c%comm, chisq_fullsky=chisq(i))
                   end if
                   if (c%apply_jeffreys) chisq(i) = chisq(i) + chisq_jeffreys

                   ! Accept/reject test
                   ! Reset accepted bool
                   accepted = .false.
                   if (info%myid == 0) then

                      chisq(i) = chisq(i) + chisq_prior

                      diff = chisq(i-1)-chisq(i)
                      if ( chisq(i) > chisq(i-1) ) then             
                         ! Small chance of accepting this too
                         ! Avoid getting stuck in local mminimum
                         accepted = (rand_uni(handle) < exp(0.5d0*diff))
                      else
                         accepted = .true.
                      end if

                      ! Count accepted and assign chisq values
                      if (accepted) then
                         num_accepted = num_accepted + 1
                         ar_tag = " accepted"
                      else
                         chisq(i) = chisq(i-1)
                         ar_tag = " rejected"
                      end if

                      write(*,fmt='(a,i6, a, f14.2, a, f10.2, a, f7.4, a)') tag, i, " - chisq: " , chisq(i), " - diff: ", diff, " - a00-prop: ", alms(i,0,pl)/sqrt(4.d0*PI), ar_tag
                   end if

                   ! Broadcast result of accept/reject test
                   call mpi_bcast(accepted, 1, MPI_LOGICAL, 0, c%comm, ierr)

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

                      ! Write to screen every out_every'th
                      if (mod(i,out_every) == 0) then
                         diff = chisq(i-out_every) - chisq(i) ! Output diff
                         write(*,fmt='(a,i6, a, f14.2, a, f10.2, a, f7.4)') tag, i, " - chisq: " , chisq(i), " - diff: ", diff, " - a00-curr: ", alms(i,0,pl)/sqrt(4.d0*PI)
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
                         write(*, fmt='(a, i6, a, i6, a, f8.2, a, f5.3, a, f6.2)') tag, i, " - diff last ", check_every, " ", diff, " - accept rate: ", accept_rate, " - time/sample: ", ts
                         call wall_time(t1)

                         ! Adjust steplen in tuning iteration
                         if (iter < burnin) then !( .not. c%L_read(j) .and. iter == 1) then ! Only adjust if tuning

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
                      
                      exit
                   end if
                end do
             end do

             if (info%myid == 0) close(58)
             
             ! Calculate correlation length and cholesky matrix 
             ! (Only if first iteration and not initialized from previous)
             if (info%myid == 0 .and. maxval(c%corrlen(j,:)) == 0) then
                if (c%L_read(j)) then
                   write(*,*) 'Calculating correlation function'
                   ! Calculate Correlation length
                   delta = 100
                   allocate(C_(delta))
                   allocate(N(delta))
                   open(58, file=trim(cpar%outdir)//'/correlation_function_'//trim(c%label)//'_par'//trim(jtext)//'.dat', recl=10000)
                   do pl = 1, c%theta(j)%p%info%nmaps

                      ! Skip signals with poltype tag
                      if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                      if (pl > c%poltype(j)) cycle

                      ! Calculate correlation function per alm
                      do p = 0, c%nalm_tot-1
                         N(:) = 0
                         C_(:) = 0.d0
                         alms_mean = mean(alms(:,p,pl))
                         alms_var = variance(alms(:,p,pl))
                         do q = 1, maxit(pl)
                            do k = 1, delta
                               if (q+k > maxit(pl)) cycle
                               C_(k) = C_(k) + (alms(q,p,pl)-alms_mean)*(alms(q+k,p,pl)-alms_mean)
                               N(k) = N(k) + 1 ! Less samples every q
                            end do
                         end do

                         where (N>0) C_ = C_/N
                         if ( alms_var > 0 ) C_ = C_/alms_var

                         write(58,*) p, C_ ! Write to file

                         ! Find correlation length
                         c%corrlen(j,pl) = corrlen_init
                         do k = corrlen_init, delta
                            if (C_(k) > 0.1) c%corrlen(j,pl) = k
                         end do

                      end do
                      write(*,*) "Correlation length (< 0.1): ", c%corrlen(j,pl) 
                   end do
                   close(58)
                   deallocate(C_, N)

                   ! If both corrlen and L have been calulated then output
                   filename = trim(cpar%outdir)//'/init_alm_cholesky_'//trim(c%label)//'_par'//trim(jtext)//'.dat'
                   open(58, file=filename, recl=10000)
                   write(58,*) c%corrlen(j,:)
                   write(58,*) c%L(:,:,:,j)
                   close(58)
                else if (iter == cholesky_calc) then !if (.false.) then 
                   ! If L does not exist yet, calculate
                   write(*,*) 'Calculating cholesky matrix'
                   do pl = 1, c%theta(j)%p%info%nmaps
                      if (maxit(pl) == 0) cycle ! Cycle if not sampled
                      call compute_covariance_matrix(alms(INT(maxit(pl)/2):maxit(pl),0:c%nalm_tot-1,pl), c%L(0:c%nalm_tot-1,0:c%nalm_tot-1,pl,j), .true.)
                   end do
                   c%steplen(:,j) = 0.3d0
                   c%L_read(j) = .true. ! L now exists!
                end if

             end if

             if (info%myid == 0) close(69)   

             deallocate(alms, rgs, chisq, maxit)

             ! Clean up
             if (c%apply_jeffreys) then
                do k = 1, numband
                   call df(k)%p%dealloc()
                end do
                deallocate(df)
             end if

          end do ! End of j
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
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind
    real(dp)     :: t1, t2, ts, dalm, fwhm_prior
    real(dp)     :: mu, sigma, par, accept_rate, diff, chisq_prior
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext, jtext
    logical :: accepted, exist, doexit, skip
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_comp),    pointer :: c    => null()
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
          ELSE
             do i = 1, numband
                call int2string(i, itext)
                call data(i)%res%writeFITS("res"//itext//".fits")
             end do
!!$             call mpi_finalize(ierr)
!!$             stop
          end if

          ! Set up smoothed data
          select type (c)
          class is (comm_line_comp)
             if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label), ' ', trim(c%indlabel(j))
          class is (comm_ptsrc_comp)
             if (cpar%myid == 0) write(*,*) '   Sampling ', trim(c%label)
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
