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

    integer(i4b) :: i, j, k, q, pl, np, nlm, out_every, num_accepted, smooth_scale, id_native, ierr, l_, m_, ind
    real(dp)     :: t1, t2, lr
    real(dp)     :: mu, sigma, par, chisq, chisq_old, chisq_d, accept_rate, diff
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    character(len=2) :: itext, jtext
    logical :: accepted
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_comp),    pointer :: c    => null()
    real(dp),          allocatable, dimension(:,:) :: m, alm_old, alm, lm, buffer2
    real(dp),          allocatable, dimension(:)   :: buffer
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
       
       do j = 1, c%npar

          if (c%p_gauss(2,j) == 0.d0) cycle


          ! Add current component back into residual
!!$          call int2string(j,itext)
!!$          call data(j)%res%writeFITS('res_before'//itext//'.fits')

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

!!$                call int2string(i,itext)
!!$                call res_smooth(i)%p%writeFITS('res'//itext//'.fits')

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

             !c%x_smooth%map = c%x_smooth%map * c%F(1)%p%map

!!$             call c%x_smooth%writeFITS('x.fits')
!!$             call mpi_finalize(i)
!!$             stop
             if (c%lmax_ind > 0) then
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
             end if
                
          ! Trygve spectral index sampling
          call wall_time(t1)
          if (trim(c%operation) == 'optimize') then
             ! Calls comm_diffuse_comp_mod
             call c%sampleSpecInd(handle, j)

          else if (trim(c%operation) == 'sample' .and. c%lmax_ind > 0) then
             if (info%myid == 0) open(69, file='nonlin-samples.dat', recl=10000)
             allocate(buffer(0:c%theta(j)%p%info%nalm))
             
             lr = 0.001 ! Init learning rate for proposals
             num_accepted = 0
             out_every = 10

             ! Cholesky decompose covariance to get faster sampling.
             ! call cholesky_decompose_single(s, ierr=ierr)
             ! Save all values in buffer. Only accepted values?
             ! After 1000? samples, calculate covariance
           
             ! MCMC loop
             do i = 1, 500
                call wall_time(t1)

                ! Reset accepted bool
                accepted = .false.

                ! Save old alms
                alm_old = c%theta(j)%p%alm(:,:)

                ! Calculate current chisq
                call compute_chisq(c%comm, chisq_fullsky=chisq_old)

                ! Sample new alms (Account for poltype)
                do pl = 1, c%theta(j)%p%info%nmaps

                   ! if sample only pol, skip T
                   if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 

                   ! p already calculated if larger than poltype ( smart ;) )
                   if (pl > c%poltype(j)) cycle

                   ! Propose new alms
                   do k = 0, c%theta(j)%p%info%nalm-1   
                      buffer(k) = c%theta(j)%p%alm(k,pl) + lr*rand_gauss(handle)
                   end do
                      
                   ! Save to correct poltypes
                   if (c%poltype(j) == 1) then      ! {T+E+B}
                      do q = 1, c%theta(j)%p%info%nmaps
                         c%theta(j)%p%alm(:,q)  = buffer
                      end do
                   else if (c%poltype(j) == 2) then ! {T,E+B}
                      if (pl == 1) then
                         c%theta(j)%p%alm(:,1)  = buffer
                      else
                         do q = 2, c%theta(j)%p%info%nmaps
                            c%theta(j)%p%alm(:,q)  = buffer
                         end do 
                      end if
                   else if (c%poltype(j) == 3) then ! {T,E,B}
                         c%theta(j)%p%alm(:,pl)  = buffer
                   end if
                end do

                ! Update mixing matrix with new alms
                call c%updateMixmat

                ! Calculate proposed chisq
                call compute_chisq(c%comm, chisq_fullsky=chisq)
                diff = chisq_old-chisq

                ! Accept/reject test
                if (info%myid == 0) then
                   if ( chisq > chisq_old ) then                 
                      ! Small chance of accepting this too
                      ! Avoid getting stuck in local mminimum
                      accepted = (rand_uni(handle) < exp(0.5d0*diff))
                      !write(*,*) "exp", exp(0.5d0*diff), accepted, diff
                   else
                      accepted = .true.
                   end if
                end if

                ! Broadcast result of accept/reject test
                call mpi_bcast(accepted, 1, MPI_LOGICAL, 0, c%comm, ierr)
                
                ! Count accepted, if rejected revert changes.
                if (accepted) then
                   !if (info%myid == 0) write(*,fmt='(i8, a, f10.2, a, f8.2)') i, " chisq ", chisq, " diff ", diff
                   num_accepted = num_accepted + 1
                else
                   ! If rejected, restore old values
                   !if (info%myid == 0) write(*,*) "Rejected"
                   c%theta(j)%p%alm(:,:) = alm_old
                   call c%updateMixmat                
                   chisq = chisq_old
                end if

                ! Write to screen first iter
                if (i == 1) then
                   if (info%myid == 0) then 
                      write(*,fmt='(a,i6, a, f16.2)') "- sample ", i, " -- chisq " , chisq
                      write(69, *) ' sample       chisq           alms'
                   end if
                   chisq_d = chisq
                end if

                ! Write to screen
                if (mod(i,out_every) == 0) then
                   call wall_time(t2)
                   if (info%myid == 0) write(*,fmt='(a,i6, a, f16.2, a, f10.5, a, f7.3)') "- sample ", i, " -- chisq " , chisq, " diff ", diff, " time/sample ", t2-t1 
                   chisq_d = chisq
                end if

                ! Adjust learning rate every 100th
                if (mod(i, 100) == 0) then 
                   ! Accept rate
                   accept_rate = num_accepted/100.d0
                   num_accepted = 0

                   ! Write to screen
                   if (info%myid == 0) write(*, fmt='(a, f5.3)') "- accept rate ", accept_rate
                   
                   ! Adjusrt lr
                   if (accept_rate < 0.2) then                 
                      lr = lr*0.5d0
                      if (info%myid == 0) write(*,fmt='(a,f10.5)') "Reducing lr -> ", lr
                   else if (accept_rate > 0.8) then
                      lr = lr*2.0d0
                      if (info%myid == 0) write(*,fmt='(a,f10.5)') "Increasing lr -> ", lr
                   end if
                end if

                ! Output samples, chisq, alms to file
                if (info%myid == 0) then 
                   allocate(alm(0:(c%lmax_ind+1)**2-1,info%nmaps))

                   ! Dumping whatever is in proc = 0
                   do k = 0, c%theta(j)%p%info%nalm-1
                         l_ = c%theta(j)%p%info%lm(1,k)
                         m_ = c%theta(j)%p%info%nmapslm(2,k)
                         ind = l_**2 + l_ + m_
                         alm(ind,:) = c%theta(j)%p%alm(k,:)
                   end do

                   do np = 1, info%nprocs-1
                      call mpi_recv(nlm, 1, MPI_INTEGER, np, 420, c%comm, mpistat, ierr)
                      allocate(lm(2,0:nlm-1))
                      call mpi_recv(lm, size(lm), MPI_INTEGER, i, 98, c%comm, mpistat, ierr)
                      allocate(buffer2(0:nlm, c%theta(j)%p%info%nmaps))
                      call mpi_recv(buffer2, size(buffer2), MPI_DOUBLE_PRECISION, np, 420, c%comm, mpistat, ierr)
                      do k = 0, nlm-1
                         l_ = lm(1,k)
                         m_ = lm(2,k)
                         ind = l_**2 + l_ + m_
                         alm(ind,:) = buffer2(k,:)
                      end do
                      deallocate(lm, buffer2)
                   end do
                   write(69, fmt='(i6, f16.2, 999f6.2)') i, chisq, alm
                   deallocate(alm)
                else
                   call mpi_send(c%theta(j)%p%info%nalm, 1, MPI_INTEGER, 0, 420, c%comm, mpistat, ierr)
                   call mpi_send(c%theta(j)%p%info%lm, size(c%theta(j)%p%info%lm), MPI_INTEGER, 0, 420, c%comm, mpistat, ierr)
                   call mpi_recv(c%theta(j)%p%alm, size(c%theta(j)%p%alm), MPI_DOUBLE_PRECISION, 0, 420, c%comm, mpistat, ierr)
                end if
              
             end do
             deallocate(alm_old, buffer)
          end if
          end select

          ! Clean up temporary data structures
          select type (c)
          class is (comm_line_comp)
          class is (comm_diffuse_comp)
             
             if (associated(c%x_smooth)) then
                call c%x_smooth%dealloc()
                nullify(c%x_smooth)
             end if
             do k = 1, c%npar
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

!!$          call int2string(j,itext)
!!$          call data(j)%res%writeFITS('res_after'//itext//'.fits')


       end do

       ! Loop to next component
       c => c%next()
    end do
    deallocate(status_fit)

!!$    call mpi_finalize(i)
!!$    stop

    ! Sample calibration factors
    do i = 1, numband
       if (.not. data(i)%sample_gain) cycle
       call sample_gain(cpar%operation, i, cpar%outdir, cpar%mychain, iter, handle)
    end do


    ! Update mixing matrices
    if (any(data%sample_gain)) then
       c => compList
       do while (associated(c))
          call c%updateMixmat
          c => c%next()
       end do
    end if

    call wall_time(t2)
    if (cpar%myid == 0) write(*,*) 'CPU time specind = ', real(t2-t1,sp)
    
  end subroutine sample_nonlin_params

end module comm_nonlin_mod
