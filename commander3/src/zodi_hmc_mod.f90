module zodi_hmc_mod
   use comm_zodi_mod
   use comm_chisq_mod
   use comm_zodi_samp_mod
   use hmc_mod
   implicit none

   private
   public minimize_zodi

 contains

   subroutine minimize_zodi(cpar, iter, handle, samp_group)
      implicit none 
      type(comm_params), intent(in)      :: cpar
      integer(i4b),      intent(in)      :: iter
      type(planck_rng),  intent(inout)   :: handle
      integer(i4b),      intent(in)      :: samp_group

      logical(lgt) :: accept
      real(dp), allocatable :: theta(:), theta_new(:), theta_old(:), scale(:)
      real(dp), allocatable, dimension(:) :: theta_prev, chisq_prev
      integer(i4b) :: i, j, k, ind, ierr, flag, ntot, npar, unit
      real(dp) :: chisq_old, chisq_new
      character(len=6) :: iter_string
      character(len=6) :: sgroup
      character(len=512) :: filename

      if (cpar%myid == 0) print *, "minimizing zodi parameters with", trim(cpar%zs_samp_method(samp_group)), ", samp_group =", samp_group

      npar = count(zodi_model%theta_stat(:,samp_group)==0)
      ntot = zodi_model%npar_tot
      allocate(theta_old(npar), theta_new(npar), theta(npar))
      allocate(theta_prev(npar), chisq_prev(numband), scale(npar))

      ! Initialize active monopoles
      do i = 1, numband
         ind = zodi_model%get_par_ind(mono_band=i)
         band_monopole(i) = get_monopole_amp(data(i)%label)
      end do

      if (cpar%myid_chain == 0) then
         call int2string(samp_group, sgroup)
         call int2string(iter, iter_string)
         filename = trim(cpar%outdir)//'/zodi_powell_sg'//sgroup//'_k'//iter_string//'.dat'
         unit     = getlun()
         open(unit, file=trim(filename), recl=10000)
         write(unit, '(a)', advance="no") "# chisq_red "
         do i = 1, zodi_model%npar_tot
            if (zodi_model%theta_stat(i,samp_group)==0) then
               write(unit, "(a,a)", advance="no") trim(adjustl(zodi_model%par_labels_full(i))), " "
            end if
         end do
         write(unit,*)
      end if

      ! Get chisq of old point
      scale = pack(zodi_model%theta_scale(:,1), zodi_model%theta_stat(:,samp_group)==0)
      call model_to_params(zodi_model, theta_old, samp_group)
   !!$      if (cpar%myid == cpar%root) then
   !!$         do i = 1, npar
   !!$            write(*,*) i, theta_old(i)
   !!$         end do
   !!$      end if
   !!$      call mpi_finalize(ierr)
   !!$      stop

      ! Enforce priors; rms = 0
      call randomize_zodi_init(theta_old, samp_group, cpar, handle, rms=0.d0)

      theta_prev = 0.d0
      chisq_prev = 0.d0
      if (cpar%myid == cpar%root) then
         theta = theta_old/scale
         chisq_old  = lnL_zodi(theta)
      else
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
         chisq_old  = lnL_zodi()
      end if
      if (cpar%zs_output_tod_res) then
         call mpi_finalize(ierr)
         stop
      end if

      ! Initialize new point
      theta_new = theta_old
      call randomize_zodi_init(theta_new, samp_group, cpar, handle)

      ! Rescale 
      theta = theta_new/scale

      if (cpar%myid == cpar%root) then
         
         ! Perform search
         if (trim(cpar%zs_samp_method(samp_group))=='powell') then
            !call powell(theta, lnL_zodi, ierr, tolerance=1d-5)  
            write (*,*) 'yee powell'
            stop
         else if (trim(cpar%zs_samp_method(samp_group))=='hmc') then 
            !call hmc(theta, lnL_zodi, grad_lnL_zodi, ierr, handle) ! tolerance=1d-5)
            write (*,*) 'yee hmc'
            stop
         else
            write(*,*) 'Unsupported zs_samp_method=', trim(cpar%zs_samp_method(samp_group))
            stop
         end if

         chisq_new = lnL_zodi(theta)
         flag = 0
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)

         ! Apply approximate Metropolis rule, using reduced chisq instead of chisq
         if (chisq_new < chisq_old) then
            accept = .true.
         else
            accept = rand_uni(handle) < exp(-0.5d0*(chisq_new-chisq_old)/0.02d0)
         end if
         if (accept) then
            ! Accept new point; update
            if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample accepted, chisq_new =', chisq_new, ', chisq_old = ', chisq_old
            theta_new = theta*scale ! Convert to physical units
         else
            ! Reject new solution, reset to previous solution
            if (cpar%myid == cpar%root) write(*,fmt='(a,f8.2,a,f8.2)') 'Zodi sample rejected, chisq_new =', chisq_new, ', chisq_old = ', chisq_old
            theta_new = theta_old
         end if
      else
         do while (.true.)
             call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
             if (flag == 1) then 
                 chisq_new = lnL_zodi()
             else
                 exit
              end if
           end do
          !call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)          
      end if
      if (cpar%myid_chain == 0) close(unit)

      ! Distribute final solution
      call mpi_bcast(theta_new, size(theta_new), MPI_DOUBLE_PRECISION, cpar%root, cpar%comm_chain, ierr)

      ! update model with final parameters
      call params_to_model(zodi_model, theta_new, samp_group)

      ! Update monopole for requested bands
      do i = 1, numband
         if (band_update_monopole(i,samp_group)) call set_monopole_amp(data(i)%label, band_monopole(i))
      end do

      deallocate(theta_old, theta_new, theta, theta_prev, chisq_prev, scale)

    contains

      function lnL_zodi(p)
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p
         real(dp)                                     :: lnL_zodi

         real(dp), allocatable :: theta(:)
         real(dp) :: chisq, chisq_tot, box_width, t1, t2, t3, t4, mono
         integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, ndof, ndof_tot
         logical(lgt) :: accept
         logical(lgt), dimension(numband) :: update_band
         character(len=4) :: scan_str
         type(hdf_file) :: tod_file

         call wall_time(t1)   !what is this?
         
         allocate(theta(npar))
         if (cpar%myid_chain == 0) then
            flag = 1
            call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
            theta = p*scale
         end if
         call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)

         ! Check which parameters have changed
         update_band = .true.
   
         ! Check priors
         if (cpar%myid_chain == 0) then
            chisq_tot = get_chisq_priors(theta, samp_group)
            accept    = chisq_tot < 1.d30
         else
            chisq_tot = 0.d0
         end if
         call mpi_bcast(accept, 1, MPI_LOGICAL, 0, data(1)%tod%comm, ierr)

         if (.not. accept) then
            deallocate(theta)
            lnL_zodi = 1.d30
            return
         end if

         !overwrite to the background table
         call params_to_model(zodi_model, theta, samp_group)
         if (cpar%myid_chain == 0) call print_zodi_model(theta, samp_group)

         ndof = 0
         do i = 1, numband
            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle
            if (.not. zodi_model%sampgroup_active_band(i,samp_group)) cycle
            ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
            if (chisq_tot >= 1.d30) exit

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan

            if (.not. update_band(i)) then
               if (cpar%myid_chain == 0) write(*,*) 'skipping band ', i, chisq_prev(i)
               chisq_tot = chisq_tot + chisq_prev(i)
               do scan = 1, nscan
                  do j = 1, ndet
                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle
                     ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  end do
               end do
               cycle
            end if

            ! Get monopole
            mono = band_monopole(i) !get_monopole_amp(data(i)%label)

            box_width = get_boxwidth(data(i)%tod%samprate_lowres, data(i)%tod%samprate)

            ! Make sure that the zodi cache is cleared before each new band
            call data(i)%tod%clear_zodi_cache()

            ! Evaluate zodi model with newly proposed values for each band and calculate chisq
            chisq_prev(i) = 0.d0
            do scan = 1, nscan
               ! Skip scan if no accepted data
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) then
                     cycle
                  end if

                  !difference between get_zodi_emission and get_s_zodi?
                  !I will need to evaluate zodi with only n0 as a free parameter, I don't know which function to use
                  call wall_time(t3)
                  call get_zodi_emission(&
                      & tod=data(i)%tod, &
                      & pix=data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                      & scan=scan, &
                      & det=j, &
                      & s_zodi_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &  !by changeing these arguments can I get zodi(n0; k98)?
                      & s_zodi_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                      & model=zodi_model, &
                      & use_lowres_pointing=.true. &
                      &)

                  call wall_time(t4)
                  call get_s_zodi(i, &
                      & s_therm=data(i)%tod%scans(scan)%d(j)%downsamp_therm, &
                      & s_scat=data(i)%tod%scans(scan)%d(j)%downsamp_scat, &
                      & s_zodi=data(i)%tod%scans(scan)%d(j)%downsamp_zodi &
                      &)

                  call wall_time(t3)
                  chisq = sum( &    !what does this do? because I will need single components of this sum, can I just use them or do I havw to sum something?
                     & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                     &   - data(i)%tod%scans(scan)%d(j)%gain * &
                     &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                     &   + data(i)%tod%scans(scan)%d(j)%downsamp_zodi &   !this is the zodi part of the model, right?
                     &   + mono) &
                     & )/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0/sqrt(box_width)))**2 &   !this is sigma, right? 
                     &)
                  chisq_tot     = chisq_tot     + chisq
                  chisq_prev(i) = chisq_prev(i) + chisq
                  call wall_time(t4)

                  ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  if (chisq_tot >= 1.d30) exit
                  if (cpar%zs_output_tod_res) then
                     call int2string(data(i)%tod%scanid(scan), scan_str)
                     open(58,file=trim(cpar%outdir)//'/todres_'//trim(data(i)%tod%freq)//'_'//scan_str//'.dat', recl=2048)
                     do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                        write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_point(k,:), data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
                             &   - data(i)%tod%scans(scan)%d(j)%gain * &
                             &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) &
                             &   + data(i)%tod%scans(scan)%d(j)%downsamp_zodi(k) &
                             &   + mono)
                     end do
                     close(58)
                  end if
                  call wall_time(t3)
                  
               end do
            end do
         end do

         ! Reduce chisq to root process
         call mpi_reduce(chisq_tot, chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
         call mpi_reduce(ndof,      ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)
         call wall_time(t4)
         
         if (cpar%myid_chain == 0) then
            lnL_zodi = chisq/ndof_tot  !why chisq and not chisq_tot?
            call wall_time(t2)
            if (ndof_tot > 0) write(*,fmt='(a,e16.8,a,f10.4,a,f8.3)') "chisq_zodi = ", chisq, ", chisq_red = ", chisq/ndof_tot, ", time = ", t2-t1
            write(*,*)
            write(unit,*) chisq/ndof_tot, real(theta,sp)
         end if

         theta_prev = theta
      end function lnL_zodi

      !function grad_lnL_zodi(p)

      !Currently the only free parameter is n0 which is linear in the model, this means 
      !grad(lnL_zodi)=d(lnL_zodi)/dno=d(chisq)/dno=-(2/sigma)*sqrt(chisq)*zodi_model(n0=1)
      !and here lnL_zodi is the k98 model with only the density parameter n0 as something that I sample

      

      
      !end function grad_lnL_zodi

   end subroutine minimize_zodi

end module zodi_hmc_mod
