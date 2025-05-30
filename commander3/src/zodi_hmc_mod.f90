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
      real(dp), allocatable :: theta(:), theta_new(:), theta_old(:), scale(:), grad_new(:), hmc_mass(:)
      real(dp), allocatable, dimension(:) :: theta_prev, chisq_prev
      integer(i4b) :: i, j, k, ind, ierr, flag, ntot, npar, unit
      real(dp) :: chisq_old, chisq_new
      character(len=6) :: iter_string
      character(len=6) :: sgroup
      character(len=512) :: filename

      if (cpar%myid == 0) write(*,*) "minimizing zodi parameters with ", trim(cpar%zs_samp_method(samp_group)), ", samp_group = ", samp_group

      npar = count(zodi_model%theta_stat(:,samp_group)==0)
      ntot = zodi_model%npar_tot
      allocate(theta_old(npar), theta_new(npar), theta(npar))
      allocate(theta_prev(npar), chisq_prev(numband), scale(npar))
      allocate(grad_new(npar))

      ! Initialize active monopoles
      do i = 1, numband
         ind = zodi_model%get_par_ind(mono_band=i)
         band_monopole(i) = get_monopole_amp(data(i)%label)
      end do

      if (cpar%myid_chain == 0) then
         call int2string(samp_group, sgroup)
         call int2string(iter, iter_string)
         filename = trim(cpar%outdir)//'/zodi_' // trim(cpar%zs_samp_method(samp_group)) // '_sg'//sgroup//'_k'//iter_string//'.dat'
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
      call zodi_model%model_to_params(theta_old, samp_group)
      
      ! Enforce priors; rms = 0
      call randomize_zodi_init(theta_new, samp_group, cpar, handle)

      theta_prev = 0.d0
      chisq_prev = 0.d0
      if (cpar%myid == cpar%root) then
         theta = theta_old/scale
         chisq_old  = lnL_zodi_hmc(theta, .true.)
      else
         call mpi_bcast(flag, 1, MPI_INTEGER, cpar%root, cpar%comm_chain, ierr)
         chisq_old  = lnL_zodi_hmc()
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
            call powell(theta, lnL_zodi_hmc, ierr, tolerance=1d-5) 
         else if (trim(cpar%zs_samp_method(samp_group))=='hmc') then
            write(*,*) "calling nuts"
            !allocate(hmc_mass(npar))
            !hmc_mass = 0.0001
            !call nuts(theta, lnL_zodi_hmc, grad_lnL_zodi_hmc, ierr, handle, M=hmc_mass) !does it need tolerance=1d-5? ierr=number of hmc iteration? 
            !deallocate(hmc_mass)
            call nuts(theta, lnL_zodi_hmc, grad_lnL_zodi_hmc, ierr, handle)
         else
            write(*,*) 'Unsupported zs_samp_method=', trim(cpar%zs_samp_method(samp_group))
            stop
         end if

         chisq_new = lnL_zodi_hmc(theta, .true.)
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
                 chisq_new = lnL_zodi_hmc()
             else if (flag == 2) then
                 call grad_lnL_zodi_hmc(grad_lnL_zodi=grad_new)
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
      call zodi_model%params_to_model(theta_new, samp_group)

      ! Update monopole for requested bands
      do i = 1, numband
         if (band_update_monopole(i,samp_group)) call set_monopole_amp(data(i)%label, band_monopole(i))
      end do

      deallocate(theta_old, theta_new, theta, theta_prev, chisq_prev, scale, grad_new)

    contains

    function lnL_zodi_hmc(p, checks)
      !same as lnL_zodi
         use healpix_types
         implicit none
         real(dp), dimension(:), intent(in), optional :: p !optional so that can be called on threads not root that just take theta from the root thread
         logical(lgt), intent(in), optional  :: checks !when true means lnL is being used to actually produce a part of the analysis so the results shouls make sense, 
                                                       !otherwise the function is being used to tune hmc so it is not necessay to check the results.

         real(dp)                                     :: lnL_zodi_hmc

         real(dp), allocatable :: theta(:)
         real(dp) :: chisq, chisq_tot, box_width, t1, t2, t3, t4, mono
         integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, ndof, ndof_tot
         logical(lgt) :: accept, flag_checks
         !logical(lgt), dimension(numband) :: update_band
         character(len=4) :: scan_str
         type(hdf_file) :: tod_file
         real(sp), allocatable, dimension(:) :: s_zodi

         real(dp) :: norm !a constant used to make the target distribution more similar to the actual log_like instead of the chisquare

         call wall_time(t1)

         allocate(theta(npar))
         
         if (cpar%myid_chain == 0) then
            if(present(checks)) then
               flag_checks = checks 
            else
               flag_checks = .true. !easyer to specify when it needs to be false 
            end if

            flag = 1
            call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)

            theta = p*scale
            
         end if

         call mpi_bcast(flag_checks, 1, MPI_LOGICAL, 0, data(1)%tod%comm, ierr)
         call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)

         if(flag_checks) then
            !Check priors
            if (cpar%myid_chain == 0) then
               chisq_tot = get_chisq_priors(theta, samp_group)
               accept    = chisq_tot < 1.d30
            else
               chisq_tot = 0.d0
            end if
            call mpi_bcast(accept, 1, MPI_LOGICAL, 0, data(1)%tod%comm, ierr)

            if (.not. accept) then
               deallocate(theta)
               lnL_zodi_hmc = 1.d30
               return
            end if
         end if

         call zodi_model%params_to_model(theta, samp_group) 

         if ((flag_checks) .and. (cpar%myid_chain == 0)) call print_zodi_model(theta, samp_group)

         ndof = 0
         do i = 1, numband
            if (data(i)%tod_type == "none") cycle
            if (.not. data(i)%tod%sample_zodi) cycle 
            if (.not. zodi_model%sampgroup_active_band(i,samp_group)) cycle 
            if (flag_checks) then
            ! If chisq is already too large, skip rest of the evaluation and go directly to rejection
               if (chisq_tot >= 1.d30) exit
            end if

            ndet = data(i)%tod%ndet
            nscan = data(i)%tod%nscan

            ! Get monopole
            mono = band_monopole(i) !get_monopole_amp(data(i)%label)
            
            do scan = 1, nscan
               ! Skip scan if no accepted data
               do j = 1, ndet
                  if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle 

                  allocate(s_zodi(size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)))

                  call wall_time(t3)
                  call get_zodi_emission(data(i)%tod, data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                     & scan, j, zodi_model, s_zodi, use_lowres_pointing=.true.)
                  call wall_time(t4)
                  !if (data(1)%tod%myid == 10) write(*,*) ' CPU2 = ', t3-t4

                  !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'tod  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_tod))
                  !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'sky  ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_sky))
                  !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'zodi ', minval((data(i)%tod%scans(scan)%d(j)%downsamp_tod)), maxval((data(i)%tod%scans(scan)%d(j)%downsamp_zodi))
                  !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'mono ', mono
                  !write(*,*) data(i)%tod%scanid(scan), data(1)%tod%myid, 'noise', data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
                  
                  chisq = sum( &
                    & ((data(i)%tod%scans(scan)%d(j)%downsamp_tod - data(i)%tod%scans(scan)%d(j)%gain * &
                    & (data(i)%tod%scans(scan)%d(j)%downsamp_sky + s_zodi + mono))/data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2)
                  chisq_tot = chisq_tot + chisq

                  call wall_time(t4)
                  !if (data(1)%tod%myid == 10) write(*,*) ' CPU3 = ', t4-t3

!!$                  write(*,*) 'a', i, scan!, allocated(data(i)%tod%scans(scan)%d(j)%downsa1mp_tod)
!!$                  write(*,*) 'b', j!, allocated(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  !write(*,*) 'do not remove -- memory corruption bug "fix"', ndof!, allocated(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  if ((flag_checks) .and. (chisq_tot >= 1.d30)) exit
                  ! call int2string(data(i)%tod%scanid(scan), scan_str)
                  ! call open_hdf_file(trim(adjustl("/mn/stornext/u3/metins/dirbe/chains/chains_downsamp/dtodlnl_"//scan_str//".h5")), tod_file, 'w')
                  ! call write_hdf(tod_file, '/dtod', data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  ! call write_hdf(tod_file, '/dzodi', data(i)%tod%scans(scan)%d(j)%downsamp_zodi)
                  ! call write_hdf(tod_file, '/dsky', data(i)%tod%scans(scan)%d(j)%downsamp_sky)
                  ! call write_hdf(tod_file, '/dpix', data(i)%tod%scans(scan)%d(j)%downsamp_pix)
                  ! call close_hdf_file(tod_file)

                  !if (.false. .and. data(1)%tod%myid == 0 .and. scan == 1) then
                  if (flag_checks) then
                     if (cpar%zs_output_tod_res) then
                        call int2string(data(i)%tod%scanid(scan), scan_str)
                        !write(*,*) "scan = ", data(i)%tod%scanid(scan), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_tod)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_sky)), sum(abs(data(i)%tod%scans(scan)%d(j)%downsamp_zodi)), data(i)%tod%scans(scan)%d(j)%N_psd%sigma0
                        open(58,file=trim(cpar%outdir)//'/todres_'//trim(data(i)%tod%freq)//'_'//scan_str//'.dat', recl=2048)
                        do k = 1, size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                           write(58,*) data(i)%tod%scans(scan)%d(j)%downsamp_point(k,:), data(i)%tod%scans(scan)%d(j)%downsamp_tod(k) &
                                &   - data(i)%tod%scans(scan)%d(j)%gain * &
                                &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky(k) &
                                &   + s_zodi(k) &
                                &   + mono)
                        end do
                        close(58)
                     end if
                  end if 
                  call wall_time(t3)
                  !if (data(1)%tod%myid == 10) write(*,*) ' CPU4 = ', t3-t4

                  deallocate(s_zodi)
               end do
            end do
         end do

         !write(*,*) data(1)%tod%myid, ' -- recomputed =', chisq_prev

         ! Reduce chisq to root process
         call mpi_reduce(chisq_tot, chisq,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
         call mpi_reduce(ndof,      ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)

         norm = 1.0_dp/ndof_tot
         !norm = 1.0_dp

         call wall_time(t4)
         !if (data(1)%tod%myid == 0) write(*,*) ' CPU5 = ', t4-t3

         !if (data(1)%tod%myid == 0) write(*,*) chisq, "|", ndof_tot

         !this too is useful only after tuning
         if (cpar%myid_chain == 0) then
            lnL_zodi_hmc = chisq*norm
            call wall_time(t2)
            if (ndof_tot > 0) write(*,fmt='(a,e16.8,a,i5,a,f10.4,a,f8.3)') "chisq_zodi = ", chisq, ", ndof_tot = ", ndof_tot, ", chisq_red = ", chisq/ndof_tot, ", time = ", t2-t1
            !write(*,*) 
            if (flag_checks) write(unit,*) chisq*norm, real(theta,sp)
         end if

         theta_prev = theta
      end function lnL_zodi_hmc

     subroutine grad_lnL_zodi_hmc(p, grad_lnL_zodi)
     !Currently the only free parameter is comp:n0 which is linear in the model, this means 
     !grad(lnL_zodi)=d(lnL_zodi)/dno=d(chisq)/dno=-(2/(dof*sigma**2))*(data-models)*zodi_model_comp(n0=1)
     !attention: the sum is outside everything
     !first I use the usual lnL_zodi, then I force comp:n0=1 setting its 'stat'=-3 and other n0s=0 with stat=-2 then I re-evaluate zodi
     
        use healpix_types
        implicit none
        real(dp), dimension(:), intent(in), optional       :: p
        real(dp), dimension(:), intent(inout) :: grad_lnL_zodi

         real(dp), dimension(npar) :: theta
         real(dp), allocatable :: zodi_mod(:), zodi_n0(:)
         real(dp) :: grad, grad_tot
         character(len=:), allocatable :: free_comp, free_par 
         integer(i4b), dimension(6) :: par_inds, prev_vals   !number of components =6 (cloud, band1, band2, band3, ring, feature)
         integer(i4b) :: g, pos, c

         real(dp) :: box_width, t1, t2, t3, t4, mono
         integer(i4b) :: i, j, k, scan, ntod, ndet, nscan, flag, ierr, ndof, ndof_tot
         logical(lgt) :: accept
         character(len=4) :: scan_str
         type(hdf_file) :: tod_file
         real(sp), allocatable, dimension(:) :: s_zodi

         real(dp) :: norm !a constant used to make the target distribution more similar to the actual log_like instead of the chisquare

         call wall_time(t1)  

         !write(*,*) 'h1'
         
         !#grad_lnL is called only when performing the minimization which is something that only 
         !#the root process performs but then all the threads are needed when evaluating lnL or grad because I need all the data
         if (cpar%myid_chain == 0) then
            flag = 2
            call mpi_bcast(flag, 1, MPI_INTEGER, 0, data(1)%tod%comm, ierr)
            theta = p*scale
         end if
         call mpi_bcast(theta, size(theta), MPI_DOUBLE_PRECISION, 0, data(1)%tod%comm, ierr)

         call params_to_model(zodi_model, theta, samp_group) !check this

         !write(*,*) 'h2'

         do g=1, size(theta)   !loop over all free parameters i.e. all components of the gradient
            ndof = 0
            grad_tot = 0 !check this
            do i = 1, numband
               ndet = data(i)%tod%ndet
               nscan = data(i)%tod%nscan

               ! Get monopole
               mono = band_monopole(i) !get_monopole_amp(data(i)%label)

               !Evaluate zodi model with newly proposed values for each band and calculate the gradient
               do scan = 1, nscan
                  ! Skip scan if no accepted data
                  do j = 1, ndet
                     if (.not. data(i)%tod%scans(scan)%d(j)%accept) cycle

                     !Actual zodi model value
                     allocate(s_zodi(size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)))
               
                     call wall_time(t3)
                     call get_zodi_emission(data(i)%tod, data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                          & scan, j, zodi_model, s_zodi, use_lowres_pointing=.true.)
                     call wall_time(t4)
                     zodi_mod = s_zodi

                     !write(*,*) 'zodi_mod'

                     !now I want zodi_n0 i.e., zodi model with n0=1

                     !look for the free parameter index in theta_stat (probably there's a more intelligent way to do this)
                     !cpar%zs_samp_groups contains all free parameters in the format component:parameter
                     pos = index(cpar%zs_samp_groups(g), ":")
                     if(pos<=0) write(*,*) "Error: cpar%zs_samp_groups(1) is not in the form comp:param"
                     allocate(character(len=pos-1) :: free_comp)
                     allocate(character(len=len(cpar%zs_samp_groups(1))-pos) :: free_par)
                     free_comp = cpar%zs_samp_groups(1)(1:pos-1)
                     free_par = cpar%zs_samp_groups(1)(pos+1:len(cpar%zs_samp_groups(1)))
                     par_inds(1) = zodi_model%get_par_ind(comp_str = free_comp, param = free_par)
                     !cpar%zs_samp_groups stores all the string that identify a free parameter, the idea will be to put a for loop on those to compute the gradient for each of them 
                     !set n0=1 i.e., theta_stat=-3  
                     if(zodi_model%theta_stat(par_inds(1),samp_group)/=0) write(*,*) "Wrong index extracted in grad_lnL"
                     prev_vals(1) = 0
                     zodi_model%theta_stat(par_inds(1),samp_group) = -3

                     !all the other n0 must be set to zero because the derivative acts only on one component
                     select case (free_comp)
                     case ('cloud')
                        !write(*,*) 'calling case cloud'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'band1', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'band2', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'band3', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'ring', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'feature', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2

                     case ('band1')
                        !write(*,*) 'calling case band1'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'cloud', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'band2', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'band3', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'ring', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'feature', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2

                     case ('band2')
                        !write(*,*) 'calling case band2'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'band1', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'cloud', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'band3', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'ring', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'feature', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2

                     case ('band3')
                        !write(*,*) 'calling case band3'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'band1', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'band2', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'cloud', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'ring', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'feature', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2

                     case ('ring')
                        !write(*,*) 'calling case ring'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'band1', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'band2', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'band3', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'cloud', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'feature', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2

                     case ('feature')
                        !write(*,*) 'calling case feature'
                        par_inds(2) = zodi_model%get_par_ind(comp_str = 'band1', param = 'n_0')
                        prev_vals(2) = zodi_model%theta_stat(par_inds(2),samp_group)
                        zodi_model%theta_stat(par_inds(2),samp_group) = -2

                        par_inds(3) = zodi_model%get_par_ind(comp_str = 'band2', param = 'n_0')
                        prev_vals(3) = zodi_model%theta_stat(par_inds(3),samp_group)
                        zodi_model%theta_stat(par_inds(3),samp_group) = -2

                        par_inds(4) = zodi_model%get_par_ind(comp_str = 'band3', param = 'n_0')
                        prev_vals(4) = zodi_model%theta_stat(par_inds(4),samp_group)
                        zodi_model%theta_stat(par_inds(4),samp_group) = -2

                        par_inds(5) = zodi_model%get_par_ind(comp_str = 'ring', param = 'n_0')
                        prev_vals(5) = zodi_model%theta_stat(par_inds(5),samp_group)
                        zodi_model%theta_stat(par_inds(5),samp_group) = -2

                        par_inds(6) = zodi_model%get_par_ind(comp_str = 'cloud', param = 'n_0')
                        prev_vals(6) = zodi_model%theta_stat(par_inds(6),samp_group)
                        zodi_model%theta_stat(par_inds(6),samp_group) = -2
                     end select
                     deallocate(free_comp, free_par)

                     call params_to_model(zodi_model, theta, samp_group)

                     !evaluate zodi with this condition
                     call wall_time(t3)
                     call get_zodi_emission(data(i)%tod, data(i)%tod%scans(scan)%d(j)%downsamp_pix, &
                          & scan, j, zodi_model, s_zodi, use_lowres_pointing=.true.)
                     call wall_time(t4)
                     zodi_n0 = s_zodi
                     !write(*,*) 'zodi_n0'

                     !set theta_stat[n0]=prev_vals again so that in the next run is going to optimize it
                     do c=1, 6
                        zodi_model%theta_stat(par_inds(c),samp_group) = prev_vals(c)
                     end do

                     call params_to_model(zodi_model, theta, samp_group)

                     !compute gradient
                     grad = sum( &
                        !-(data-tot_model)
                        & -(data(i)%tod%scans(scan)%d(j)%downsamp_tod &
                        &   - data(i)%tod%scans(scan)%d(j)%gain * &
                        &     (data(i)%tod%scans(scan)%d(j)%downsamp_sky &
                        &   + zodi_mod & 
                        &   + mono)) &
                        !*2/sigma**2
                        & *(2.d0/(data(i)%tod%scans(scan)%d(j)%N_psd%sigma0)**2) & 
                        !*zodi_mod with n0=1
                        & *zodi_n0)

                     grad_tot = grad_tot + grad !ha senso?
                     !call wall_time(t4)

                     ndof = ndof + size(data(i)%tod%scans(scan)%d(j)%downsamp_tod)
                  end do
               end do
            end do

            deallocate(s_zodi)

            !write(*,*) 'h4'

            ! Reduce grad to root process
            call mpi_reduce(grad_tot, grad, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, data(1)%tod%comm, ierr)
            call mpi_reduce(ndof, ndof_tot, 1, MPI_INTEGER, MPI_SUM, 0, data(1)%tod%comm, ierr)

            call wall_time(t4)

            norm = 1.0_dp/ndof_tot
            !norm = 1.0_dp

            if (cpar%myid_chain == 0) then
               grad_lnL_zodi(g) = grad*norm     !statistical scaling
               !write(*,*) 'unscaled'
               !grad is computed using the physical parameter theta*scale but hmc works better with just theta
               !so the idea is to apply the chain rule to get the dl/d(t*s)=s*dl/dt
               grad_lnL_zodi = scale*grad_lnL_zodi    !physical scaling
               write(*,*) 'grad_lnL_zodi(',g,')=', grad_lnL_zodi(g)
               call wall_time(t2)
            end if
            
            call wall_time(t2)

         end do 
 
      end subroutine grad_lnL_zodi_hmc

   end subroutine minimize_zodi

end module zodi_hmc_mod
