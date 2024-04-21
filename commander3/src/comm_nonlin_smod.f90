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
submodule (comm_nonlin_mod) comm_nonlin_smod

contains

  module subroutine sample_nonlin_params(cpar, iter, handle, handle_noise)
    !
    ! Routine that loops through all components and samples the spectral parameters that are defined with
    ! non-zero RMS values
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! handle_noise: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handles are updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle, handle_noise    

    integer(i4b) :: i, j, p
    real(dp)     :: t1, t2
    logical(lgt) :: samp_cg
    class(comm_comp),    pointer :: c    => null()
    character(len=512), allocatable, dimension(:) :: comp_labels

    call wall_time(t1)

    call update_status(status, "sample_nonlin_params")

    ! Sample calibration factors
    do i = 1, numband
       if (.not. data(i)%sample_gain .or. index(data(i)%gain_comp, 'firas') .ne. 0) cycle
       call sample_gain(cpar%operation, i, cpar%outdir, cpar%mychain, iter, mod(iter,cpar%resamp_hard_gain_prior_nth_iter)==0, handle)
    end do

    ! Update mixing matrices if gains have been sampled
    if (any(data%sample_gain)) then
       c => compList
       do while (associated(c))
          call c%updateMixmat
          c => c%nextComp()
       end do
    end if

    !call output_FITS_sample(cpar, 200+iter, .true.)


    !return
                    
    c => compList
    do while (associated(c))
       if (c%npar == 0) then
          c => c%nextComp()
          cycle
       end if
                    
       do j = 1, c%npar
          if (c%p_gauss(2,j) == 0.d0) cycle
          select type (c)
          class is (comm_diffuse_comp)

             !check if any poltype has been sampled in a way that breaks the Gibbs chain
             samp_cg = .false.

             !lmax_ind_pol is the lmax of poltype index p, for spec. ind. j 
             if (any(c%lmax_ind_pol(1:c%poltype(j),j) >= 0)) then
                call sample_specind_alm(cpar, iter, handle, c%id, j)
                if (cpar%almsamp_pixreg) then
                  do p = 1,min(c%nmaps, c%poltype(j))
                   if (cpar%almsamp_priorsamp_frozen .and. &
                        & any(c%fix_pixreg(:c%npixreg(p,j),p,j) .eqv. .true.)) then
                      samp_cg = .true.
                   end if
                 end do
                end if
             end if
             

          class is (comm_line_comp) !these codes should (maybe) not need to change
             call sample_specind_local(cpar, iter, handle, c%id, j)
          class is (comm_ptsrc_comp)
             if (j == 1) then
                call sample_specind_local(cpar, iter, handle, c%id, j)
             end if
          end select

       end do
       
       !go to next component
       c => c%nextComp()
           
    end do


    call wall_time(t2)
    if (cpar%myid_chain == 0) write(*,*) '| CPU time specind = ', real(t2-t1,sp)

    
  end subroutine sample_nonlin_params


  module subroutine sample_specind_alm(cpar, iter, handle, comp_id, par_id)
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, r, q, p, pl, np, nlm, l_, m_, idx, delta, burnin
    integer(i4b) :: nsamp, out_every, check_every, num_accepted, smooth_scale, id_native, ierr, ind, nalm_tot_reg
    integer(i4b) :: p_min, p_max, nalm_tot, pix, region
    real(dp)     :: t1, t2, ts, dalm, thresh
    real(dp)     :: mu, par, accept_rate, diff, chisq_prior, chisq_jeffreys, chisq_temp
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    character(len=3) :: tag
    character(len=15) :: regfmt
    character(len=9) :: ar_tag
    character(len=1000) :: outmessage
    character(len=512) :: filename

    logical :: accepted, doexit, optimize, apply_prior
    class(comm_mapinfo), pointer :: info => null()
    class(comm_mapinfo), pointer :: info_theta => null()
    class(comm_map),     pointer :: theta => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_map),     pointer :: theta_smooth => null() ! Spectral parameter of one poltype index (too be smoothed)
    class(comm_comp),    pointer :: c    => null()
    type(map_ptr),     allocatable, dimension(:) :: df

    real(dp),          allocatable, dimension(:,:,:)  :: alms, regs, buffer3
    real(dp),          allocatable, dimension(:)      :: buffer, buffer2, rgs, chisq, theta_pixreg_prop, theta_delta_prop
    integer(c_int),    allocatable, dimension(:)      :: maxit
    real(dp)     :: theta_max, theta_min
    logical      :: outside_limit


    !  Subroutine to sample the (non-linear) diffuse component spectral parameters
    !  using an MCMC alm sampler, rather than pixel-by-pixel (or local) sampling.
    !
    !  Some specifications of the behaviour of the alm-sampler is defined in the 
    !  Commander parameter file, see documentation.
    !
    !  Returns the sampled alms of the diffuse component's spectral parameter.
    !  This is done internally through updating the component's alm directly.
    !  There are no return arguments in this routine, except for the RNG handle.
    !
    !  Arguments (fixed):
    !  ------------------
    !  cpar: comm_params
    !     a class containing all parameters read in from the Commander parameter file 
    !  iter: integer(i4b)
    !     Gibbs chain sample number.
    !  handle: planck_rng
    !     Random number generator handle (or current seed)
    !  comp_id: integer(i4b)
    !     Component id number of the component being sampled. Reference in the compList
    !  par_id: integer(i4b)
    !     id number for the spectral parameter to be sampled in the given component.
    

    ! Sample spectral parameter (parid) for the given signal component
    c => compList
    do while (c%id /= comp_id)
       c => c%nextComp()
    end do

    if (c%p_gauss(2,par_id) == 0.d0) return
    allocate(status_fit(numband))  

    select type (c)
    class is (comm_diffuse_comp)
       
       j = par_id !quick fix to only sample spec. ind. parameter par_id
             
       ! Set up smoothed data
       if (cpar%myid_chain == 0) write(*,*) '|   Sampling '//trim(c%label)//' '//trim(c%indlabel(j))
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
            & 3*c%x%info%nside, 1, .false.)
       theta => comm_map(info_theta)

       ! Params
       out_every = 10
       check_every = 25 !100
       nsamp = cpar%almsamp_nsamp !2000
       burnin = cpar%almsamp_burnin ! Gibbs iter burnin. Tunes steplen.
       optimize = cpar%almsamp_optimize
       apply_prior = cpar%almsamp_apply_prior
       thresh = FLOAT(check_every)*0.8d0 !40.d0 ! 40.d0
       theta_min = c%p_uni(1,par_id) !hard lower prior on theta (i.e. parameter) 
       theta_max = c%p_uni(2,par_id) !hard upper prior on theta (i.e. parameter) 

       if (info%myid == 0 .and. c%L_read(j)) then
          write(*,*) "| Sampling with cholesky matrix"
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
          allocate(buffer(c%nalm_tot), buffer2(c%nalm_tot))
          buffer2 = alms(0,:,pl)
          call mpi_allreduce(buffer2, buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
          alms(0,:,pl) = buffer
          if (cpar%almsamp_pixreg) regs(0,:,pl) = c%theta_pixreg(:,pl,j)
          deallocate(buffer, buffer2)
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
          !if (pl==1) cycle 

          ! p already calculated if larger than poltype 
          if (pl > c%poltype(j)) cycle

          ! p to be sampled with a local sampler 
          if (c%lmax_ind_pol(pl,j) < 0) cycle


          if (cpar%almsamp_pixreg) then
             allocate(theta_pixreg_prop(0:c%npixreg(pl,j))) 
             allocate(theta_delta_prop(0:c%npixreg(pl,j))) 
             allocate(rgs(0:c%npixreg(pl,j))) ! Allocate random vecto
             
             ! Formatter for region output
             write(regfmt,'(I0)') size(c%theta_pixreg(1:,pl,j))
             regfmt = '(a,'//adjustl(trim(regfmt))//'(f7.3))'
             if (cpar%myid_chain == 0) then
               allocate(buffer(c%npixreg(pl,j)))
               buffer = c%pixreg_priors(:c%npixreg(pl,j),pl,j)
               write(*,regfmt) ' | using region priors', buffer
               deallocate(buffer)
             end if
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
             !write(*,*) 'a -- mask'
             !call compute_chisq(c%comm, chisq_fullsky=chisq(0), mask=c%indmask, lowres_eval=.true., evalpol=.true.)
             call compute_chisq(c%comm, chisq_fullsky=chisq(0), mask=c%indmask)
          else
             !write(*,*) 'a -- no mask'
             !call compute_chisq(c%comm, chisq_fullsky=chisq(0), lowres_eval=.true., evalpol=.true.)
             call compute_chisq(c%comm, chisq_fullsky=chisq(0))
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
                      !chisq_prior = chisq_prior + (((c%theta_pixreg(p,pl,j) - c%p_gauss(1,j))/c%p_gauss(2,j))**2)
                      chisq_prior = chisq_prior + (((c%theta_pixreg(p,pl,j) - c%pixreg_priors(p,pl,j))/c%p_gauss(2,j))**2)
                   end do
                end if

                ! Output init sample
                write(*,*) " chisq: " , chisq(0)
                write(*,fmt='(a, i6, a, e12.2, a, f6.2, a, 3f7.2)') " | # sample: ", 0, " - chisq: " , chisq(0), " prior: ", chisq_prior,  " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
                if (cpar%almsamp_pixreg) write(*,fmt=regfmt) " | regs:", real(c%theta_pixreg(1:,pl,j), sp)

                chisq(0) = chisq(0) + chisq_prior 
                !chisq(0) = chisq_prior ! test2

             else 
                write(*,fmt='(a, i6, a, e12.2, a, 3f7.2)') " |# sample: ", 0, " - chisq: " , chisq(0),  " - a_00: ", alms(0,0,:)/sqrt(4.d0*PI)
             end if

          end if

          ! Sample new alms (Account for poltype)
          num_accepted = 0

          nalm_tot = (c%lmax_ind_pol(pl,j) + 1)**2
          do i = 1, nsamp                   

             ! Gather alms from threads to alms array with correct indices
             call gather_alms(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)

             ! Send all alms to 0 (Dont allreduce because only root will do calculation)
             allocate(buffer(c%nalm_tot), buffer2(c%nalm_tot))
             buffer2 = alms(i,:,pl)
             call mpi_reduce(buffer2, buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)
             alms(i,:,pl) = buffer
             deallocate(buffer, buffer2)

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

                ! Broadcast proposed alms from root
                allocate(buffer(c%nalm_tot))
                buffer = alms(i,:,pl)
                call mpi_bcast(buffer, c%nalm_tot, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)                   
                alms(i,:,pl) = buffer
                deallocate(buffer)

             else
                ! --------- region sampling start
                !c%theta_pixreg(c%npixreg(pl,j),pl,j) = 0.d0 ! Just remove the last one for safe measure
                if (info%myid == 0) then
                   q = 0 
                   outside_limit = .true.
                   do while (outside_limit) 
                      q = q + 1
                      ! Save old values
                      theta_pixreg_prop = c%theta_pixreg(:c%npixreg(pl,j),pl,j)
                   
                      rgs = 0.d0
                      do p = 1, c%npixreg(pl,j)
                         rgs(p) = c%steplen(pl,j)*rand_gauss(handle)     
                      end do
                   
                      ! Only propose change to regions not frozen
                      theta_delta_prop = matmul(c%L(:c%npixreg(pl,j), :c%npixreg(pl,j), pl, j), rgs)  !0.05d0*rgs
                      do p = 1, c%npixreg(pl,j)
                         if (.not. c%fix_pixreg(p,pl,j)) theta_pixreg_prop(p) = theta_pixreg_prop(p) + theta_delta_prop(p)
                      end do
                      
                      if (all(theta_pixreg_prop < theta_max) .and. all(theta_pixreg_prop > theta_min)) outside_limit = .false.
                      
                      if (q >= 1000) then !just to not get stucked close to a hard limit
                         theta_pixreg_prop = c%theta_pixreg(:c%npixreg(pl,j),pl,j) !no proposed change
                         outside_limit = .false.
                      end if
                   end do
                end if

                call mpi_bcast(theta_pixreg_prop, c%npixreg(pl,j)+1, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)

                if (.not. (any(c%ind_pixreg_arr(:,pl,j) > 1) ) .and. (c%npixreg(pl,j)>1) ) write(*,*) "Bug in pixreg init"
                ! Loop over pixels in region
                !if (info%myid==0) write(*,*) size(c%ind_pixreg_arr(1,:,j)), size(c%ind_pixreg_arr(1,pl,:)), pl, j
                do pix = 0, theta%info%np-1
                   !if (info%myid==0) write(*,*) c%ind_pixreg_arr(pix,pl,j), theta_pixreg_prop(c%ind_pixreg_arr(pix,pl,j)), theta%map(pix,1), info%myid, pl, pix
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

                !threshold theta map on uniform priors (in case of ringing; done after smoothing)
                theta_smooth%map = min(c%p_uni(2,j),max(c%p_uni(1,j),theta_smooth%map)) 

                call theta_smooth%YtW_scalar
                call mpi_allreduce(theta_smooth%info%nalm, nalm_tot_reg, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
                allocate(buffer3(0:1, 0:nalm_tot_reg-1, 2)) ! Denne er nalm for 1! trenger alle
                
                call gather_alms(theta_smooth%alm, buffer3, theta_smooth%info%nalm, theta_smooth%info%lm, 0, 1, 1)
                call mpi_allreduce(MPI_IN_PLACE, buffer3, nalm_tot_reg, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                alms(i,:,pl) = buffer3(0,0:c%nalm_tot-1,1)
                deallocate(buffer3)

                call theta_smooth%dealloc(); deallocate(theta_smooth)
                ! ------- region sampling end
             end if

             ! Save to correct poltypes
             if (c%poltype(j) == 1) then      ! {T+E+B}
                do q = 1, c%theta(j)%p%info%nmaps
                   alms(i,:,q) = alms(i,:,pl) ! Save to all maps
                   call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)
                end do
             else if (c%poltype(j) == 2) then ! {T,E+B}
                if (pl == 1) then
                   call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                else
                   do q = 2, c%theta(j)%p%info%nmaps
                      alms(i,:,q) = alms(i,:,pl)
                      call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)
                   end do
                end if
             else if (c%poltype(j) == 3) then ! {T,E,B}
                call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
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
                !call compute_chisq(c%comm, chisq_fullsky=chisq_temp, mask=c%indmask, lowres_eval=.true., evalpol=.true.)
                call compute_chisq(c%comm, chisq_fullsky=chisq_temp, mask=c%indmask)
             else
                !call compute_chisq(c%comm, chisq_fullsky=chisq_temp, lowres_eval=.true., evalpol=.true.)
                call compute_chisq(c%comm, chisq_fullsky=chisq_temp)
             end if

             chisq(i) = chisq_temp

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
                         !chisq_prior = chisq_prior + (((theta_pixreg_prop(p) - c%p_gauss(1,j))/c%p_gauss(2,j))**2)
                         chisq_prior = chisq_prior + (((theta_pixreg_prop(p) - c%pixreg_priors(p,pl,j))/c%p_gauss(2,j))**2)
                      end do
                      !write(*,*) "prior ", chisq_prior
                   end if

                   chisq(i) = chisq(i) + chisq_prior
                   !chisq(i) = chisq_prior !test2
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

                ! Reject if proposed values are outside of range
                if (any(theta_pixreg_prop(1:) > c%p_uni(2,j)) .or. any(theta_pixreg_prop(1:) < c%p_uni(1,j))) then
                   accepted = .false.
                   write(*,fmt='(a, f7.3, f7.3, a)') " | Proposed value outside range: ", c%p_uni(1,j), c%p_uni(2,j), ", rejected."
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
                
                ! Output chisq and diff and mean alm
                write(*,*) 'chisq = ', chisq(i), chisq_prior
                write(outmessage,fmt='(a, i6, a, e16.8, a, e16.8, a, e16.8)') "| "//tag, i, " - chisq: " , chisq(i)-chisq_prior, " ", chisq_prior, " diff: ", diff
                write(*,*) adjustl(trim(ar_tag)//trim(outmessage)//trim(achar(27)//'[0m'))
                !write(*,*) trim(outmessage)

                ! Output region information
                if (cpar%almsamp_pixreg) then
                   regs(i,:,pl) = c%theta_pixreg(:,pl,j)
                   write(outmessage,fmt=regfmt) "| regs:", theta_pixreg_prop(1:)
                   write(*,*) adjustl(trim(ar_tag)//trim(outmessage)//trim(achar(27)//'[0m'))
                   !write(*,*) trim(outmessage)
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
                      call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, q, q)
                   end do
                else if (c%poltype(j) == 2) then ! {T,E+B}
                   if (pl == 1) then
                      alms(i,:,pl) = alms(i-1,:,pl) 
                      call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, pl, 1)
                   else
                      do q = 2, c%theta(j)%p%info%nmaps
                         alms(i,:,q) = alms(i-1,:,pl)
                         call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, q, q)                            
                      end do
                   end if
                else if (c%poltype(j) == 3) then ! {T,E,B}
                   alms(i,:,pl) = alms(i-1,:,pl)
                   call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i-1, pl, pl)
                end if

                ! Revert results
                call c%updateMixmat                                   
             end if


             if (info%myid == 0) then 
                ! Output log to file
                allocate(buffer2(c%nalm_tot))
                buffer2 = alms(i,:,pl)
                ! [Maksym]: 
                ! Comment this out, because otherwise I am getting
                ! At line 735 of file comm_nonlin_mod.f90 (unit = 69, file =
                ! 'chains_maksymb/nonlin-samples_synch_beta.dat')
                ! "Fortran runtime error: End of record"
                !write(69, *) iter, tag, i, chisq(i), buffer2
                deallocate(buffer2)
                write(66, *) iter, tag, i, chisq(i), c%theta_pixreg(:, pl, j)

                ! Write to screen every out_every'th
                if (mod(i,out_every) == 0) then
                   diff = chisq(i-out_every) - chisq(i) ! Output diff
                   write(*,fmt='(a, i3, a, e12.2, a, e9.2, a, f7.2, a, e11.4)') " | "//tag, i, " - chisq: " , chisq(i)-chisq_prior, " ", chisq_prior, " diff: ", diff, " - a00: ", alms(i,0,pl)/sqrt(4.d0*PI)

                   ! Format region info
                   if (cpar%almsamp_pixreg) write(*,fmt=regfmt) " | regs:", real(c%theta_pixreg(1:,pl,j), sp)
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
                   write(*, fmt='(a, i3, a, i4, a, f8.2, a, f5.3, a, f5.2)') " | "//tag, i, " - diff last ", check_every, " ", diff, " - accept rate: ", accept_rate, " - time/sample: ", ts
                   call wall_time(t1)

                   ! Adjust steplen in tuning iteration
                   if (iter <= burnin) then
                      if (accept_rate < 0.2) then                 
                         c%steplen(pl,j) = c%steplen(pl,j)*0.5d0
                         write(*,fmt='(a,f10.5)') " | Reducing steplen -> ", c%steplen(pl,j)
                      !else if (accept_rate > 0.45 .and. accept_rate < 0.55) then
                      !   c%steplen(pl,j) = c%steplen(pl,j)*2.0d0
                      !   write(*,fmt='(a,f10.5)') "Equilibrium - Increasing steplen -> ", c%steplen(pl,j)
                      else if (accept_rate > 0.8) then
                         c%steplen(pl,j) = c%steplen(pl,j)*2.d0
                         write(*,fmt='(a,f10.5)') " | Increasing steplen -> ", c%steplen(pl,j)
                      end if
                   end if

                   ! Exit if threshold in tuning stage (First 2 iterations if not initialized on L)
                   if (c%corrlen(j,pl) == 0 .and. diff < thresh .and. accept_rate > 0.2 .and. i>=500) then
                      doexit = .true.
                      write(*,*) "| Chisq threshold and accept rate reached for tuning iteration", thresh
                   end if
                end if
             end if

             call mpi_bcast(doexit, 1, MPI_LOGICAL, 0, c%comm, ierr)
             if (doexit .or. i == nsamp) then
                if (optimize) c%chisq_min(j,pl) = chisq(i) ! Stop increase in chisq
                if (info%myid == 0 .and. i == nsamp) write(*,*) "| nsamp samples reached", nsamp

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

          if (cpar%almsamp_pixreg) then
             if (cpar%almsamp_priorsamp_frozen .and. &
                  & any(c%fix_pixreg(:c%npixreg(pl,j),pl,j) .eqv. .true.)) then
                !Sample frozen regions using component prior
                if (info%myid == 0) then
                   ! Save old values
                   theta_pixreg_prop = c%theta_pixreg(:,pl,j)

                   do p = 1, c%npixreg(pl,j)
                      !if (c%fix_pixreg(p,pl,j)) theta_pixreg_prop(p) = c%p_gauss(1,j) + rand_gauss(handle)*c%p_gauss(2,j)
                      if (c%fix_pixreg(p,pl,j)) then
                         q = 0
                         outside_limit=.true.
                         do while (outside_limit)
                            q = q + 1
                            !draw a new pixel region value from prior
                            theta_pixreg_prop(p) = c%pixreg_priors(p,pl,j) + rand_gauss(handle)*c%p_gauss(2,j)
                            !check if we are outside hard priors, if so, draw new sample
                            if (theta_pixreg_prop(p) < theta_max .and. theta_pixreg_prop(p) > theta_min) outside_limit = .false.
                            if (q > 1000) then !in case the prior RMS is high and we constantly end up outside hard limits
                               theta_pixreg_prop(p) = c%pixreg_priors(p,pl,j) !set to prior value
                               outside_limit = .false.
                            end if
                         end do
                      end if
                   end do
                end if

                call mpi_bcast(theta_pixreg_prop, c%npixreg(pl,j)+1, MPI_DOUBLE_PRECISION, 0, c%comm, ierr)

                !assign new thetas to frozen
                c%theta_pixreg(:c%npixreg(pl,j),pl,j) = theta_pixreg_prop

                ! Loop over pixels in region
                !if (info%myid==0) write(*,*) size(c%ind_pixreg_arr(1,:,j)), size(c%ind_pixreg_arr(1,pl,:)), pl, j
                do pix = 0, theta%info%np-1
                   !if (info%myid==0) write(*,*) c%ind_pixreg_arr(pix,pl,j), theta_pixreg_prop(c%ind_pixreg_arr(pix,pl,j)), theta%map(pix,1), info%myid, pl, pix
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

                !threshold theta map on uniform priors (in case of ringing; done after smoothing)
                theta_smooth%map = min(c%p_uni(2,j),max(c%p_uni(1,j),theta_smooth%map)) 

                call theta_smooth%YtW_scalar
                call mpi_allreduce(theta_smooth%info%nalm, nalm_tot_reg, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
                allocate(buffer3(0:1, 0:nalm_tot_reg-1, 2)) ! Denne er nalm for 1! trenger alle

                call gather_alms(theta_smooth%alm, buffer3, theta_smooth%info%nalm, theta_smooth%info%lm, 0, 1, 1)
                call mpi_allreduce(MPI_IN_PLACE, buffer3, nalm_tot_reg, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                alms(i,:,pl) = buffer3(0,:c%nalm_tot-1,1)
                deallocate(buffer3)

                call theta_smooth%dealloc(); deallocate(theta_smooth)
             
                ! Save to correct poltypes
                if (c%poltype(j) == 1) then      ! {T+E+B}
                   do q = 1, c%theta(j)%p%info%nmaps
                      alms(i,:,q) = alms(i,:,pl) ! Save to all maps
                      call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)
                   end do
                else if (c%poltype(j) == 2) then ! {T,E+B}
                   if (pl == 1) then
                      call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                   else
                      do q = 2, c%theta(j)%p%info%nmaps
                         alms(i,:,q) = alms(i,:,pl)
                         call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, q, q)
                      end do
                   end if
                else if (c%poltype(j) == 3) then ! {T,E,B}
                   call distribute_alms_nonlin(c%theta(j)%p%alm, alms, c%theta(j)%p%info%nalm, c%theta(j)%p%info%lm, i, pl, pl)
                end if

                ! Update mixing matrix with new alms
                if (c%apply_jeffreys) then
                   call c%updateMixmat(df=df, par=j)
                   call compute_jeffreys_prior(c, df, pl, j, chisq_jeffreys)
                else
                   call c%updateMixmat
                end if
             end if !almsamp_priorsamp_frozen .and. any frozen pixregs
          end if !almsamp_pixreg

          deallocate(rgs)
          if (cpar%almsamp_pixreg) deallocate(theta_pixreg_prop, theta_delta_prop) 
       end do ! End pl


       ! Calculate correlation length and cholesky matrix 
       ! (Only if no longer tuning and not initialized from file)
       if (info%myid == 0 .and. .not. c%L_read(j)) then
          if (c%L_calculated(j)  .and. iter > burnin) then
             write(*,*) "| Computing correlation function"
             do pl = 1, c%theta(j)%p%info%nmaps
                ! Skip signals with poltype tag
                if (c%poltype(j) > 1 .and. cpar%only_pol .and. pl == 1) cycle 
                if (pl > c%poltype(j)) cycle
                if (c%lmax_ind_pol(pl,j) < 0) cycle

                if (cpar%almsamp_pixreg) then
                   call compute_corrlen(regs(:,1:,pl), c%fix_pixreg(:,pl,j), c%npixreg(pl,j), maxit(pl), c%corrlen(j,pl))
                else
                   call compute_corrlen(alms(:,:,pl), c%fix_pixreg(:,pl,j), nalm_tot, maxit(pl), c%corrlen(j,pl))
                end if

                c%L_read(j) = .true.  ! L now exist
                write(*,*) "| Correlation length (< 0.1): ", c%corrlen(j,pl) 
             end do
          else
             ! If L does not exist yet, calculate
             write(*,*) '| Calculating cholesky matrix'

             do pl = 1, c%theta(j)%p%info%nmaps
                if (maxit(pl) == 0) cycle ! Cycle if not sampled
                if (cpar%almsamp_pixreg) then
                   call compute_covariance_matrix(regs(INT(maxit(pl)/2):maxit(pl),:,pl), c%L(:, :, pl, j), .true.)
                else
                   call compute_covariance_matrix(alms(INT(maxit(pl)/2):maxit(pl),0:c%nalm_tot-1,pl), c%L(0:c%nalm_tot-1,0:c%nalm_tot-1,pl,j), .true.)
                end if
             end do
             c%steplen(:,j) = 1.d0 ! Revert steplength, because of new proposal matrix
             c%L_calculated(j) = .true.  ! L now exist
          end if

          ! Write almsamp tuning parameters to file
          filename = trim(cpar%outdir)//'/init_alm_'//trim(c%label)//'_'//trim(c%indlabel(j))//'.dat'
          write(*,*) "| Writing tuning parameters to file: ", trim(filename)
          open(58, file=filename, recl=10000)

          if (maxval(c%corrlen(j,:)) > 0) then
             write(58,*) c%corrlen(j,:)
          else
             write(58,*) nsamp, nsamp, nsamp
          end if
          ! Formatting proposal array
          do p = 1, info%nmaps
             write(58,*) "Proposal matrix L for signal", p
             do q = 0, size(c%L(:,1,p,j))-1
                write(58,fmt='(*(f10.5))') c%L(q,:,p,j)
             end do
          end do

          close(58)

       end if

       ! Clean up
       if (info%myid == 0) close(69)   
       if (info%myid == 0) close(66)   
       deallocate(alms, regs, chisq, maxit)
       call theta%dealloc(); deallocate(theta)

       if (c%apply_jeffreys) then
          do k = 1, numband
             call df(k)%p%dealloc(); deallocate(df(k)%p)
          end do
          deallocate(df)
       end if
    end select
    
    deallocate(status_fit)

  end subroutine sample_specind_alm

  module subroutine sample_specind_local(cpar, iter, handle, comp_id, par_id)
    !
    ! Routine that sets up the sampling using the local sampling routine  for the spectral parameter given by
    ! par_id for the component given by the comp_id parameter. 
    ! Then it calls on the specific sampling routine and finally updates the components spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    integer(i4b),       intent(in)    :: comp_id     !component id, only doing one (!) component 
    integer(i4b),       intent(in)    :: par_id      !parameter index, 1 -> npar (per component)

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx, p_ind, p_min, p_max
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind, ind_pol
    real(dp)     :: t1, t2, ts, dalm
    real(dp)     :: mu, par
    integer(i4b), allocatable, dimension(:) :: status_fit   ! 0 = excluded, 1 = native, 2 = smooth
    integer(i4b)                            :: status_amp   !               1 = native, 2 = smooth
    class(comm_mapinfo), pointer :: info => null()
    class(comm_N),       pointer :: tmp  => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_map),     pointer :: temp_res  => null()
    class(comm_map),     pointer :: temp_map  => null()
    class(comm_map),     pointer :: temp_map2  => null()
    class(comm_comp),    pointer :: c    => null()
    class(comm_map),     pointer :: theta_single_pol => null() ! Spectral parameter of one poltype index (too be smoothed)
    real(dp),          allocatable, dimension(:,:) :: m


    c           => compList     
    do while (c%id /= comp_id)
       c => c%nextComp()
    end do

    ! HKE: Only perform joint Powell search
    if (par_id /= 1 .and. trim(c%type) == "diffuse") then
       if (cpar%myid == 0) write(*,*) 'Performing Powell'
       return
    end if

    call wall_time(t1)
    
    ! Initialize residual maps
    do i = 1, numband
       res             => compute_residual(i)
       data(i)%res%map =  res%map
       call res%dealloc(); deallocate(res)
       nullify(res)
!!$       call data(i)%res%writeFITS('res_'//trim(data(i)%label)//'.fits')
!!$       call mpi_finalize(ierr)
!!$       stop
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
       if (cpar%myid == 0) write(*,*) '|    Sampling ', trim(c%label), ' ', trim(c%indlabel(par_id))

       ! Sample spectral parameters
       call c%sampleSpecInd(cpar, handle, par_id, iter)

    class is (comm_ptsrc_comp)
       if (cpar%myid == 0) write(*,*) '|    Sampling ', trim(c%label)

       ! Sample spectral parameters
       call c%sampleSpecInd(cpar, handle, par_id, iter)

    class is (comm_diffuse_comp)
       if (cpar%myid == 0) write(*,*) '|    Sampling ', trim(c%label), ' ', trim(c%indlabel(par_id))
       call update_status(status, "nonlin start " // trim(c%label)// ' ' // trim(c%indlabel(par_id)))

       ! Set up type of smoothing scale
       id_native    = 0

       ! Compute smoothed residuals
       nullify(info)
       status_amp   = 0
       status_fit   = 0
       smooth_scale = c%smooth_scale(par_id)
       do i = 1, numband
          ! Chooses an index that is polarized so that smoothing can be done
          ! correctly later on.
          if (data(i)%info%nmaps == 3) ind_pol = i
          if (cpar%num_smooth_scales == 0) then
             status_fit(i)   = 1    ! Native
          else
             if (.not. associated(data(i)%N_smooth(smooth_scale)%p) .or. &
                  & maxval(data(i)%bp(0)%p%nu) < c%nu_min_ind(par_id) .or. &
                  & minval(data(i)%bp(0)%p%nu) > c%nu_max_ind(par_id)) then
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
             rms_smooth(i)%p    => data(i)%N
          else if (status_fit(i) == 2) then
             ! Fit is done with downgraded data
             ! Needs rewrite! We should smooth at Nside and lmax of the data band, then downgrade to low resolution Nside. This to limit chances for ringing. 
             info  => comm_mapinfo(data(i)%res%info%comm, data(i)%res%info%nside, &
                  & data(i)%res%info%lmax, data(i)%res%info%nmaps, data(i)%res%info%pol)

             !now we have changed the lmax of B_smooth (and nside) to match that of the data band
             !Smooth the residual map at full resolution (i.e. at data%info%nside) with lmax equal to data%info%lmax
             call smooth_map(info, .false., data(i)%B(0)%p%b_l, data(i)%res, &
                  & data(i)%B_smooth(smooth_scale)%p%b_l, temp_map)

             !create comm_map_info for low resolution residual
             info  => comm_mapinfo(data(i)%res%info%comm, cpar%nside_smooth(smooth_scale), cpar%lmax_smooth(smooth_scale), &
                  & data(i)%res%info%nmaps, data(i)%res%info%pol)
             !down_grade smoothed residual map to smoothing scale Nside
             res_smooth(i)%p => comm_map(info)
             call temp_map%udgrade(res_smooth(i)%p)

             !deallocate fullres smooth rms
             call temp_map%dealloc(); deallocate(temp_map)
             nullify(temp_map)

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
          ! Also needs rewriting:
          ! We must smooth at full resolution, then down-grade.
          ! We add beams to component (B_smooth and B_postproc)
          ! using lmax of component (at full resolution/component Nside)

          info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, c%x%info%lmax, &
               & c%x%info%nmaps, c%x%info%pol)

          !smooth at component Nside and lmax of amplitude
          call smooth_map(info, .true., &
               & c%B_smooth_amp(par_id)%p%b_l*0.d0+1.d0, c%x, &  
               & c%B_smooth_amp(par_id)%p%b_l,           temp_map)

          info  => comm_mapinfo(c%x%info%comm, cpar%nside_smooth(smooth_scale), &
               & cpar%lmax_smooth(smooth_scale), &
               & c%x%info%nmaps, c%x%info%pol)
          c%x_smooth => comm_map(info)
          call temp_map%udgrade(c%x_smooth)

          !deallocate fullres smooth map (temp_map)
          call temp_map%dealloc(); deallocate(temp_map)
          nullify(temp_map)

       end if

       ! Compute smoothed spectral index maps
       allocate(c%theta_smooth(c%npar))
       do k = 1, c%npar
          !if (k == par_id) cycle
          if (status_amp == 1) then ! Native resolution
             info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, &
                  & c%x%info%lmax, c%x%info%nmaps, c%x%info%pol)
             call smooth_map(info, .false., &
                  & data(id_native)%B(0)%p%b_l*0.d0+1.d0, c%theta(k)%p, &  
                  & data(id_native)%B(0)%p%b_l,           c%theta_smooth(k)%p)
          else if (status_amp == 2) then ! Common FWHM resolution

             !smooth at full resolution with the smoothing scale FWHM parameter beam
             !allocate a temporary full resolution smoothed theta map
             info  => comm_mapinfo(c%theta(k)%p%info%comm, &
                  & c%B_smooth_specpar(par_id)%p%info%nside, &
                  & c%B_smooth_specpar(par_id)%p%info%lmax, &
                  & c%theta(k)%p%info%nmaps, c%theta(k)%p%info%nmaps==3)
             temp_map => comm_map(info)

             ! To make sure we smooth spec param as spin zero map, 
             ! we smooth one polarization/map at a time
             ! allocate a single map 
             info  => comm_mapinfo(c%theta(k)%p%info%comm, &
                  & c%B_smooth_specpar(par_id)%p%info%nside, &
                  & c%B_smooth_specpar(par_id)%p%info%lmax, &
                  & 1, .false.)
             temp_res => comm_map(info)
             do j = 1,c%theta(k)%p%info%nmaps
                !copy polarization map to single map
                temp_res%map(:,1)=c%theta(k)%p%map(:,j)
                !smooth single map
                call smooth_map(info, .false., &
                     & c%B_smooth_specpar(par_id)%p%b_l*0.d0+1.d0, temp_res, &  
                     & c%B_smooth_specpar(par_id)%p%b_l,           temp_map2)

                !copy smoothed single map to temporary smoothed map, correct polarization 
                temp_map%map(:,j)=temp_map2%map(:,1)
                !deallocate smoothed single map
                call temp_map2%dealloc(); deallocate(temp_map2)
                nullify(temp_map2)
             end do

             !create the smoothing scale Nside theta map and downgrade parameter map
             info  => comm_mapinfo(c%theta(k)%p%info%comm, &
                  & cpar%nside_smooth(smooth_scale), &
                  & cpar%lmax_smooth(smooth_scale), &
                  & c%theta(k)%p%info%nmaps, c%theta(k)%p%info%pol)
             c%theta_smooth(k)%p => comm_map(info)
             !call temp_map%udgrade(c%theta_smooth(k)%p)
             call c%theta(k)%p%udgrade(c%theta_smooth(k)%p)

             !deallocate fullres smooth parameter maps (temp_map,temp_map2,temp_res)
             call temp_map%dealloc(); deallocate(temp_map)
             nullify(temp_map)
             call temp_res%dealloc(); deallocate(temp_res)
             nullify(temp_res)

          end if
       end do

       ! Sample spectral parameters for diffuse component
       !call sampleDiffuseSpecInd_nonlin(cpar, handle, c%id, par_id, iter)
       !call sampleDiffuseSpecInd_simple(cpar, handle, c%id, par_id, iter)
       call sampleDiffuseSpecInd_powell(cpar, handle, c%id, iter)

    end select
          
          
    ! Clean up temporary data structures
    select type (c)
    class is (comm_line_comp)
    class is (comm_diffuse_comp)
       
       if (associated(c%x_smooth)) then
          call c%x_smooth%dealloc(); deallocate(c%x_smooth)
          nullify(c%x_smooth)
       end if
       do k =1, c%npar
          if (k == par_id) cycle
          if (allocated(c%theta_smooth)) then
             if (associated(c%theta_smooth(k)%p)) then
                call c%theta_smooth(k)%p%dealloc(); deallocate(c%theta_smooth(k)%p)
             end if
          end if
       end do
       if (allocated(c%theta_smooth)) deallocate(c%theta_smooth)
       do i = 1, numband
          if (.not. associated(rms_smooth(i)%p)) cycle
          if (status_fit(i) == 2) then
             call res_smooth(i)%p%dealloc(); deallocate(res_smooth(i)%p)
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
             info  => comm_mapinfo(c%theta(par_id)%p%info%comm, &
                  & c%theta(par_id)%p%info%nside, &
                  & c%B_pp_fr(par_id)%p%info%lmax, 1, .false.) !only want 1 map

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

                   c%theta_smooth(par_id)%p%map = min(c%p_uni(2,par_id), &
                        & max(c%p_uni(1,par_id), &
                        & c%theta_smooth(par_id)%p%map)) ! threshold smoothed map on uniform limits
                   
                   do p = p_min,p_max
                      ! assign smoothed theta map to relevant polarizations
                      c%theta(par_id)%p%map(:,p) = c%theta_smooth(par_id)%p%map(:,1)
                   end do
                   call c%theta_smooth(par_id)%p%dealloc(); deallocate(c%theta_smooth(par_id)%p)
                   deallocate(c%theta_smooth)
                end if
             end do
             call theta_single_pol%dealloc(); deallocate(theta_single_pol)
             theta_single_pol => null()

          end if
       end if

       call update_status(status, "nonlin stop " // trim(c%label)// ' ' // trim(c%indlabel(par_id)))
       if (c%theta(par_id)%p%info%myid == 0 .and. cpar%verbosity > 2) &
            & write(*,*) '| Updating Mixing matrix after sampling '// &
            & trim(c%label)// ' ' // trim(c%indlabel(par_id))

    end select

    ! Update mixing matrix 
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


  module subroutine sample_specind_multi(cpar, iter, handle, comp_labels)
    !
    ! Routine that sets up the sampling using the local sampling routine  for the spectral parameter given by
    ! par_id for the component given by the comp_id parameter. 
    ! Then it calls on the specific sampling routine and finally updates the components spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),  intent(in)    :: cpar
    integer(i4b),       intent(in)    :: iter
    type(planck_rng),   intent(inout) :: handle    
    character(len=512), dimension(:), intent(in)    :: comp_labels

    integer(i4b) :: i, j, k, q, p, pl, np, nlm, l_, m_, idx, p_min, p_max, id, ncomp, npar, pix_out
    integer(i4b) :: nsamp, out_every, num_accepted, smooth_scale, id_native, ierr, ind, pol, pix
    real(dp)     :: t1, t2, ts, dalm
    real(dp)     :: mu
    real(dp)     :: lnL_old, lnL_new, eps
    real(dp), allocatable, dimension(:) :: theta, theta_old
    class(comm_mapinfo), pointer :: info => null()
    class(comm_mapinfo), pointer :: info_lowres => null()
    class(comm_map),     pointer :: res  => null()
    class(comm_map),     pointer :: raw  => null()
    class(comm_map),     pointer :: raw2  => null()
    class(comm_map),     pointer :: raw_lowres  => null()
    class(comm_map),     pointer :: temp_map  => null()
    class(comm_map),     pointer :: chisq_lowres  => null()
    class(comm_map),     pointer :: chisq_old  => null()
    class(comm_comp),    pointer :: c_in    => null()
    class(comm_diffuse_comp),    pointer :: c    => null()
    type(diff_ptr),    allocatable, dimension(:) :: comps

    call wall_time(t1)
    eps = 1d-9
    pix_out = 370000
    
    ncomp = size(comp_labels)
    npar  = 0
    allocate(comps(ncomp))
    do j = 1, size(comp_labels)
       c_in => compList     
       do while (trim(c_in%label) /= trim(comp_labels(j)))
          c_in => c_in%nextComp()
       end do
       select type (c_in)
       class is (comm_diffuse_comp)
          comps(j)%p => c_in 
          npar       = npar + c_in%npar + 1
       end select
    end do

    !create comm_map_info for low resolution residual
    c => comps(1)%p
    smooth_scale =  c%smooth_scale(1) ! Require same smoothing scale for all par
    info_lowres  => comm_mapinfo(c%x%info%comm, cpar%nside_smooth(smooth_scale), &
         & cpar%lmax_smooth(smooth_scale), c%x%info%nmaps, c%x%info%pol)

    ! Initialize residual maps
    !call output_FITS_sample(cpar, 1000+iter, .true.)
    do i = 1, numband
       if (.not. associated(data(i)%N_smooth(smooth_scale)%p)) cycle
       if (data(i)%bp(0)%p%nu_c < c%nu_min_ind(1) .or. &
            & data(i)%bp(0)%p%nu_c > c%nu_max_ind(1)) cycle

       !res             => compute_residual(i,exclude_comps=comp_labels)
       res             => compute_residual(i)
!       call res%writeFITS("res_"//trim(data(i)%label)//"_sig.fits")
       data(i)%res%map =  res%map
       call res%dealloc(); deallocate(res)
       nullify(res)


       raw             => compute_residual(i)
       raw_lowres      => comm_map(info_lowres)
       call raw%writeFITS("res_"//trim(data(i)%label)//"_hires_raw.fits")
       info  => comm_mapinfo(data(i)%res%info%comm, data(i)%res%info%nside, &
            & data(i)%res%info%lmax, data(i)%res%info%nmaps, data(i)%res%info%pol)
       raw2 => comm_map(info)
       call smooth_map(info, .false., data(i)%B(0)%p%b_l, raw, &
            & data(i)%B_smooth(smooth_scale)%p%b_l, raw2)
       call raw2%udgrade(raw_lowres)
       call raw_lowres%writeFITS("res_"//trim(data(i)%label)//"_lowres_raw.fits")

       call raw%dealloc(); deallocate(raw)
       nullify(raw)
       call raw2%dealloc(); deallocate(raw2)
       nullify(raw2)
       call raw_lowres%dealloc(); deallocate(raw_lowres)
       nullify(raw_lowres)
    end do

    ! Add components back into residual
!!$    do j = 1, ncomp
!!$       do i = 1, numband
!!$          allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
!!$          m               = comps(j)%p%getBand(i)
!!$          data(i)%res%map = data(i)%res%map + m
!!$          deallocate(m)
!!$       end do
!!$    end do

    ! Compute smoothed residuals
    nullify(info)
    do i = 1, numband
       if (.not. associated(data(i)%N_smooth(smooth_scale)%p)) cycle
       if (data(i)%bp(0)%p%nu_c < c%nu_min_ind(1) .or. &
            & data(i)%bp(0)%p%nu_c > c%nu_max_ind(1)) cycle

       !Smooth the residual map at full resolution 
       info  => comm_mapinfo(data(i)%res%info%comm, data(i)%res%info%nside, &
            & data(i)%res%info%lmax, data(i)%res%info%nmaps, data(i)%res%info%pol)
       call smooth_map(info, .false., data(i)%B(0)%p%b_l, data(i)%res, &
            & data(i)%B_smooth(smooth_scale)%p%b_l, temp_map)

       !downgrade smoothed residual map to smoothing scale Nside
       rms_smooth(i)%p => data(i)%N_smooth(smooth_scale)%p
       res_smooth(i)%p => comm_map(info_lowres)
       res_lowres(i)%p => comm_map(info_lowres)
       dust_lowres(i)%p => comm_map(info_lowres)
       hotpah_lowres(i)%p => comm_map(info_lowres)
       call temp_map%udgrade(res_smooth(i)%p)
       call res_smooth(i)%p%writeFITS("res_"//trim(data(i)%label)//"_lowres_in.fits")
!       rms_smooth(i)%p => data(i)%N_smooth(smooth_scale)%p
       rms_smooth(i)%p => data(i)%N_smooth(smooth_scale)%p
       !call rms_smooth(i)%writeFITS("rmssmooth_"//trim(data(i)%label)//".fits")

       ! Clean up
       call temp_map%dealloc(); deallocate(temp_map); nullify(temp_map)

       !temp_map => comm_map(info)
       !temp_map%map = comps(1)%p%getBand(i)
       !call temp_map%writeFITS("dust_"//trim(data(i)%label)//"_fullres.fits")
       !call temp_map%dealloc(); deallocate(temp_map); nullify(temp_map)
    end do

    ! Check that there are relevant data
    if (.not. associated(info)) then
       write(*,*) 'Error: No bands contribute to index fit!'
       call mpi_finalize(i)
       stop
    end if

    ! Compute smoothed component maps
    do j = 1, ncomp
       c => comps(j)%p

       ! Compute smoothed amplitude maps
       info  => comm_mapinfo(c%x%info%comm, c%x%info%nside, c%x%info%lmax, &
            & c%x%info%nmaps, c%x%info%pol)
       call smooth_map(info, .true., &
            & c%B_smooth_amp(1)%p%b_l*0.d0+1.d0, c%x, &  
            & c%B_smooth_amp(1)%p%b_l,           temp_map)
       c%x_smooth => comm_map(info_lowres)
       call temp_map%udgrade(c%x_smooth)
       call temp_map%dealloc(); deallocate(temp_map); nullify(temp_map)
       call c%x_smooth%writeFITS("smooth_amp_"//trim(c%label)//".fits")

       ! Compute smoothed spectral index maps; spin zero
       allocate(c%theta_smooth(c%npar))
       do k = 1, c%npar
!!$          info  => comm_mapinfo(c%theta(k)%p%info%comm, &
!!$               & c%B_smooth_specpar(1)%p%info%nside, &
!!$               & c%B_smooth_specpar(1)%p%info%lmax, &
!!$               & c%theta(k)%p%info%nmaps, c%theta(k)%p%info%nmaps==3)
!!$          temp_map => comm_map(info)
!!$          temp_res => comm_map(info)
!!$          temp_res%map=c%theta(k)%p%map
!!$
!!$          call smooth_map(info, .false., &
!!$               & c%B_smooth_specpar(1)%p%b_l*0.d0+1.d0, temp_res, &  
!!$               & c%B_smooth_specpar(1)%p%b_l,           temp_map, &
!!$               & spinzero=.true.)

          c%theta_smooth(k)%p => comm_map(info_lowres)
          !call temp_map%udgrade(c%theta_smooth(k)%p)
          call c%theta(k)%p%udgrade(c%theta_smooth(k)%p)

         call c%theta_smooth(k)%p%writeFITS("smooth_ind_"//trim(c%label)//"_"//trim(c%indlabel(k))//".fits")

!!$          call temp_map%dealloc(); deallocate(temp_map); nullify(temp_map)
!!$          call temp_res%dealloc(); deallocate(temp_res); nullify(temp_res)
       end do
    end do

    ! Output initial residual maps
    pol = 1
    allocate(theta(npar), theta_old(npar))
    do pix = 0, info_lowres%np-1
       i = 1
       do j = 1, ncomp
          c => comps(j)%p          
          theta_old(i) = c%x_smooth%map(pix,pol)
          do k = 1, c%npar
             theta_old(i+k) = min(c%theta_smooth(k)%p%map(pix,pol), c%p_uni(2,k)-eps) 
             theta_old(i+k) = max(theta_old(i+k),                   c%p_uni(1,k)+eps) 
          end do
          i = i + c%npar + 1
       end do
       call compute_lowres_residuals(theta_old)
    end do
    do i = 1, numband
       if (associated(rms_smooth(i)%p)) then
          res_smooth(i)%p%map = res_smooth(i)%p%map + dust_lowres(i)%p%map
          res_smooth(i)%p%map = res_smooth(i)%p%map + hotpah_lowres(i)%p%map
       end if
    end do

    do i = 1, numband
       if (associated(rms_smooth(i)%p)) then
          call res_lowres(i)%p%writeFITS("res_"//trim(data(i)%label)//"_lowres_pre.fits")
          call hotpah_lowres(i)%p%writeFITS("hotpah_"//trim(data(i)%label)//"_lowres_pre.fits")
          call dust_lowres(i)%p%writeFITS("dust_"//trim(data(i)%label)//"_lowres_pre.fits")
       end if
    end do


    ! Perform the actual fit, pixel-by-pixel
    pol = 1
    chisq_lowres => comm_map(info_lowres)
    chisq_old => comm_map(info_lowres)
    do pix = 0, info_lowres%np-1
       
       if (info_lowres%pix(pix+1) == pix_out) then
          open(58,file="pix.dat")
          do k = 1, numband
             if (.not. associated(rms_smooth(k)%p)) cycle
             if (trim(data(k)%unit) == 'uK_cmb') then
                write(58,*) data(k)%bp(0)%p%nu_c, res_smooth(k)%p%map(pix,pol)/data(k)%bp(0)%p%a2t, rms_smooth(k)%p%rms_pix(pix,pol)/data(k)%bp(0)%p%a2t
             else 
                write(58,*) data(k)%bp(0)%p%nu_c, res_smooth(k)%p%map(pix,pol)/data(k)%bp(0)%p%a2f, rms_smooth(k)%p%rms_pix(pix,pol)/data(k)%bp(0)%p%a2f
             end if
          end do
          close(58)
       end if

       i = 1
       do j = 1, ncomp
          c => comps(j)%p          
          theta_old(i) = c%x_smooth%map(pix,pol)
          do k = 1, c%npar
             theta_old(i+k) = min(c%theta_smooth(k)%p%map(pix,pol), c%p_uni(2,k)-eps) 
             theta_old(i+k) = max(theta_old(i+k),                   c%p_uni(1,k)+eps) 
          end do
          i = i + c%npar + 1
       end do
       theta = theta_old

       if (info_lowres%pix(pix+1) == pix_out) lnL_old = lnL_multi(theta_old)
       chisq_old%map(pix,pol) = 2*lnL_multi(theta_old)
       do i = 1, 10
          call powell(theta, lnL_multi, ierr)
       end do
       chisq_lowres%map(pix,pol) = 2*lnL_multi(theta)
       if (info_lowres%pix(pix+1) == pix_out) then
          lnL_new = lnL_multi(theta)
          write(*,*) info_lowres%pix(pix+1), theta_old
          write(*,*) info_lowres%pix(pix+1), theta
          write(*,fmt='(i8,2e16.8)') info_lowres%pix(pix+1), lnL_old, lnL_new
       end if
       call compute_lowres_residuals(theta)

       i = 1
       do j = 1, ncomp
          c => comps(j)%p
          c%x_smooth%map(pix,pol) = theta(i)
          do k = 1, c%npar
             c%theta_smooth(k)%p%map(pix,pol) = theta(i+k)
          end do
          i = i + c%npar + 1
       end do
    end do
    call chisq_lowres%writeFITS("chisq_lowres.fits")
    call chisq_old%writeFITS("chisq_old.fits")
    do i = 1, numband
       if (associated(rms_smooth(i)%p)) then
          call res_lowres(i)%p%writeFITS("res_"//trim(data(i)%label)//"_lowres_post.fits")
       end if
    end do

    do j = 1, ncomp
       c => comps(j)%p
       do k = 1, c%npar
          call c%theta_smooth(k)%p%udgrade(c%theta(k)%p)
          call c%theta(k)%p%writeFITS(trim(c%label)//trim(c%indlabel(k))//".fits")
       end do
       call c%x_smooth%writeFITS(trim(c%label)//"_amp_lowres.fits")
!!$       call c%x_smooth%udgrade(c%x)
!!$       call c%x%YtW
    end do

    deallocate(theta, theta_old)

    ! Smooth if requested
    if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
       !write(*,*) 'skal ikke vaere her'
       do j = 1, ncomp
          c => comps(j)%p
          do k = 1, c%npar
             !info  => c%theta(k)%p%info
             info  => comm_mapinfo(c%theta(k)%p%info%comm, c%theta(k)%p%info%nside, &
                  & 3*c%theta(k)%p%info%nside, c%theta(k)%p%info%nmaps, c%theta(k)%p%info%pol)

             temp_map => comm_map(info)
             call smooth_map(info, .false., &
                  & c%B_pp_fr(1)%p%b_l*0.d0+1.d0, c%theta(k)%p, &  
                  & c%B_pp_fr(1)%p%b_l, temp_map, spinzero=.true.)
             c%theta(k)%p%map = temp_map%map
             call temp_map%dealloc(); deallocate(temp_map); nullify(temp_map)
          end do
       end do
    end if


    ! Update mixing matrix and residuals
    do j = 1, ncomp
       c => comps(j)%p
       call c%updateMixmat 
!!$       do i = 1, numband
!!$          allocate(m(0:data(i)%info%np-1,data(i)%info%nmaps))
!!$          m               = c%getBand(i)
!!$          data(i)%res%map = data(i)%res%map - m
!!$          deallocate(m)
!!$       end do
    end do

    !call output_FITS_sample(cpar, 1100+iter, .true.)

!!$    do i = 1, numband
!!$       res             => compute_residual(i)
!!$       call res%writeFITS("res_"//trim(data(i)%label)//"_new.fits")
!!$       call res%dealloc(); deallocate(res)
!!$       nullify(res)
!!$    end do


    ! Clean up temporary data structures
    do j = 1, ncomp
       c => comps(j)%p
       if (associated(c%x_smooth)) then
          call c%x_smooth%dealloc(); deallocate(c%x_smooth); nullify(c%x_smooth)
       end if
       do k =1, c%npar
          if (associated(c%theta_smooth(k)%p)) then
             call c%theta_smooth(k)%p%dealloc(); deallocate(c%theta_smooth(k)%p)
          end if
       end do
       if (allocated(c%theta_smooth)) deallocate(c%theta_smooth)
    end do

    do i = 1, numband
       if (associated(rms_smooth(i)%p)) then
          call res_smooth(i)%p%dealloc()
          deallocate(res_smooth(i)%p); nullify(res_smooth(i)%p)
       end if
    end do
    deallocate(comps)   

!!$    call mpi_finalize(ierr)
!!$    stop

  contains

    function lnL_multi(x)
      use healpix_types
      implicit none
      real(dp), dimension(:), intent(in), optional :: x
      real(dp)             :: lnL_multi
      
      integer(i4b) :: i, j, k, n
      real(dp)     :: s, res, sigma, amp
      
      ! Check priors
      i = 1
      do j = 1, ncomp
         c => comps(j)%p
         do k = 1, c%npar
            if (x(i+k) < c%p_uni(1,k) .or. x(i+k) > c%p_uni(2,k)) then
               lnL_multi = 1.d30
               return
            end if
         end do
         i = i + c%npar + 1
      end do

      ! Compute chi-square term
      lnL_multi      = 0.d0
      do k = 1, numband
         if (.not. associated(rms_smooth(k)%p)) cycle
         s = 0
         i = 1
         do j = 1, ncomp
            c => comps(j)%p
            amp = x(i) !c%x_smooth%map(pix,pol)
            s   = s + amp * c%F_int(1,k,0)%p%eval(x(i+1:i+c%npar))&
                 & * data(k)%gain * c%cg_scale(pol)
            if (info_lowres%pix(pix+1) == pix_out) then
               write(*,fmt='(a,a,i3,f16.4)') '      sig', trim(data(k)%label), j, s
            end if
            i = i + c%npar + 1
         end do
         sigma = rms_smooth(k)%p%rms_pix(pix,pol)
         !sigma = 0.03d0 * res_smooth(k)%p%map(pix,pol) 
         res = res_smooth(k)%p%map(pix,pol) - s
         lnL_multi = lnL_multi - 0.5d0 * res**2 / sigma**2
          if (info_lowres%pix(pix+1) == pix_out) then
            write(*,fmt='(a,i4,2f10.3,f16.3)') '  chisq', k, res, sigma, res**2 / sigma**2
         end if
      end do

      ! Add Gaussian prior
      i = 1
      do j = 1, ncomp
         c => comps(j)%p
         do k = 1, c%npar
            if (c%p_gauss(2,k) > 0.d0) then
               lnL_multi = lnL_multi -0.5d0*((x(i+k)-c%p_gauss(1,k))/c%p_gauss(2,k))**2
               if (info_lowres%pix(pix+1) == pix_out) then
                  write(*,*) 'prior', j, k, x(i+k), c%p_gauss(1,k), c%p_gauss(2,k), ((x(i+k)-c%p_gauss(1,k))/c%p_gauss(2,k))**2
               end if
            end if
         end do
         i = i + c%npar + 1
      end do

      ! Switch sign, since powell is a minimization routine
      lnL_multi = -lnL_multi

      if (info_lowres%pix(pix+1) == pix_out) then
         write(*,fmt='(a,7f7.2)') 'multipowell', real(x,sp), 2*lnL_multi
      end if
         
    end function lnL_multi

    subroutine compute_lowres_residuals(x)
      use healpix_types
      implicit none
      real(dp), dimension(:), intent(in), optional :: x
      
      real(dp)             :: lnL_multi
      
      integer(i4b) :: i, j, k, n
      real(dp)     :: s, s_tot, res, sigma, amp
      

      ! Compute chi-square term
      lnL_multi      = 0.d0
      do k = 1, numband
         if (.not. associated(rms_smooth(k)%p)) cycle
         s_tot = 0
         i = 1
         do j = 1, ncomp
            c => comps(j)%p
            amp = x(i) !c%x_smooth%map(pix,pol)
            !if (j == 1) then
               !write(*,*) k, pix, real(x(i+1:i+c%npar),sp)
               s   = amp * c%F_int(1,k,0)%p%eval(x(i+1:i+c%npar))&
                    & * data(k)%gain * c%cg_scale(pol)
               if (j == 1) then
                  dust_lowres(k)%p%map(pix,pol) = s
               else
                  hotpah_lowres(k)%p%map(pix,pol) = s
               end if
               s_tot = s_tot + s
            !end if
            i = i + c%npar + 1
         end do
         !res_lowres(k)%p%map(pix,pol) = res_smooth(k)%p%map(pix,pol) - s_tot
      end do

    end subroutine compute_lowres_residuals


  end subroutine sample_specind_multi
  

  module subroutine gather_alms(alm, alms, nalm, lm, i, pl, pl_tar)
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

  module subroutine distribute_alms_nonlin(alm, alms, nalm, lm, i, pl, pl_tar)
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

  end subroutine distribute_alms_nonlin

  module subroutine compute_corrlen(x, fix, n, maxit, corrlen)
    implicit none

    real(dp), dimension(:,:),    intent(in)    :: x
    logical(lgt), dimension(:),  intent(in)      :: fix        
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
    corrlen = corrlen_init
    allocate(C_(delta))
    allocate(N_(delta))
    
    !open(58, file='correlation_function.dat', recl=10000)
          
    ! Calculate correlation function per parameter
    do p = 1, n
       if (fix(p)) cycle ! Skip fixed regions

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
    !close(58)
  end subroutine compute_corrlen


  ! Sample spectral parameters using straight Powell or inversion sampler
  module subroutine sampleDiffuseSpecInd_simple(cpar, handle, comp_id, par_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: par_id, comp_id, iter

    integer(i4b) :: pix, ierr, pol, status
    real(dp)     :: x(3), lnL_old, lnL_new, eps, theta_old
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()

    eps = 1d-9

    c           => compList     
    do while (comp_id /= c%id)
       c => c%nextComp()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select


!    do pol = 1, c_lnL%theta_smooth(par_id)%p%info%nmaps
    pol = 1
    do pix = 0, c_lnL%theta_smooth(par_id)%p%info%np-1
       if (c_lnL%theta_smooth(par_id)%p%map(pix,pol) <= c_lnL%p_uni(1,par_id)+eps) then
          x(1) = c_lnL%p_uni(1,par_id)+eps
          x(3) = min(x(1) + c_lnL%p_gauss(2,par_id), c_lnL%p_uni(2,par_id)-eps)
          x(2) = 0.5d0*(x(1)+x(3))
       else if (c_lnL%theta_smooth(par_id)%p%map(pix,pol) >= c_lnL%p_uni(2,par_id)-eps) then
          x(3) = c_lnL%p_uni(2,par_id)-eps
          x(1) = max(x(3) - c_lnL%p_gauss(2,par_id), c_lnL%p_uni(1,par_id)+eps)
          x(2) = 0.5d0*(x(1)+x(3))
       else
          x(2) = c_lnL%theta_smooth(par_id)%p%map(pix,pol)
          x(1) = max(x(2) - c_lnL%p_gauss(2,par_id), 0.5d0*(x(2) + c_lnL%p_uni(1,par_id)))
          x(3) = min(x(2) + c_lnL%p_gauss(2,par_id), 0.5d0*(x(2) + c_lnL%p_uni(2,par_id)))
          if (x(2) == x(1) .or. x(2) == x(3)) x(2) = 0.5d0*(x(1)+x(3))
       end if
       theta_old = x(2)

       if (mod(pix,1000) == 0) lnL_old = lnL_simple(c_lnL%theta_smooth(par_id)%p%map(pix,pol))
       c_lnL%theta_smooth(par_id)%p%map(pix,pol) = &
            & sample_InvSamp(handle, x, lnL_simple, prior=c_lnL%p_uni(:,par_id), &
            & optimize=(trim(cpar%operation)=='optimize'), status=status)
       if (mod(pix,1000) == 0) then
          lnL_new = lnL_simple(c_lnL%theta_smooth(par_id)%p%map(pix,pol))
          write(*,fmt='(i8,2f8.3,2e16.8)') c_lnL%theta_smooth(par_id)%p%info%pix(pix+1), theta_old, c_lnL%theta_smooth(par_id)%p%map(pix,pol), lnL_old, lnL_new
       end if
    end do
    !end do
!!$call mpi_finalize(ierr)
!!$    stop
    call c_lnL%theta_smooth(par_id)%p%udgrade(c_lnL%theta(par_id)%p)

  contains

    function lnL_simple(x)
      use healpix_types
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: lnL_simple
      
      integer(i4b) :: i, j, n
      real(dp)     :: s, res, theta(6)
      real(dp), allocatable, dimension(:,:) :: f_precomp
      
      ! Check priors
      if (x < c_lnL%p_uni(1,par_id) .or. x > c_lnL%p_uni(2,par_id)) then
         lnL_simple = 1.d30
         return
      end if
      
      ! Set up internal parameter array
      n = c_lnL%npar
      do i = 1, n
         theta(i) = c_lnL%theta_smooth(i)%p%map(pix,pol)
         theta(i) = max(theta(i), c_lnL%p_uni(1,i)+eps)
         theta(i) = min(theta(i), c_lnL%p_uni(2,i)-eps)
      end do
      theta(par_id) = x
      
      ! Compute chi-square term
      lnL_simple      = 0.d0
      do j = 1, numband
         if (.not. associated(rms_smooth(j)%p)) cycle
         s   = c_lnL%x_smooth%map(pix,pol) * c_lnL%F_int(1,j,0)%p%eval(theta(1:n))&
              & * data(j)%gain * c_lnL%cg_scale(pol)
         res = res_smooth(j)%p%map(pix,pol) - s
         lnL_simple = lnL_simple - 0.5d0 * res**2 / rms_smooth(j)%p%rms_pix(pix,pol)**2
      end do

      ! Add Gaussian prior
      if (c_lnL%p_gauss(2,par_id) > 0.d0) then
         lnL_simple = lnL_simple -0.5d0*((x-c_lnL%p_gauss(1,par_id))/c_lnL%p_gauss(2,par_id))**2
      end if

    end function lnL_simple
    
  end subroutine sampleDiffuseSpecInd_simple


  ! Sample spectral parameters using straight Powell or inversion sampler
  module subroutine sampleDiffuseSpecInd_powell(cpar, handle, comp_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, iter

    integer(i4b) :: i, pix, ierr, pol, status
    real(dp)     :: lnL_old, lnL_new, eps
    real(dp), allocatable, dimension(:) :: theta, theta_old
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()

    eps = 1d-9

    c           => compList     
    do while (comp_id /= c%id)
       c => c%nextComp()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select


!    do pol = 1, c_lnL%theta_smooth(par_id)%p%info%nmaps
    pol = 1
    allocate(theta(c_lnL%npar), theta_old(c_lnL%npar))
    do pix = 0, c_lnL%theta_smooth(1)%p%info%np-1
       
       do i = 1, c_lnL%npar
          theta_old(i) = min(c_lnL%theta_smooth(i)%p%map(pix,pol), c_lnL%p_uni(2,i)-eps) 
          theta_old(i) = max(theta_old(i),                         c_lnL%p_uni(1,i)+eps) 
       end do
       theta = theta_old

       if (mod(pix,1000) == 0) lnL_old = lnL_simple(theta_old)
       call powell(theta, lnL_simple, ierr)
       if (mod(pix,1000) == 0) then
          lnL_new = lnL_simple(theta)
          write(*,fmt='(i8,2e16.8)') c_lnL%theta_smooth(1)%p%info%pix(pix+1), lnL_old, lnL_new
       end if

       do i = 1, c_lnL%npar
          c_lnL%theta_smooth(i)%p%map(pix,pol) = theta(i)
       end do
    end do
    !end do
!!$call mpi_finalize(ierr)
!!$    stop
    do i = 1, c_lnL%npar
       call c_lnL%theta_smooth(i)%p%udgrade(c_lnL%theta(i)%p)
    end do

    deallocate(theta, theta_old)

  contains

    function lnL_simple(x)
      use healpix_types
      implicit none
      real(dp), dimension(:), intent(in), optional :: x
      real(dp)             :: lnL_simple
      
      integer(i4b) :: i, j, n
      real(dp)     :: s, res
      real(dp), allocatable, dimension(:,:) :: f_precomp
      
      ! Check priors
      do i = 1, c_lnL%npar
         if (x(i) < c_lnL%p_uni(1,i) .or. x(i) > c_lnL%p_uni(2,i)) then
            lnL_simple = 1.d30
            return
         end if
      end do

      ! Compute chi-square term
      lnL_simple      = 0.d0
      do j = 1, numband
         if (.not. associated(rms_smooth(j)%p)) cycle
         s   = c_lnL%x_smooth%map(pix,pol) * c_lnL%F_int(1,j,0)%p%eval(x)&
              & * data(j)%gain * c_lnL%cg_scale(pol)
         res = res_smooth(j)%p%map(pix,pol) - s
         lnL_simple = lnL_simple - 0.5d0 * res**2 / rms_smooth(j)%p%rms_pix(pix,pol)**2
      end do

      ! Add Gaussian prior
      do i = 1, c_lnL%npar
         if (c_lnL%p_gauss(2,i) > 0.d0) then
            lnL_simple = lnL_simple -0.5d0*((x(i)-c_lnL%p_gauss(1,i))/c_lnL%p_gauss(2,i))**2
         end if
      end do

      ! Switch sign, since powell is a minimization routine
      lnL_simple = -lnL_simple

      if (c_lnL%theta_smooth(1)%p%info%pix(pix+1) == 100000) then
         write(*,*) 'powell', real(x,sp), 2*lnL_simple
      end if
         
    end function lnL_simple
    
  end subroutine sampleDiffuseSpecInd_powell


  !Here comes all subroutines for sampling diffuse components locally
  ! Sample spectral parameters
  module subroutine sampleDiffuseSpecInd_nonlin(cpar, handle, comp_id, par_id, iter)
    !
    ! Overarching routine that sets up the sampling of diffuse type component spectral parameters
    ! 
    ! Calls on the specific sampling routine and finally updates the component's spectral parameter map
    ! 
    !
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       a parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! Returns:
    !       No explicit parameter is returned.
    !       The RNG handle is updated as they are used and returned from the routine
    !       All other changes are done internally
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: par_id, comp_id, iter

    integer(i4b) :: i, j, k, n, p, q, pix, ierr, ind(1), n_ok, id
    integer(i4b) :: npar, np, nmaps
    real(dp)     :: t1, t2, delta_lnL_threshold
    real(dp)     :: mu
    real(dp)     :: theta_min, theta_max
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    real(dp),     allocatable, dimension(:,:) :: buffer_lnL
    !Following is for the local sampler
    integer(i4b) :: p_min, p_max
    class(comm_mapinfo),            pointer :: info => null()

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    first_call          = .false.


    id = par_id
    c           => compList     
    do while (comp_id /= c%id)
       c => c%nextComp()
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


       allocate(buffer_lnL(0:c_lnL%theta(id)%p%info%np-1,c_lnL%theta(id)%p%info%nmaps))
       buffer_lnL=max(min(c_lnL%theta(id)%p%map,theta_max),theta_min) 
       
       
       do p = 1,c_lnL%poltype(id)
          if (c_lnL%lmax_ind_pol(p,id) >= 0) cycle !this set of polarizations are not to be local sampled (is checked before this point)
          if (c_lnL%poltype(id) > 1 .and. cpar%only_pol .and. p == 1) cycle !only polarization (poltype > 1)
          if (p > c_lnL%nmaps) cycle ! poltype > number of maps
          ! Return if all prior RMS's are zero

          call wall_time(t1)
          if (c_lnL%pol_pixreg_type(p,id) /= 0) then
             if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) '| Sampling poltype index', p, &
                  & 'of ', c_lnL%poltype(id) !Needed?
             call sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
             call wall_time(t2)
             if (info%myid == 0 .and. cpar%verbosity > 1) write(*,*) '| poltype:',c_lnL%poltype(id),' pol:', &
                  & p,'CPU time specind = ', real(t2-t1,sp)
          else
             write(*,*) '| Undefined spectral index sample region'
             write(*,*) '| Component: ',trim(c_lnL%label),', ind: ',trim(c_lnL%indlabel(id))
             stop
          end if



       end do

       !after sampling is done we assign the spectral index its new value(s)
       c_lnL%theta(id)%p%map = buffer_lnL !unsmoothed map (not post_proc smoothed)

       ! Ask for CG preconditioner update
       if (c_lnL%cg_unique_sampgroup > 0) recompute_diffuse_precond = .true.

       ! deallocate

       deallocate(buffer_lnL)

    end select

  end subroutine sampleDiffuseSpecInd_nonlin



  module subroutine sampleDiffuseSpecIndPixReg_nonlin(cpar, buffer_lnL, handle, comp_id, par_id, p, iter)
    !
    ! Routine that samples diffuse type component spectral parameters in pixel regions
    ! 
    ! Arguments:
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! iter: integer
    !       Gibb's sample counter
    !
    ! handle: planck_rng type
    !       A parameter for the RNG to produce random numbers
    !
    ! comp_id: integer
    !       Integer ID for the specific component to be sampled (in the list of the active components)
    !
    ! par_id: integer
    !       Integer ID for the specific spectral parameter to be sampled in the component given by 'comp_id'
    !
    ! buffer_lnL: double precision array
    !       An array copy of the current spectral parameter map, truncated by the absolute parameter limits.
    !       Array with dimension (0:npix-1,nmaps), where npix is the number of pixels given by the components 
    !       resolution parameter (Nside, see HEALPix), and nmaps is the number of polarizations (1 if only Temperature;
    !       3 if polarization is included, TQU)
    !
    ! p: integer
    !       Index counter for the polarization type that is to be sampled. Sets what map polarization(s) to be sampled.
    !
    ! Returns:
    !       No explicit parameter is returned, except for the sampled spectral parameter through the 
    !       'buffer_lnL' parameter.
    !       The RNG handle is updated as it is used and returned from the routine.
    !       All other changes are done internally.
    !
    implicit none
    type(comm_params),                       intent(in)           :: cpar
    real(dp),               dimension(0:,:), intent(inout)        :: buffer_lnL
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: comp_id, par_id
    integer(i4b),                            intent(in)           :: p       !incoming polarization
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b) :: i, j, k, l, n, q, pr, pix, ierr, ind(1), n_ok, id
    integer(i4b) :: i_mono,l_mono,m_mono, N_theta_MC, l_count
    integer(i4b) :: npar, np, nmaps, nsamp
    real(dp)     :: a, b, t0, t1, t2, delta_lnL_threshold
    real(dp)     :: mu, sigma, mu_p, sigma_p
    real(dp)     :: theta_min, theta_max, running_accept, running_dlnL
    real(dp)     :: running_correlation, correlation_limit
    logical(lgt) :: outside, sample_mono
    logical(lgt), save :: first_call = .true.
    class(comm_comp),         pointer :: c => null()
    class(comm_comp),         pointer :: c2 => null()
    class(comm_diffuse_comp), pointer :: c_lnL => null()
    class(comm_md_comp),      pointer :: c_mono => null()
    !Following is for the local sampler
    real(dp)     :: mixing_old, mixing_new, lnL_new, lnL_old, delta_lnL, lnL_prior, lnL_init
    real(dp)     :: accept_rate, accept_scale, proplen, avg_dlnL, lnL_total_init
    real(dp)     :: old_theta, new_theta, prior_rms_scaling
    integer(i4b) :: p_min, p_max, pixreg_nprop, band_count, pix_count, burn_in
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len, out_every
    integer(i4b) :: npixreg, smooth_scale, arr_ind, np_lr, np_fr, myid_pix, unit
    logical(lgt) :: first_sample, loop_exit, use_det, burned_in, sampled_nprop, sampled_proplen, first_nprop, sample_accepted
    character(len=512) :: filename, postfix, monocorr_type
    character(len=6) :: itext
    character(len=4) :: ctext
    character(len=2) :: pind_txt
    character(len=512), dimension(1000) :: tokens
    real(dp),      allocatable, dimension(:) :: all_thetas, data_arr, invN_arr, mixing_new_arr
    real(dp),      allocatable, dimension(:) :: old_thetas, new_thetas, init_thetas
    real(dp),      allocatable, dimension(:) :: theta_corr_arr
    real(dp),      allocatable, dimension(:) :: new_theta_smooth, dlnL_arr
    real(dp),  allocatable, dimension(:,:,:) :: theta_MC_arr
    integer(i4b),  allocatable, dimension(:) :: band_i, pol_j, accept_arr
    class(comm_mapinfo),             pointer :: info_fr => null() !full resolution
    class(comm_mapinfo),             pointer :: info_fr_single => null() !full resolution, nmaps=1
    class(comm_mapinfo),             pointer :: info_lr => null() !low resolution
    class(comm_mapinfo),             pointer :: info_lr_single => null() !low resolution, nmaps=1
    class(comm_map),                 pointer :: theta_single_fr => null() ! Spectral parameter of polt_id index
    class(comm_map),                 pointer :: theta_single_lr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_lr_hole => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: theta_fr => null() ! smoothed spec. ind. of polt_id index
    class(comm_map),                 pointer :: mask_lr => null() ! lowres mask
    class(comm_map),                 pointer :: mask_mono => null() ! lowres mask
    class(comm_map),                 pointer :: res_map => null() ! lowres  residual map
    class(comm_map),                 pointer :: ones_map => null() ! Constant sky
    class(comm_map),                 pointer :: Ninv_map => null() ! Inverse covariance diagonal
    class(comm_map),                 pointer :: temp_map => null() 
    class(comm_map),                 pointer :: temp_chisq => null() 
    type(map_ptr), allocatable, dimension(:) :: df
    real(dp),      allocatable, dimension(:) :: monopole_val, monopole_rms, monopole_mu, monopole_mixing
    real(dp),      allocatable, dimension(:) :: old_mono, new_mono
    real(dp), allocatable, dimension(:,:,:,:) :: multipoles_trace !trace for proposed and accepted multipoles. (2,Nsamp,numband,0:3), the final dimension can be adjusted if only monopoles are estimated
    logical(lgt),  allocatable, dimension(:) :: monopole_active
    real(dp),    allocatable, dimension(:,:) :: reduced_data
    real(dp),    allocatable, dimension(:,:) :: harmonics, harmonics2
    real(dp),                 dimension(0:3) :: multipoles, md_b   
    real(dp),             dimension(0:3,0:3) :: md_A
    real(dp),                 dimension(3)   :: vector
    integer(i4b) :: i_md, j_md, k_md, max_prop

    c           => compList     
    do while (comp_id /= c%id)
       c => c%nextComp()
    end do
    select type (c)
    class is (comm_diffuse_comp)
       c_lnL => c !to be able to access all diffuse comp parameters through c_lnL
    end select

    id = par_id !hack to not rewrite too much from diffuse_comp_mod

    call update_status(status, "nonlin pixreg samling start " // trim(c_lnL%label)// ' ' // trim(c_lnL%indlabel(par_id)))

    info_fr  => comm_mapinfo(c_lnL%theta(id)%p%info%comm, c_lnL%theta(id)%p%info%nside, &
         & c_lnL%B_pp_fr(id)%p%info%lmax, c_lnL%theta(id)%p%info%nmaps, c_lnL%theta(id)%p%info%pol)

    info_fr_single  => comm_mapinfo(c_lnL%theta(id)%p%info%comm, &
         & c_lnL%theta(id)%p%info%nside, &
         & c_lnL%B_pp_fr(id)%p%info%lmax, 1, .false.)

    info_lr  => comm_mapinfo(c_lnL%x_smooth%info%comm, c_lnL%x_smooth%info%nside, &
         & c_lnL%x_smooth%info%lmax, c_lnL%x_smooth%info%nmaps, c_lnL%x_smooth%info%pol)

    info_lr_single  => comm_mapinfo(c_lnL%x_smooth%info%comm, c_lnL%x_smooth%info%nside, &
         & c_lnL%x_smooth%info%lmax, 1, .false.)
                
    myid_pix = info_fr%myid ! using this proc id and info_fr%comm in all mpi calls to ensure that values are properly dispersed between processors 

    delta_lnL_threshold = 25.d0
    n                   = 101
    n_ok                = 50
    burn_in             = 50
    burned_in           = .false.
    first_call          = .false.    
    n_prop_limit        = 100
    n_corr_limit        = n_prop_limit
    out_every           = 100

    if (c_lnL%spec_corr_convergence(id)) then
       correlation_limit = c_lnL%spec_corr_limit(id)
    else
       correlation_limit = 0.1d0 !setting a default correlation limit if non have been defined for the spectral parameter
    end if
    !     buffer_lnL (theta, limited by min/max)

    ! Check to see which polarization types we care about
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

    !###################################################################################################
    ! Here we must collect monopoles if spectral parameter is to be sampled with monopoles
    ! we do this while checking which bands to include in the sampling

    !set up which bands and polarizations to include
    allocate(band_i(3*numband),pol_j(3*numband))

    !set up the monopoles if applicable
    if (c_lnL%spec_mono_combined(par_id)) then
       allocate(monopole_val(numband))
       allocate(monopole_mu(numband))
       allocate(monopole_rms(numband))
       allocate(monopole_mixing(numband))
       allocate(monopole_active(numband))
       monocorr_type=trim(c_lnL%spec_mono_type(par_id))
       monopole_active=.false.
       monopole_val=0.d0
       monopole_mu=0.d0
       monopole_rms=0.d0
       monopole_mixing=0.d0

       !ud_grade monopole mask (if it is not same nside as smoothing scale)
       mask_mono => comm_map(info_lr)
       mask_mono%map = 0d0
       call c_lnL%spec_mono_mask(par_id)%p%udgrade(mask_mono) !ud_grade monopole mask to nside of smoothing scale
       where (mask_mono%map > 0.5d0)
          mask_mono%map=1.d0
       elsewhere
          mask_mono%map=0.d0
       end where

       ! Creating map of all ones. Useful in general.
       ones_map  => comm_map(info_lr)
       call c_lnL%spec_mono_mask(par_id)%p%udgrade(ones_map)
       ones_map%map      = 0.d0
       ones_map%map(:,1) =  1.d0

       ! Creating chisq map
       temp_chisq => comm_map(info_lr)
       call c_lnL%spec_mono_mask(par_id)%p%udgrade(temp_chisq)
       temp_chisq%map = 0.d0

       !set up harmonics matrices for solving mono- and dipole estimates
       allocate(harmonics(0:np_lr-1,0:3)) !harmonics without mixing scaling (will not chainge)
       allocate(harmonics2(0:np_lr-1,0:3)) !harmonics with mixing (will change)
       do i = 0, np_lr-1
          call pix2vec_ring(info_lr%nside, info_lr%pix(i+1), vector) !important to get the correct pixel number, i.e. "info_lr%pix(+1i)", not "i". The +1 is because the pix array starts from index 1 (not 0)
          
          if (mask_mono%map(i,1) > 0.5d0) then
             harmonics(i,0) = 1.d0
             harmonics(i,1) = vector(1)
             harmonics(i,2) = vector(2)
             harmonics(i,3) = vector(3)
          else
             harmonics(i,:) = 0.d0
          end if

       end do

    end if

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

       !if monopole sampling, read in all neccessary values from the monopole of the band
       if (c_lnL%spec_mono_combined(par_id)) then
          if (p_min > 1) cycle !only polarization maps active, i.e. no monopoles are to be sampled

          c2           => compList     

          do while (associated(c2))
             select type (c2)
             class is (comm_md_comp)
                c_mono => c2 !to be able to access all diffuse comp parameters through c_lnL
                if (trim(c_mono%label) == trim(data(k)%label)) then
                   do i_mono = 0, c_mono%x%info%nalm-1
                      call c_mono%x%info%i2lm(i_mono,l_mono,m_mono)
                      if (l_mono == 0) then ! Monopole
                         monopole_val(k)=c_mono%x%alm(i_mono,1) / sqrt(4.d0*pi) 
                         monopole_mu(k)=c_mono%mu%alm(i_mono,1) / sqrt(4.d0*pi)
                      end if
                   end do

                   monopole_rms(k)=sqrt(c_mono%Cl%Dl(0,1) / (4.d0*pi))
                   monopole_mixing(k)=c_mono%RJ2unit_(1)
                   monopole_active(k)=.true.
                   if (c_mono%mono_from_prior) monopole_active(k)=.false. !we should not sample monopoles that are sampled from prior

                   call get_tokens(trim(c_lnL%spec_mono_freeze(par_id)), ",", tokens, n)

                   do i_mono = 1, n
                      if (trim(data(k)%label) == trim(tokens(i_mono)) ) then
                         monopole_active(k)=.false. !we freeze monopoles the user wants to freeze
                         exit
                      end if
                   end do
                   exit !exit the while-assiciated-loop, we have found the necessary information
                end if
             end select

             c2 => c2%nextComp()
          end do
          
       end if

    end do

    if (c_lnL%spec_mono_combined(par_id)) then
       ! since the monopole value and prior mean alms are spread on the different processors, only one of them 
       ! has the monopole value, We therefor perfom a mpi_allreduce with a sum so that all get the correct value
       call mpi_allreduce(MPI_IN_PLACE, monopole_val, numband, MPI_DOUBLE_PRECISION, MPI_SUM, info_fr%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, monopole_mu, numband, MPI_DOUBLE_PRECISION, MPI_SUM, info_fr%comm, ierr)

       if (.true.) then !debugging
          if (cpar%verbosity>2 .and. myid_pix==0) then
             write(*,*) '|   Monopoles of active bands in sampling '
             write(*,*) '|   label         band_number  mono[uK_RJ]  mu[uK_RJ]  rms[uK_RJ]    mixing'
             do i = 1,numband
                if (monopole_active(i)) write(*,fmt='(a,a15,i13,e13.3,e11.3,e11.3,e11.4,i6)') &
                     & ' |',trim(data(i)%label),i,monopole_val(i),monopole_mu(i),monopole_rms(i),monopole_mixing(i)
             end do

             write(*,*) '|   '
             write(*,*) '|   And in band units'
             write(*,*) '|   label         band_number  mono        mu        rms'
             do i = 1,numband
                if (monopole_active(i)) write(*,fmt='(a,a15,i13,e13.3,e11.3,e11.3, i6)') &
                     & ' |',trim(data(i)%label),i,monopole_val(i)*monopole_mixing(i),monopole_mu(i)*monopole_mixing(i), &
                     & monopole_rms(i)*monopole_mixing(i)
             end do

          end if
       end if

    end if

    !###################################################################################################

    if (trim(c_lnL%pol_lnLtype(p,id)) == 'prior') then !Prior sample
       allocate(new_thetas(0:npixreg))
       new_thetas = c_lnL%theta_pixreg(:npixreg,p,id)
       do pr = 1,npixreg
          ! Don't even think about sampling regions which are held fixed
          if (c_lnL%pol_pixreg_type(p,id) == 3) then
             if (c_lnL%fix_pixreg(pr,p,id)) cycle
          end if
          if (myid_pix == 0) then
             new_theta = c_lnL%theta_prior(1,p,id) + c_lnL%theta_prior(2,p,id)*rand_gauss(handle)
             new_theta = min(theta_max,max(theta_min,new_theta))
          end if
          
          !broadcast new_theta
          call mpi_bcast(new_theta, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
          if (myid_pix == 0) then
             write(*,*) '| Sampled new value using Gaussian prior'
             write(*,fmt='(a,i6)') ' |  Pixel region = ',pr
             write(*,fmt='(a,f7.3)') ' |  Prior mean = ',c_lnL%theta_prior(1,p,id)
             write(*,fmt='(a,f10.6)') ' |  Prior RMS  = ',c_lnL%theta_prior(2,p,id)
             write(*,fmt='(a,f10.6)') ' |  Old value  = ',c_lnL%theta_pixreg(pr,p,id)
             write(*,fmt='(a,f10.6)') ' |  New value  = ',new_theta
          end if
          c_lnL%theta_pixreg(pr,p,id) = new_theta
          new_thetas(pr) = new_theta

       end do
       do pix=0,np_fr-1
          do i = p_min, p_max
             buffer_lnL(pix,i)=new_thetas(c_lnL%ind_pixreg_arr(pix,p,id))
          end do
       end do
       deallocate(new_thetas)

       !###################################################################################################
       ! If we sample with monopoles, we must calculate the new band monopoles, given the new theta
       sample_mono=.true.
       !if all pixel regions are frozen (or no RMS), we do not resample monopoles as theta does not change
       if (c_lnL%pol_pixreg_type(p,id) == 3) then
          if (c_lnL%theta_prior(2,p,id) == 0.d0) then
             sample_mono = .false.
          else if (all(c_lnL%fix_pixreg(1:npixreg,p,id))) then
             sample_mono = .false.
          end if
       end if

       if (c_lnL%spec_mono_combined(par_id) .and. sample_mono) then
          !allocate and set new full resolution theta map
          theta_fr => comm_map(info_fr_single)
          theta_fr%map(:,1) = buffer_lnL(:,1)
          
          !then we smooth new theta map to lower resolution
          smooth_scale = c_lnL%smooth_scale(id)
          if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
             if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then !smooth to correct resolution
                call smooth_map(info_fr_single, .false., &
                     & c_lnL%B_pp_fr(id)%p%b_l*0.d0+1.d0, theta_fr, &  
                     & c_lnL%B_pp_fr(id)%p%b_l, temp_map)
                theta_lr_hole => comm_map(info_lr_single)
                call temp_map%udgrade(theta_lr_hole)
                call temp_map%dealloc(); deallocate(temp_map)
                nullify(temp_map)

             else !no postproc smoothing, ud_grade to correct resolution
                theta_lr_hole => comm_map(info_lr_single)
                call theta_fr%udgrade(theta_lr_hole)
             end if
          else !native resolution, fr = lr
             theta_lr_hole => comm_map(info_lr_single)
             theta_lr_hole%map(:,1) = theta_fr%map(:,1)
          end if


          ! produce reduced data sets of the active bands in Temperature with active monopole
          ! then sample a new monopole
          allocate(reduced_data(0:np_lr-1,1))
          allocate(all_thetas(npar))
          reduced_data=0.d0
          if (cpar%verbosity>2 .and. myid_pix==0) then
             !print info on the change in monopoles from initial monopoles 
             write(*,*) '|   Change in monopoles of active bands in sampling. In band units. '
             write(*,*) '|   band_label    mono_in        mono_out       prior_mean     prior_rms'
          end if

          do j = 1,numband
             if (monopole_active(j)) then
                do pix = 0,np_lr-1
                   do i = 1, npar
                      if (i == id) cycle
                      all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                   end do

                   ! get conversion factor from amplitude to data (i.e. mixing matrix element)         
                   all_thetas(id)=theta_lr_hole%map(pix,1)
                   mixing_old = c_lnL%F_int(1,j,0)%p%eval(all_thetas) * &
                        & data(j)%gain * c_lnL%cg_scale(1)

                   !remove the parameter signal (given new theta) and add monopole back into map
                   reduced_data(pix,1) = res_smooth(j)%p%map(pix,1) - mixing_old* &
                        & c_lnL%x_smooth%map(pix,1) + monopole_mixing(j)*monopole_val(j)
                end do

                if (trim(monocorr_type) == 'monopole+dipole' .or. &
                     & trim(monocorr_type) == 'monopole-dipole') then
                   ! resample the monopole (and dipole) given the new reduced data

                   md_A = 0.d0
                   md_b = 0.d0
                   harmonics2=harmonics*monopole_mixing(j)
                   do j_md = 0, 3
                      do k_md = 0, 3
                         md_A(j_md,k_md) = sum(harmonics2(:,j_md) * harmonics2(:,k_md)) 
                      end do

                      md_b(j_md) = sum(reduced_data(:,1) * harmonics2(:,j_md)) !is to be set later, this will change
                   end do

                   multipoles=0.d0
                   !we need to run an MPI reduce to get all harmonics for md_A and md_b
                   call mpi_allreduce(MPI_IN_PLACE, md_A, 16, MPI_DOUBLE_PRECISION, MPI_SUM, info_lr%comm, ierr)
                   call mpi_allreduce(MPI_IN_PLACE, md_b, 4, MPI_DOUBLE_PRECISION, MPI_SUM, info_lr%comm, ierr)

                   !solve the mono-/dipole system
                   call solve_system_real(md_A, multipoles, md_b) 
                   ! Need to get the statistical power for when adding the monopole prior

                   ! Not all cores get low resolution pixels. These if tests
                   ! avoid that
                   if (np_lr > 0) then
                       ones_map%map = 0d0
                       ones_map%map(:,1)  = 1d0
                       ones_map%map(:,1) = ones_map%map(:,1) * mask_mono%map(:,1)
                       ones_map%map(:,2:) = 0d0
                       call rms_smooth(j)%p%sqrtInvN(ones_map)
                       ones_map%map = ones_map%map * monopole_mixing(j)
                       a = sum(ones_map%map(:,1)**2)
                   else 
                       a = 0d0
                   end if
                   call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, & 
                        & MPI_SUM, info_lr%comm, ierr)

                else if (trim(monocorr_type) == 'monopole') then

                   a=0.d0
                   b=0.d0
                   if (np_lr > 0) then
                       ones_map%map = 0d0
                       ones_map%map(:,1) = 1d0
                       ones_map%map(:,1) = ones_map%map(:,1) * mask_mono%map(:,1)
                       call rms_smooth(j)%p%sqrtInvN(ones_map)
                       ones_map%map = ones_map%map * monopole_mixing(j)
                       a = sum(ones_map%map(:,1)**2)

                       ones_map%map = reduced_data
                       call rms_smooth(j)%p%sqrtInvN(ones_map)
                       ones_map%map = ones_map%map * monopole_mixing(j)
                       b = sum(ones_map%map(:,1)**2)
                   end if
                   call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, & 
                        & MPI_SUM, info_lr%comm, ierr)
                   call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, & 
                        & MPI_SUM, info_lr%comm, ierr)


                   if (a > 0.d0) then
                      multipoles(0) = b/a
                   else
                      multipoles(0) = 0.d0
                   end if

                end if

                if (info_lr%myid == 0) then

                   if (a > 0.d0) then !we have statistical power to estimate a monopole
                      sigma=sqrt(a)
                      mu = multipoles(0)
                   else if (monopole_rms(j) > 0.d0) then !this will effectively set monopole to prior mean
                      mu=0.d0
                      sigma = 0.d0
                   else
                      mu=0.d0
                      sigma = 0.d0
                   end if

                   ! applying prior 
                   if (monopole_rms(j) > 0.d0) then
                      sigma_p = 1.d0/monopole_rms(j)
                      mu_p = monopole_mu(j)
                      mu = (mu*sigma**2 + mu_p*sigma_p**2)/(sigma**2 + sigma_p**2)
                   end if
                end if
                
                ! bcast new monopole to other processors
                call mpi_bcast(mu, 1, MPI_DOUBLE_PRECISION, 0, info_lr%comm, ierr)
                if (cpar%verbosity>2 .and. myid_pix==0) then
                   if (monopole_active(j)) write(*,fmt='(a,a15,e15.5,e15.5,e15.5,e15.5)') &
                        & ' |',trim(data(j)%label),monopole_val(j)*monopole_mixing(j), &
                        & mu*monopole_mixing(j),monopole_mu(j)*monopole_mixing(j), &
                        & monopole_rms(j)*monopole_mixing(j)
                end if
                monopole_val(j)=mu
                
                !find the monopole component and update the monopole
                c2           => compList     
                do while (associated(c2))
                   select type (c2)
                   class is (comm_md_comp)
                      c_mono => c2 !to be able to access all diffuse comp parameters through c_lnL
                      if (trim(c_mono%label) == data(j)%label) then
                         do i_mono = 0, c_mono%x%info%nalm-1
                            call c_mono%x%info%i2lm(i_mono,l_mono,m_mono)
                            if (l_mono == 0) then ! Monopole
                               c_mono%x%alm(i_mono,1) = sqrt(4.d0*pi) * monopole_val(j) 
                            end if
                         end do
                         exit !exit loop
                      end if
                   end select

                   c2 => c2%nextComp()
                end do

             end if
          end do

          call theta_lr_hole%dealloc(); deallocate(theta_lr_hole)
          theta_lr_hole => null()
          call theta_fr%dealloc(); deallocate(theta_fr)
          theta_fr => null()
          
       end if

       if (allocated(monopole_val)) deallocate(monopole_val)
       if (allocated(monopole_mu)) deallocate(monopole_mu)
       if (allocated(monopole_rms)) deallocate(monopole_rms)
       if (allocated(monopole_mixing)) deallocate(monopole_mixing)
       if (allocated(monopole_active)) deallocate(monopole_active)
       if (allocated(reduced_data)) deallocate(reduced_data)
       if (allocated(all_thetas)) deallocate(all_thetas)
       if (allocated(harmonics)) deallocate(harmonics)
       if (allocated(harmonics2)) deallocate(harmonics2)


       !###################################################################################################

       return
    end if
    ! Now the prior sampling section is done
    !###################################################################################################

    if (band_count==0) then
       buffer_lnL(:,p_min:p_max)=c_lnL%p_gauss(1,id) !set theta to prior, as no bands are valid, no data
       deallocate(band_i,pol_j)
       if (myid_pix == 0 .and. cpar%verbosity>1)  write(*,*) '| no data bands available for sampling of spec ind'
       return
    else
       if (myid_pix==0 .and. cpar%verbosity>2) write(*,*) '| ### Using '//trim(c_lnL%pol_lnLtype(p,id))//' lnL evaluation ###'
       if (cpar%verbosity>3 .and. myid_pix == 0) then
          write(*,fmt='(a)') ' |  Active bands'
          do k = 1,band_count
             write(*,fmt='(a,i1)') ' |   band: '//trim(data(band_i(k))%label)//', -- polarization: ',pol_j(k)
          end do

          if (c_lnL%spec_mono_combined(par_id)) then
             write(*,fmt='(a)') ' |  Active monopoles'
             do k = 1,numband
                if (monopole_active(k)) write(*,fmt='(a,i1)') ' |   band: '//trim(data(k)%label)
             end do
          end if
       end if
    end if


    !allocate and assign low resolution reduced data maps
    allocate(reduced_data(0:np_lr-1,band_count))
    do k = 1,band_count
       reduced_data(:,k) = res_smooth(band_i(k))%p%map(:,pol_j(k))
    end do

    ! init full resolution theta maps for smoothing
    theta_single_fr => comm_map(info_fr_single)
    theta_fr => comm_map(info_fr_single)

    !ud_grade mask
    mask_lr => comm_map(info_lr)
    call c_lnL%pol_ind_mask(id)%p%udgrade(mask_lr)
    where (mask_lr%map < 0.5d0)
       mask_lr%map = 0.d0
    elsewhere
       mask_lr%map = 1.d0
    end where

    !init lowres residual map
    res_map => comm_map(info_lr)

    ! Init inverse noise diagonal map
    Ninv_map => comm_map(info_lr)


    ! This is used for marginal/ridge sampling
    allocate(all_thetas(npar))
    allocate(mixing_new_arr(band_count),invN_arr(band_count),data_arr(band_count))

    !This is used for all
    allocate(old_thetas(0:npixreg),new_thetas(0:npixreg),init_thetas(0:npixreg))
    allocate(new_theta_smooth(0:np_lr-1))
    allocate(accept_arr(n_prop_limit), dlnL_arr(n_prop_limit))
    allocate(theta_corr_arr(n_prop_limit))
    
    !if monopole sampling
    if (c_lnL%spec_mono_combined(par_id)) then
       allocate(old_mono(numband),new_mono(numband))
       old_mono=monopole_val
       new_mono=old_mono
    end if

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

    old_thetas = c_lnL%theta_pixreg(:npixreg,p,id)
    old_thetas = min(max(old_thetas,theta_min),theta_max)
    new_thetas = old_thetas
    ! that the root processor operates on
    init_thetas = old_thetas
    if (cpar%verbosity>2 .and. info_fr%myid == 0 .and. npixreg > 0) write(*,fmt='(a, f10.5)') " |  initial (avg) spec. ind. value: ", &
         & sum(init_thetas(1:npixreg))/npixreg

    if (myid_pix==0) then
       N_theta_MC = 1000*n_prop_limit
       do pr = 1, npixreg
          if (c_lnL%nprop_pixreg(pr,p,id) > N_theta_MC ) N_theta_MC = c_lnL%nprop_pixreg(pr,p,id)
       end do
       allocate(theta_MC_arr(N_theta_MC+1,npixreg,2))
       theta_MC_arr = 0.d0
       if (c_lnL%spec_mono_combined(par_id)) then 
          allocate(multipoles_trace(2,0:N_theta_MC+1,numband,0:3))
          multipoles_trace=0.d0
       end if
    end if
    max_prop = 0

    do pr = 1,npixreg

       if (c_lnL%pol_pixreg_type(p,id) == 3) then
          if (c_lnL%fix_pixreg(pr,p,id)) cycle
       end if

       call wall_time(t0)

       nsamp=-1
       n_spec_prop = 0
       n_accept = 0
       n_corr_prop = 0 
       avg_dlnL = 0.d0
       accept_rate=0.d0
       accept_arr(:)=0
       dlnL_arr(:)=0.d0
       running_dlnL=1.d3
       running_accept=0.d0
       theta_corr_arr=0.d0
       running_correlation=-2.d0
       first_nprop = .true.

       old_theta = old_thetas(pr)

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
       if (.false.) then
          !scaling number of pixels in pixel region (from component resolution to smoothing scale resolution)
          prior_rms_scaling=1.d0*pix_count*info_lr%npix/info_fr%npix
       else
          prior_rms_scaling=1.d0 !no prior RMS scaling
       end if


       if (c_lnL%pol_sample_nprop(p,id) .or. c_lnL%pol_sample_proplen(p,id)) then
          pixreg_nprop = 1000*n_prop_limit !should be enough to find proposal/correlation length, if prompted. 
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

       !then we smooth and downgrade to lower resolution
       smooth_scale = c_lnL%smooth_scale(id)
       if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
          if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then 
             !smooth to correct resolution. First smooth at native component resolution, 
             !then downgrade to low resolution Nside
             call smooth_map(info_fr_single, .false., &
                  & c_lnL%B_pp_fr(id)%p%b_l*0.d0+1.d0, theta_fr, &  
                  & c_lnL%B_pp_fr(id)%p%b_l, temp_map)
             theta_lr_hole => comm_map(info_lr_single)
             call temp_map%udgrade(theta_lr_hole)
             call temp_map%dealloc(); deallocate(temp_map)
             nullify(temp_map)

             call smooth_map(info_fr_single, .false., &
                  & c_lnL%B_pp_fr(id)%p%b_l*0.d0+1.d0, theta_single_fr, &  
                  & c_lnL%B_pp_fr(id)%p%b_l, temp_map)
             theta_single_lr => comm_map(info_lr_single)
             call temp_map%udgrade(theta_single_lr)
             call temp_map%dealloc(); deallocate(temp_map)
             nullify(temp_map)
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
       j=-1
       do while (j < pixreg_nprop) !propose new sample n times (for faster mixing in total)
          j = j + 1
          
          nsamp = nsamp + 1
          if (myid_pix == 0) then

             if (first_sample) then
                old_theta = old_thetas(pr)
                proplen=c_lnL%proplen_pixreg(pr,p,id)
                if (cpar%verbosity>2) then                
                   write(*,fmt='(a, i3, a, f10.5)') " |  initial pix.reg. ",pr,"  -- spec. ind. value: ", init_thetas(pr)
                   write(*,fmt='(a, e14.5)') " |  initial proposal length: ",proplen
                end if
                new_thetas = old_thetas
                new_theta = old_theta
             else
                new_theta = old_theta + proplen*rand_gauss(handle)
             end if

          end if

          !broadcast new_theta
          call mpi_bcast(new_theta, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
          new_thetas(pr) = new_theta

          if (cpar%verbosity>3 .and. myid_pix==0 .and. mod(j,out_every)==0) then
             write(*,fmt='(a, i6, a, i3, a, f10.5, a, f10.5)') " |  proposal: ", j," -- Pixreg ", pr, &
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
                write(*,fmt='(a, f10.5, a, f10.5)') " |    Proposed ind outside limits.  min: ", &
                     & theta_min," -- max: ", theta_max
                write(*,fmt='(a, f10.5)') " |    Proposed ind ",new_theta
             end if
             loop_exit = .true.
          else

             !set up the new theta map
             new_theta_smooth(:)=theta_lr_hole%map(:,1)+ new_theta*theta_single_lr%map(:,1)
             ! threshold smoothed map on uniform limits
             new_theta_smooth(:) =min(theta_max,max(theta_min, new_theta_smooth(:))) 


             !###################################################################################################
             ! if monopole sampling together with spectral parameter:
             !
             ! Given the new theta, compute the best fit monopoles (using the old amplitude for the component),
             ! i.e. conditional sampling
             if (c_lnL%spec_mono_combined(par_id)) then
                ! Need to update the reduced data maps for the k's with active monopoles

                do k = 1,band_count
                   if (monopole_active(band_i(k)) .and. pol_j(k)==1) then
                      do pix = 0, np_lr-1
                         do i = 1, npar
                            if (i == id) cycle
                            all_thetas(i) = c_lnL%theta_smooth(i)%p%map(pix,p) 
                         end do

                         all_thetas(id)=new_theta_smooth(pix)

                         mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                              & data(band_i(k))%gain * c_lnL%cg_scale(pol_j(k))

                         reduced_data(pix,k) = res_smooth(band_i(k))%p%map(pix,pol_j(k)) - mixing_new* &
                              & c_lnL%x_smooth%map(pix,pol_j(k))
                      end do
                      !reduced_data is now the original residual map

                      !create reduced data map for monopole estimation from original residual
                      !i.e. a (residual+monopole)-map, add the original (input) monopole to the residual
                      reduced_data(:,k) = reduced_data(:,k) + monopole_mixing(band_i(k))*monopole_val(band_i(k))
                                            
                      if (trim(monocorr_type) == 'monopole+dipole' .or. &
                           & trim(monocorr_type) == 'monopole-dipole') then
                         ! resample the monopole (and dipole) given the new reduced data
                         md_A = 0.d0
                         md_b = 0.d0
                         harmonics2=harmonics*monopole_mixing(band_i(k))
                         do j_md = 0, 3
                            do k_md = 0, 3
                               md_A(j_md,k_md) = sum(harmonics2(:,j_md) * harmonics2(:,k_md)) 
                            end do

                            md_b(j_md) = sum(reduced_data(:,k) * harmonics2(:,j_md)) !is to be set later, this will change
                         end do
                         !we need to run an MPI reduce to get all harmonics for md_A and md_b
                         call mpi_allreduce(MPI_IN_PLACE, md_A, 16, MPI_DOUBLE_PRECISION, MPI_SUM, info_lr%comm, ierr)
                         call mpi_allreduce(MPI_IN_PLACE, md_b, 4, MPI_DOUBLE_PRECISION, MPI_SUM, info_lr%comm, ierr)

                         !solve the mono-/dipole system
                         call solve_system_real(md_A, multipoles, md_b) 

                         if (np_lr > 0) then
                             ones_map%map = 0d0
                             ones_map%map(:,1) = 1d0
                             ones_map%map(:,1) = ones_map%map(:,1) * mask_mono%map(:,1)
                             call rms_smooth(band_i(k))%p%sqrtInvN(ones_map)
                             a = sum((ones_map%map(:,1)/monopole_mixing(band_i(k)))**2)
                         else
                             a = 0d0
                         end if
                         call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, & 
                              & MPI_SUM, info_lr%comm, ierr)

                      else if (trim(monocorr_type) == 'monopole') then
                         if (np_lr > 0) then
                             ones_map%map = 0d0
                             ones_map%map(:,1) = 1d0
                             ones_map%map(:,1) = ones_map%map(:,1) * mask_mono%map(:,1)
                             call rms_smooth(band_i(k))%p%sqrtInvN(ones_map)
                             a = sum((ones_map%map(:,1)/monopole_mixing(band_i(k)))**2)
      
                             ones_map%map = reduced_data
                             call rms_smooth(band_i(k))%p%sqrtInvN(ones_map)
                             b = sum((ones_map%map(:,1)/monopole_mixing(band_i(k)))**2)
                         else
                             a = 0d0
                             b = 0d0
                         end if
                         call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, & 
                              & MPI_SUM, info_lr%comm, ierr)
                         call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, & 
                              & MPI_SUM, info_lr%comm, ierr)


                         if (a > 0.d0) then
                            multipoles(0) = b/a
                         else
                            multipoles(0) = 0.d0
                         end if
                      end if

                      
                      if (info_lr%myid == 0) then

                         if (a > 0.d0) then !we have statistical power to estimate a monopole
                            sigma=sqrt(a)
                            mu = multipoles(0)
                         else if (monopole_rms(band_i(k)) > 0.d0) then
                            mu=0.d0
                            sigma = 0.d0
                         else
                            mu=0.d0
                            sigma = 0.d0
                         end if

                         ! adding prior 
                         if (monopole_rms(band_i(k)) > 0.d0) then
                            sigma_p = 1.d0/monopole_rms(band_i(k))
                            mu_p = monopole_mu(band_i(k))
                            mu = (mu*sigma**2 + mu_p*sigma_p**2)/(sigma**2 + sigma_p**2)
                         end if
                      end if
                
                      ! bcast new monopole to other processors
                      call mpi_bcast(mu, 1, MPI_DOUBLE_PRECISION, 0, info_lr%comm, ierr)
                      new_mono(band_i(k)) = mu

                      if (myid_pix == 0) then
                         !trace multipoles
                         multipoles_trace(1,j,band_i(k),0) = mu
                         if (trim(monocorr_type) == 'monopole+dipole' .or. &
                              & trim(monocorr_type) == 'monopole-dipole') &
                              & multipoles_trace(1,j,band_i(k),1:3) = multipoles(1:3)
                      end if


                      !update the reduced data map with the new monopole
                      reduced_data(:,k) = res_smooth(band_i(k))%p%map(:,pol_j(k)) + &
                           & monopole_mixing(band_i(k))*(monopole_val(band_i(k))-new_mono(band_i(k)))
                      

                   end if !monopole_active(band_i(k) .and. pol_j(k)==1
                end do !band_count

             end if !monopole sampling
             ! now the reduced data sets are updated with new monopoles (for the relevant bands) 
             !###################################################################################################


             !lnL type should split here
             if (trim(c_lnL%pol_lnLtype(p,id))=='chisq') then
                temp_chisq%map = 0d0
                do k = 1,band_count !run over all active bands

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

                      mixing_new = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                           & data(band_i(k))%gain * c_lnL%cg_scale(pol_j(k))

                      temp_chisq%map(pix,pol_j(k)) = reduced_data(pix,k) - mixing_new* &
                           & c_lnL%x_smooth%map(pix,pol_j(k))

                   end do
                   call rms_smooth(band_i(k))%p%sqrtInvN(temp_chisq)
                   lnL_new = lnL_new -0.5d0*sum(temp_chisq%map)**2

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

                   !as we compute the maximum likelihood amplitude in the evaluation, we must split between
                   !polarizations, and compute the evaluations individually for each polarization
                   do l = p_min,p_max
                      l_count=0
                      !build mixing matrix
                      do k = 1,band_count !run over all active bands
                         if (pol_j(k) /= l) cycle
                         l_count = l_count+1
                         mixing_new_arr(l_count) = c_lnL%F_int(pol_j(k),band_i(k),0)%p%eval(all_thetas) * &
                              & data(band_i(k))%gain * c_lnL%cg_scale(pol_j(k))
                         data_arr(l_count)=reduced_data(pix,k)
                         invN_arr(l_count)=rms_smooth(band_i(k))%p%rms_pix(pix, pol_j(k), ret_invN=.true.)
                      end do

                      !compute the marginal/ridge log-likelihood for the pixel, by computing the maximum likelihood chisq
                      lnL_new = lnL_new + comp_lnL_max_chisq_diagonal(mixing_new_arr, invN_arr, data_arr, &
                           & use_det, l_count)
                      
                      !!compute the marginal/ridge log-likelihood for the pixel (This fails in combined monopole sampling!)
                      !lnL_new = lnL_new + comp_lnL_marginal_diagonal(mixing_new_arr, invN_arr, data_arr, &
                      !     & use_det, l_count)
                   end do
                end do
             else 
                write(*,*) 'invalid polarized lnL sampler type: ', p, id, trim(c_lnL%pol_lnLtype(p,id))
                stop

             end if !chisq/marginal/ridge


             !Sum lnL_new (and lnL_old) among all processors

             call mpi_allreduce(MPI_IN_PLACE, lnL_new, 1, MPI_DOUBLE_PRECISION, & 
                  & MPI_SUM, info_fr%comm, ierr)

             !add spec. ind. prior
             if (c_lnL%p_gauss(2,id) > 0.d0) then
                !Find prior "chisq" and add it to lnL
                lnL_prior = (new_thetas(pr)-c_lnL%p_gauss(1,id))**2
                !prior variance is scaled by 
                ! 1/<number of pixels in region (smooth scale resolution)>
                lnL_prior = -0.5d0 * lnl_prior/(c_lnL%p_gauss(2,id)**2 / prior_rms_scaling) 
                lnL_new = lnL_new + lnL_prior
             end if

             !first sample done
             if (first_sample) lnL_init = lnL_new
             if (first_sample .and. myid_pix == 0 .and. cpar%verbosity > 2) write(*,fmt='(a, e14.5)') " |    lnL_init = ", lnL_init

          end if !new_theta outside spec limits

          if (first_sample) then !set lnL_old equal to lnL_new, as no changes to the theta was made (new_theta==old_theta)
             lnL_old = lnL_new
             first_sample = .false.
             if (c_lnL%spec_mono_combined(par_id)) then
                old_mono=new_mono
             end if
             sample_accepted=.true.
          else if (myid_pix == 0) then
             !accept/reject new spec ind
             delta_lnL = lnL_new-lnL_old

             if (cpar%verbosity>3 .and. mod(j,out_every)==0) then
                write(*,fmt='(a, e14.5)') " |    lnL_new = ", lnL_new
                write(*,fmt='(a, e14.5)') " |    lnL_old = ", lnL_old
                write(*,fmt='(a, e14.5)') " |    lnL_new - lnL_old = ", delta_lnL
             end if

             avg_dlnL = (avg_dlnL*(n_spec_prop) + abs(delta_lnL))/(n_spec_prop+1) 

             n_spec_prop = n_spec_prop + 1
             arr_ind = n_spec_prop
             if (arr_ind > n_prop_limit) arr_ind=n_prop_limit
             if (n_spec_prop > n_prop_limit)  then !shift the array indices one step back to allow for progression
                dlnL_arr(:n_prop_limit-1)=dlnL_arr(2:n_prop_limit)
                accept_arr(:n_prop_limit-1)=accept_arr(2:n_prop_limit)
                theta_corr_arr(:n_prop_limit-1)=theta_corr_arr(2:n_prop_limit)
             end if
             dlnL_arr(arr_ind) = abs(delta_lnL)
             running_dlnL = sum(dlnL_arr(1:min(n_spec_prop,n_prop_limit)))/max(min(n_spec_prop,n_prop_limit),1)

             if (delta_lnL < 0.d0) then
                !if delta_lnL is more negative than -25.d0, limit to -25.d0
                if (abs(delta_lnL) > delta_lnL_threshold) delta_lnL = -delta_lnL_threshold 

                if (trim(cpar%operation) == 'sample') then
                   !draw random uniform number
                   a = rand_uni(handle) !draw uniform number from 0 to 1
                   if (exp(delta_lnL) > a) then
                      !accept
                      sample_accepted=.true.
                      old_theta = new_theta
                      lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                      n_accept = n_accept + 1
                      accept_arr(arr_ind) = 1 !accept
                      if (c_lnL%spec_mono_combined(par_id)) then
                         old_mono=new_mono
                      end if
                   else
                      sample_accepted=.false.
                      accept_arr(arr_ind) = 0 !reject
                   end if
                else
                   sample_accepted=.false.
                   accept_arr(arr_ind) = 0 !reject if running optimize
                end if
             else
                sample_accepted=.true.

                !accept new sample, higher likelihood
                old_theta = new_theta
                lnL_old = lnL_new !don't have to calculate this again for the next rounds of sampling
                n_accept = n_accept + 1
                accept_arr(arr_ind) = 1 !accept
                if (c_lnL%spec_mono_combined(par_id)) then
                   old_mono=new_mono
                end if
             end if
             
             if (j <= N_theta_MC) then
                if (j > max_prop) max_prop = j
                theta_MC_arr(j,pr,1) = new_theta
                theta_MC_arr(j,pr,2) = old_theta
                if (c_lnL%spec_mono_combined(par_id)) then
                   if (sample_accepted) then
                      multipoles_trace(2,j,:,:) = multipoles_trace(1,j,:,:) !accept new multipoles for the trace
                   else
                      if (j > 1) multipoles_trace(2,j,:,:) = multipoles_trace(2,j-1,:,:) !copy the previous sample multipoles as the "accepted" ones
                   end if
                end if
             end if
             !compute the running acceptance and correlation coefficient values
             running_accept = (1.d0*sum(accept_arr(1:min(n_spec_prop,n_prop_limit))))/max(min(n_spec_prop,n_prop_limit),1)
             theta_corr_arr(arr_ind) = old_theta !keeping track of accepted thetas
             running_correlation=calc_corr_coeff(theta_corr_arr(:),min(n_spec_prop,n_prop_limit))

             if (j-1 > burn_in) burned_in = .true.
             ! evaluate proposal_length, then correlation length. 
             !This should only be done the first gibbs iteration, if prompted to do so from parameter file.
             accept_rate = n_accept*1.d0/n_spec_prop
             if (c_lnL%pol_sample_proplen(p,id)) then
                if (cpar%verbosity>3 .and. mod(n_spec_prop,out_every)==0) then
                   write(*,fmt='(a, f6.4)') " |   accept rate = ", running_accept
                   write(*,fmt='(a, e14.5)') " |   avg. abs. delta_lnL = ", running_dlnL
                   write(*,fmt='(a, e14.5)') " |   correlation coeff = ", running_correlation
                   write(*,fmt='(a, e14.5)') " |    lnL_new = ", lnL_new
                   write(*,fmt='(a, e14.5)') " |    lnL_old = ", lnL_old
                   write(*,fmt='(a, e14.5)') " |    lnL_new - lnL_old = ", delta_lnL
                end if
                
                if (.not. burned_in) then 
                   !giving som time to let the first steps "burn in" if one starts far away from max lnL
                   n_spec_prop = 0
                   n_accept = 0
                   avg_dlnL = 0.d0
                else if (n_spec_prop > n_prop_limit) then
                   !make sure to go to half the correlation as given by the limit, 
                   !this ensures a good fitting for the proposal length
                   if (abs(running_correlation) > 0.5*correlation_limit) then
                      ! (slowly) drifting to a value/side, not converged
                      if (running_accept > 0.3d0) then
                         !short proposal lengths, drifting slowly, increase prop. length.
                         accept_scale = 1.5d0
                      else if (running_accept < 0.2d0) then
                         !we are taking too long prop. length, reduce step length
                         accept_scale = 0.5d0
                      else
                         !continue drifting with current value
                         accept_scale = 1.d0
                      end if
                   else if (running_accept > 0.7d0) then ! .and. running_dlnL < 10.d0) then
                      !we are at convergence, but sampling with too short prop. length
                      accept_scale = 1.5d0
                   else if (running_accept < 0.3d0) then
                      !we are at convergence, but sampling with long prop. length
                      accept_scale = 0.5d0
                   else 
                      !we have converged on a value
                      accept_scale = -1.d0
                   end if

                   if (accept_scale > 0.d0) then
                      !prop_len only ever used by root processor
                      proplen = proplen*accept_scale !set new proposal length
                      n_accept = 0 !reset with new prop length
                      n_spec_prop = 0 !reset with new prop length
                      avg_dlnL = 0.d0
                      if (cpar%verbosity>3) then
                         write(*,fmt='(a, f6.4)')  " |      accept rate =         ", running_accept
                         write(*,fmt='(a, e14.5)') " |      avg. abs. delta_lnL = ", running_dlnL
                         write(*,fmt='(a, e14.5)') " |      correlation coeff.  = ", running_correlation
                         write(*,fmt='(a, e14.5)') " |      New prop. len. =      ", proplen
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
                if (sampled_proplen .and. first_nprop) then
                   !we have been running the proplen sampling
                   !reset back to initial conditions
                   first_sample=.true.
                   n_accept = 0 
                   nsamp = 0 !reset sample counter
                   n_spec_prop = 0 
                   avg_dlnL = 0.d0
                   n_corr_prop = 0
                   first_nprop = .false. !make sure not to reset again for this pixel region
                   j = -1 !reset while loop counter
                   theta_MC_arr(:,pr,:) = 0.d0
                end if

                n_corr_prop = n_corr_prop + 1
                if (.not. burned_in) then
                   n_corr_prop = 0
                else if (n_corr_prop > n_corr_limit .and. n_spec_prop > n_prop_limit) then
                   ! check if the correlation coeff. is reasonably small,
                   ! but not as small as the correlation when sampling proposal length, in order to ensure convergence.
                   if (abs(running_correlation) < 0.75d0*correlation_limit) then
                      corr_len = n_corr_prop
                      ! set new nprop as minimum twice the number of samples needed to converge + burn_in, 
                      ! but no more/less than absolute limits. 
                      ! This ensures convergence in most cases later, without having to "push back"/add samples.
                      c_lnL%nprop_pixreg(pr,p,id) = 1.d0*min(max(2*corr_len+burn_in,c_lnL%nprop_uni(1,id)),c_lnL%nprop_uni(2,id))
                      loop_exit=.true. !we have found the correlation length
                      c_lnL%pol_sample_nprop(p,id)=.false. !making sure we dont go back into corr. len. sampling
                   end if
                end if
             else
                if (cpar%verbosity>2 .and. mod(j,out_every)==0) then
                   write(*,fmt='(a, f6.4)') " |    accept rate = ", accept_rate                   
                   if (cpar%verbosity>3) write(*,fmt='(a, f6.4)') "|    avg. abs. delta_lnL = ", avg_dlnL
                   if (cpar%verbosity>3) write(*,fmt='(a, f6.4)') "|    correlation coeff.  = ", running_correlation
                end if
             end if

          end if !end if (first_sample) !myid_pix == 0

          ! if j has been reset/changed by master proc (myid_pix==0), then the others need to know. Same with first_sample logical flag
          call mpi_bcast(j, 1, MPI_INTEGER, 0, info_fr%comm, ierr)
          call mpi_bcast(n_spec_prop, 1, MPI_INTEGER, 0, info_fr%comm, ierr)
          call mpi_bcast(first_sample, 1, MPI_LOGICAL, 0, info_fr%comm, ierr)

          !bcast the running correlation coefficient for the pushback check, else the other procs are going to run aditional times and get stuck in some bcats/mpi_allreduce call
          call mpi_bcast(running_correlation, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
          call wall_time(t2)

          if (cpar%verbosity>3 .and. myid_pix == 0 .and. mod(j,out_every)==0 .and. n_spec_prop/=0) then
             write(*,*) '|    Sample:',n_spec_prop,'  Wall time per sample:',real((t2-t1)/n_spec_prop,sp)
          end if

          call mpi_bcast(loop_exit, 1, MPI_LOGICAL, 0, info_fr%comm, ierr)


          if (loop_exit) exit

          if (j == pixreg_nprop) then
             !put in a check if correlation coeff is high. If i is and user wants to make sure of convergence, push back j with n_prop_limit samples
             if (c_lnL%spec_corr_convergence(id)) then
                if (n_spec_prop > c_lnL%nprop_uni(2,id)) then
                   if (cpar%verbosity>3 .and. myid_pix == 0) then
                      write(*,*) '| Maximum number of samples reached. Convercence criteria not reached.'
                      write(*,*) '| Current correlation coefficient: ',running_correlation
                      write(*,*) '| Convergence criteria:            ',correlation_limit
                   end if
                else if (abs(running_correlation) > correlation_limit) then
                   if (cpar%verbosity>3 .and. myid_pix == 0) & 
                        & write(*,*) '| pushing back samples for convergence. Nsamples pushed:', n_prop_limit
                   j = j - n_prop_limit
                   if (j < 0) j = 0
                end if
             end if
          end if
       end do !while j < nprop
       if (pr == 1) lnl_total_init = lnl_init

       call wall_time(t2)

       !save old theta for further pixelregions
       call mpi_bcast(old_theta, 1, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
       old_thetas(pr)=old_theta

       if (c_lnL%spec_mono_combined(par_id)) then
          !save old_mono further pixelregions. As monopoles are resampled from original monopoles and reduced data,
          !this should not be necessary. Also, it is reset at the beginning of each region. We only need to bcast
          !the "old_mono", i.e. last accepted monopoles after all pixelregions are done.
          call mpi_bcast(old_mono, numband, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)
       end if

       if (cpar%verbosity>2) then !pixreg pr
          if (myid_pix==0) then
             write(*,fmt='(a, i5)')    " |    Pixel region: ",pr
             write(*,fmt='(a, f10.5)') " |      Final spec. ind. value:                   ", old_thetas(pr)
             write(*,fmt='(a, e14.5)') " |      Difference in spec. ind., new - old:      ", &
                  & (old_thetas(pr)-init_thetas(pr))
             write(*,*) '|      Samples:',nsamp
             if (nsamp > 0) write(*,*) '|         Wall time per sample (sec):                 ',real((t2-t1)/nsamp,sp)
             write(*,*) '|         Initialization wall time pixel region (sec):',real((t1-t0),sp)
             write(*,*) '|         Total wall time pixel region (sec):         ',real((t2-t0),sp)

             if (sampled_nprop) write(*,fmt='(a, i5)') " |      Number of proposals after tuning: ",c_lnL%nprop_pixreg(pr,p,id)
             if (sampled_proplen) write(*,fmt='(a, e14.5)') " |      Proposal length after tuning: ",c_lnL%proplen_pixreg(pr,p,id)
             write(*,fmt='(a, e14.5)') " |      New Log-Likelihood:                       ", lnl_old
             write(*,fmt='(a, e14.5)') " |      Difference in Log-Likelihood (new - old): ", lnl_old-lnl_init
             write(*,*) '|'
          end if
       end if

       call theta_single_lr%dealloc(); deallocate(theta_single_lr)
       call theta_lr_hole%dealloc(); deallocate(theta_lr_hole)
       theta_single_lr => null()
       theta_lr_hole => null()
       
       !print MC multipoles to file, (partially debug)
       if (.true. .and. c_lnL%spec_mono_combined(par_id) .and. cpar%cs_output_localsamp_maps .and. myid_pix==0) then
          call int2string(iter,         itext)
          call int2string(p,         pind_txt)
          call int2string(cpar%mychain, ctext)
          postfix = 'c'//ctext//'_k'//itext//'_p'//pind_txt
          call int2string(pr,         pind_txt)
          postfix = trim(postfix)//'_pixreg'//pind_txt
          unit = getlun()
          do k_md = 1,numband
             if (.not. monopole_active(k_md)) cycle
             call int2string(k_md,         pind_txt)
             if (trim(monocorr_type) == 'monopole+dipole' .or. &
                  & trim(monocorr_type) == 'monopole-dipole') then
                filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//&
                     & trim(c_lnL%indlabel(id))//'_multipoles_MC_accepted_'//&
                     & trim(postfix)//'_band'//pind_txt//'.dat'
                open(unit,file=trim(filename))

                do i_md = 1,min(j,N_theta_MC)
                   write(unit,fmt='(i7, 4e14.5)') i_md, multipoles_trace(2,i_md,k_md,:)
                end do
                close(unit)
                filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//&
                     & trim(c_lnL%indlabel(id))//'_multipoles_MC_proposed_'//&
                     & trim(postfix)//'_band'//pind_txt//'.dat'
                open(unit,file=trim(filename))

                do i_md = 1,min(j,N_theta_MC)
                   write(unit,fmt='(i7, 4e14.5)') i_md, multipoles_trace(1,i_md,k_md,:)
                end do
                close(unit)
             else
                filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//&
                     & trim(c_lnL%indlabel(id))//'_multipoles_MC_accepted_'//&
                     & trim(postfix)//'_band'//pind_txt//'.dat'
                open(unit,file=trim(filename))

                do i_md = 1,min(j,N_theta_MC)
                   write(unit,fmt='(i7, e14.5)') i_md, multipoles_trace(2,i_md,k_md,0)
                end do
                close(unit)
                filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//&
                     & trim(c_lnL%indlabel(id))//'_multipoles_MC_proposed_'//&
                     & trim(postfix)//'_band'//pind_txt//'.dat'
                open(unit,file=trim(filename))

                do i_md = 1,min(j,N_theta_MC)
                   write(unit,fmt='(i7, e14.5)') i_md, multipoles_trace(1,i_md,k_md,0)
                end do
                close(unit)
             end if
          end do
          multipoles_trace=0.d0
       end if


    end do 
    if (allocated(theta_corr_arr)) deallocate(theta_corr_arr)

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
       write(*,*) "|    Average values from pixel region sampling"
       write(*,fmt='(a, i5)') " |   Number of proposals (after tuning):  ", &
            & sum(c_lnL%nprop_pixreg(1:npixreg,p,id))/npixreg
       write(*,fmt='(a, e14.5)') " | --  Proposal length (after tuning): ", &
            & sum(c_lnL%proplen_pixreg(1:npixreg,p,id))/npixreg
       write(*,fmt='(a, f10.5)') " |   Final spec. ind. value:           ", sum(old_thetas(1:npixreg))/npixreg
       write(*,fmt='(a, e14.5)') " |   Difference in spec. ind., new - old: ", &
            & sum(old_thetas(1:npixreg)-init_thetas(1:npixreg))/npixreg
       write(*,fmt='(a, e14.5)') " |   New Log-Likelihood:               ", &
            & lnl_old
       write(*,fmt='(a, e14.5)') " |   Difference in Log-Likelihood (new - old): ", &
            & lnl_old-lnl_total_init
       write(*,*) '|'
    end if

    call int2string(iter,         itext)
    call int2string(p,         pind_txt)
    call int2string(cpar%mychain, ctext)
    postfix = 'c'//ctext//'_k'//itext//'_p'//pind_txt


    if (iter >= c_lnL%sample_first_niter) then !only stop sampling after n'th iteration, in case first iter is far off.
       c_lnL%pol_sample_nprop(p,id) = .false. !set sampling of correlation length (number of proposals) and
       c_lnL%pol_sample_proplen(p,id) = .false. !proposal length to false (this is only to be done in the first 
                                               !gibbs iteration anyways)
    else !reset sampling to true if we have sampled
       if (sampled_nprop) c_lnL%pol_sample_nprop(p,id) = .true. 
       if (sampled_proplen) c_lnL%pol_sample_proplen(p,id) = .true. 
    end if


    !print MC theta to file, (partially debug)
    if (.true. .and. cpar%cs_output_localsamp_maps .and. myid_pix==0) then
       unit = getlun()
       do pr = 1,npixreg
          call int2string(pr,         pind_txt)
       
          filename=trim(cpar%outdir)//'/'//trim(c_lnl%label)//'_'//trim(c_lnL%indlabel(id))//&
               & '_theta_MC_'//trim(postfix)//'_pixreg'//pind_txt//'.dat'
          open(unit,file=trim(filename))
          do i = 1,min(max_prop,N_theta_MC)
             write(unit,fmt='(i7, 2e14.5)') i,theta_MC_arr(i,pr,:)
          end do
          close(unit)
       end do
       deallocate(theta_MC_arr)
    end if


    !##########################################################################################
    !bcast last valid monopoles to update monopole components
    if (c_lnL%spec_mono_combined(par_id)) then
       call mpi_bcast(old_mono, numband, MPI_DOUBLE_PRECISION, 0, info_fr%comm, ierr)


       if (cpar%verbosity>2 .and. myid_pix==0) then
          !print info on the change in monopoles from initial monopoles 
          write(*,*) '|  Change in monopoles of active bands in sampling. In band units. '
          write(*,*) '|  band_label    mono_in        mono_out       prior_mean     prior_rms'
          do i = 1,numband
             if (monopole_active(i)) write(*,fmt='(a15,e15.5,e15.5,e15.5,e15.5)') &
                  & trim(data(i)%label),monopole_val(i)*monopole_mixing(i), &
                  & old_mono(i)*monopole_mixing(i),monopole_mu(i)*monopole_mixing(i), &
                  & monopole_rms(i)*monopole_mixing(i)
          end do
       end if

       monopole_val=old_mono
       do j = 1,numband
          if (monopole_active(j)) then
             !find the monopole component and update the monopole
             c2           => compList     
             do while (associated(c2))
                select type (c2)
                class is (comm_md_comp)
                   c_mono => c2 !to be able to access all diffuse comp parameters through c_lnL
                   if (trim(c_mono%label) == data(j)%label) then
                      do i_mono = 0, c_mono%x%info%nalm-1
                         call c_mono%x%info%i2lm(i_mono,l_mono,m_mono)
                         if (l_mono == 0) then ! Monopole
                            c_mono%x%alm(i_mono,1) = sqrt(4.d0*pi) * monopole_val(j) 
                         end if
                      end do
                      exit !exit loop
                   end if
                end select

                c2 => c2%nextComp()
             end do
          end if
       end do

    end if
    !##########################################################################################

    ! debug output
    if (cpar%cs_output_localsamp_maps .and. .false.) then
       do i = 1,band_count
          filename=trim(cpar%outdir)//'/'//'reduced_data_band_'//trim(data(band_i(i))%label)//'_'// &
               & trim(c_lnl%label)//'_'//trim(c_lnL%indlabel(id))//'_'//trim(postfix)//'.fits'
          call res_smooth(band_i(i))%p%writeFITS(trim(filename))
       end do

       call update_status(status, "spec_local "//trim(c%label)//' '//trim(c%indlabel(id)) &
            & //' post sampling output, pre residual')

       !create residual maps and print these as well (for the final beta value)
       smooth_scale = c_lnL%smooth_scale(id)
       if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
          if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then 
             !smooth to correct resolution. First smooth at native component resolution, 
             !then downgrade to low resolution Nside
             call smooth_map(info_fr_single, .false., &
                  & c_lnL%B_pp_fr(id)%p%b_l*0.d0+1.d0, theta_fr, &  
                  & c_lnL%B_pp_fr(id)%p%b_l, temp_map)
             theta_lr_hole => comm_map(info_lr_single)
             call temp_map%udgrade(theta_lr_hole)
             call temp_map%dealloc(); deallocate(temp_map)
             nullify(temp_map)

          else !no postproc smoothing, ud_grade to correct resolution
             theta_lr_hole => comm_map(info_lr_single)
             call theta_fr%udgrade(theta_lr_hole)
          end if
       else !native resolution, fr = lr
          theta_lr_hole => comm_map(info_lr_single)
          theta_lr_hole%map(:,1) = theta_fr%map(:,1)
       end if

       do i = 1,band_count
          call int2string(pol_j(i),         pind_txt)

          filename=trim(cpar%outdir)//'/'//'sampled_conditional_residual_band_'//&
               & trim(data(band_i(i))%label)//'_pol'//pind_txt//'_'// &
               & trim(c_lnl%label)//'_'//trim(c_lnL%indlabel(id))//'_'//trim(postfix)//'.fits'
          temp_map => comm_map(info_lr_single)
          temp_map%map(:,1) = res_smooth(band_i(i))%p%map(:,pol_j(i))

          do pix = 0,np_lr-1
             do j = 1, npar
                if (j == id) cycle
                all_thetas(j) = c_lnL%theta_smooth(j)%p%map(pix,p) 
             end do

             all_thetas(id)=theta_lr_hole%map(pix,1)

             mixing_new = c_lnL%F_int(pol_j(i),band_i(i),0)%p%eval(all_thetas) * &
                  & data(band_i(i))%gain * c_lnL%cg_scale(pol_j(i))

             temp_map%map(pix,1) = temp_map%map(pix,1) - &
                  & mixing_new*c_lnL%x_smooth%map(pix,pol_j(i))

          end do

          ! Swap monopoles if monopoles are sampled together with parameter 
          if (c_lnL%spec_mono_combined(par_id)) then
             if (monopole_active(band_i(i)) .and. pol_j(i)==1) then
                temp_map%map(:,1) = temp_map%map(:,1) &
                     & + monopole_val(band_i(i))*monopole_mixing(band_i(i)) & !old monopole
                     & - old_mono(band_i(i))*monopole_mixing(band_i(i)) !new monopole
             end if
          end if

          call temp_map%writeFITS(trim(filename))
          call temp_map%dealloc(); deallocate(temp_map)
          nullify(temp_map)

       end do

       call theta_lr_hole%dealloc(); deallocate(theta_lr_hole)
       nullify(theta_lr_hole)

       call update_status(status, "spec_local "//trim(c%label)//' '//trim(c%indlabel(id)) &
            & //' post sampling output, post residual writing')

    end if


    !deallocate
    deallocate(mixing_new_arr,data_arr,invN_arr,all_thetas)
    deallocate(band_i,pol_j)
    deallocate(old_thetas,new_thetas,init_thetas, new_theta_smooth)
    if (c_lnL%apply_jeffreys) then
       do k = 1, numband
          call df(k)%p%dealloc()
          deallocate(df(k)%p)
       end do
       deallocate(df)
    end if
    deallocate(accept_arr,dlnL_arr)

    call theta_fr%dealloc();        deallocate(theta_fr);        theta_fr => null()
    call theta_single_fr%dealloc(); deallocate(theta_single_fr); theta_single_fr => null()
    call mask_lr%dealloc();         deallocate(mask_lr);         mask_lr => null()
    call res_map%dealloc();         deallocate(res_map);         res_map => null()
    call Ninv_map%dealloc();        deallocate(Ninv_map);        Ninv_map => null()
    
    if (c_lnL%spec_mono_combined(par_id)) then
       call mask_mono%dealloc() 
       mask_mono => null()
       call ones_map%dealloc()
       ones_map => null()
       !deallocate(ones_map,mask_mono)
    end if

    if (allocated(monopole_val)) deallocate(monopole_val)
    if (allocated(monopole_mu)) deallocate(monopole_mu)
    if (allocated(monopole_rms)) deallocate(monopole_rms)
    if (allocated(monopole_mixing)) deallocate(monopole_mixing)
    if (allocated(monopole_active)) deallocate(monopole_active)
    if (allocated(old_mono)) deallocate(old_mono)
    if (allocated(new_mono)) deallocate(new_mono)
    if (allocated(harmonics)) deallocate(harmonics)
    if (allocated(harmonics2)) deallocate(harmonics2)
    if (allocated(multipoles_trace)) deallocate(multipoles_trace)

    call update_status(status, "nonlin pixreg samling end " // trim(c_lnL%label)// ' ' // trim(c_lnL%indlabel(par_id)))

  end subroutine sampleDiffuseSpecIndPixReg_nonlin

  module function comp_lnL_marginal_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_marginal_diagonal

    integer(i4b) :: i, j
    real(dp)     :: MNd,MNM,invMNM
    real(dp), dimension(:), allocatable :: MN
    !  Function to evaluate the log-likelihood for a pixel across all bands in one polarization,
    !  returning the log-likelihood value, equivalent to the highest likelihood chisq for the given theta/mixing matrix.
    !
    !  It does so by substituting the equation for the highest likelihood amplitude into the chisq equation and
    !  making a few assumptions on the behaviour of the theta values.
    !
    !  This evaluation is equivalen to the maximum chisq likelihood evaluation, 
    !  with the exception of a constant difference term.
    !
    !  This version of the evaluation is a little faster than the maximum chisq evaluation, but it does not allow
    !  for combined monopole sampling, as one can show that it will be biased compared to the true value of the
    !  sampled parameter.
    !
    !  Note: This evaluation assumes no correlation between pixels and are only evaluating on a pixel-by-pixel basis.
    !        Also, The length of each input array must be equal, and be larger than or equal to 'arr_len' if the latter
    !        is present.
    !
    !  Arguments:
    !  ------------------
    !  mixing: double precision array, unknown length
    !     An arraycontaining the pixel specific mixing matrix values for the different bands in the evaluation.
    !  invN_arr: double precision array, unknown length
    !     An array with the pixel specific inverse noise variance values for the different bands in the evaluation.
    !  data: double precision array, unknown length
    !     An array with the pixel specific (reduced) data values for the different bands in the evaluation.
    !  use_det: logical
    !     A logical flag specifying to add the logarithm of the determinant of the inverse MNM^T matrix to the
    !     log-likelihood. This is the case for the marginal likelihood evaluation.
    !  arr_len: integer(i4b), optional
    !     If present, one will only evaluate the first arr_len elements in the input arrays.
    !
    !  Returns:
    !  -----------------
    !  comp_lnL_marginal_diagonal: double precision real number
    !     The evaluated log-likelihood value for the pixel, the ridge/marginal likelihood.
    !

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
       invMNM=1.d0/MNM !invert 1x1 matrix
    else
       comp_lnL_marginal_diagonal=-1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    comp_lnL_marginal_diagonal = 0.5d0*MNd*invMNM*MNd

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_marginal_diagonal = comp_lnL_marginal_diagonal - 0.5d0*log(invMNM) 

    deallocate(MN)
  end function comp_lnL_marginal_diagonal

  module function comp_lnL_max_chisq_diagonal(mixing,invN_arr,data,use_det,arr_len)
    implicit none
    logical(lgt),               intent(in)           :: use_det
    real(dp),     dimension(:), intent(in)           :: mixing
    real(dp),     dimension(:), intent(in)           :: invN_arr
    real(dp),     dimension(:), intent(in)           :: data
    integer(i4b),               intent(in), optional :: arr_len
    real(dp)                                         :: comp_lnL_max_chisq_diagonal

    integer(i4b) :: i, j
    real(dp)     :: MNd,MNM,invMNM, amp, chisq
    real(dp), dimension(:), allocatable :: MN
    !  
    !  Evaluates "ridge likelihood" from BP13, arXiv:2201.08188. Take equation
    !  23 from this, but replace the integral with the the maximum likelihood
    !  solution for a given beta.
    !
    !  Implicitly assumes that we are evaluating for a single component.
    !
    !  Function to evaluate the log-likelihood for a pixel across all bands in one polarization,
    !  returning the highest likelihood chisq value of the log-likelihood for the given theta/mixing matrix.
    !
    !  It does so by solving the equation for the highest likelihood amplitude and using this amplitude
    !  to evaluate the chisq.
    !
    !  This evaluation is equivalen to the ridge/marginal evaluation, 
    !  with the exception of a constant difference term.
    !
    !  The benefit of this evaluation is that it allows for combined monopole sampling, 
    !  where one can show that the ridge/marginal evaluation does not return (or are able to find) 
    !  the true value of the sampled parameter, but will be biased. 
    !
    !  Note: This evaluation assumes no correlation between pixels and are only evaluating on a pixel-by-pixel basis.
    !        Also, The length of each input array must be equal, and be larger than or equal to 'arr_len' if the latter
    !        is present.
    !
    !  Arguments:
    !  ------------------
    !  mixing: double precision array, unknown length
    !     An array containing the pixel specific mixing matrix values for the different bands in the evaluation.
    !  invN_arr: double precision array, unknown length
    !     An array with the pixel specific inverse noise variance values for the different bands in the evaluation.
    !  data: double precision array, unknown length
    !     An array with the pixel specific (reduced) data values for the different bands in the evaluation.
    !  use_det: logical
    !     A logicalflag specifying to add the logarithm of the determinant of the inverse MNM^T matrix to the
    !     log-likelihood. This is the case for the marginal likelihood equivalent.
    !  arr_len: integer(i4b), optional
    !     If present, one will only evaluate the first arr_len elements in the input arrays.
    !
    !  Returns:
    !  -----------------
    !  comp_lnL_max_chisq_diagonal: double precision real number
    !     The evaluated log-likelihood value for the pixel, assuming the maximum likelihood chisq given the mixing matrix.
    !
    !  

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

    comp_lnL_max_chisq_diagonal = 0.d0

    if (MNM /= 0.d0) then 
       invMNM=1.d0/MNM !invert 1x1 matrix
       amp = MNd*invMNM
    else
       comp_lnL_max_chisq_diagonal=-1.d30 !MNM = 0.d0, i.e. no inversion possible 
       deallocate(MN)
       return
    end if

    if (present(arr_len)) then
       chisq=sum((data(1:arr_len)-mixing(1:arr_len)*amp)**2 * invN_arr(1:arr_len))
    else
       chisq=sum((data-mixing*amp)**2 * invN_arr)
    end if

    comp_lnL_max_chisq_diagonal = -0.5d0*chisq

    !determinant of 1x1 matrix is the value of the matrix itself
    if (use_det) comp_lnL_max_chisq_diagonal = comp_lnL_max_chisq_diagonal - 0.5d0*log(invMNM) 

    deallocate(MN)
  end function comp_lnL_max_chisq_diagonal


end submodule comm_nonlin_smod
