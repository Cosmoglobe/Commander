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
submodule (comm_diffuse_comp_mod) comm_diffuse_comp_smod
  
contains

  module subroutine initDiffuse(self, cpar, id, id_abs)
    !
    ! Routine that initializes a diffuse type component. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! id: integer
    !       Integer ID of the diffuse component with respect to the active components
    !
    ! id_abs: integer
    !       Integer ID of the diffuse component with respect to all components defined in the parameter file
    !       (and also the id in the 'cpar' parameter)
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    character(len=512) :: filename
    character(len=512) :: temp_filename, temp2
    character(len=512), dimension(1000) :: tokens
    integer(i4b) :: i, j, k, l, m, ntot, nloc, p
    real(dp) :: fwhm_prior, sigma_prior, param_dp
    logical(lgt) :: exist
    type(comm_mapinfo), pointer :: info => null(), info_def => null(), info_ud, info_tempfit
    class(comm_map), pointer :: indmask, mask_ud
    
    call self%initComp(cpar, id, id_abs)

    call update_status(status, "init_diffuse_start")

    ! Initialize variables specific to diffuse source type
    self%pol           = cpar%cs_polarization(id_abs)
    self%nside         = cpar%cs_nside(id_abs)
    self%lmin_amp      = cpar%cs_lmin_amp(id_abs)
    self%lmax_amp      = cpar%cs_lmax_amp(id_abs)
    self%lmax_prior    = cpar%cs_lmax_amp_prior(id_abs)
    self%l_apod        = cpar%cs_l_apod(id_abs)
    self%nu_min        = cpar%cs_nu_min(id_abs)
    self%nu_max        = cpar%cs_nu_max(id_abs)

    if(self%npar == 0) then
       self%lmax_ind = 0 !default
       allocate(self%lmax_ind_mix(3,1))
       self%lmax_ind_mix = 0
    end if

    self%cltype        = cpar%cs_cltype(id_abs)
    self%cg_scale(1:3) = cpar%cs_cg_scale(1:3,id_abs)
    self%nmaps         = 1; if (self%pol) self%nmaps = 3
    self%output_mixmat = cpar%output_mixmat
    self%latmask       = cpar%cs_latmask(id_abs)
    self%apply_jeffreys = .false.
    self%sample_first_niter = cpar%cs_local_burn_in
    self%output_localsamp_maps = cpar%cs_output_localsamp_maps
    self%x_scale       = 1.d0

    only_pol            = cpar%only_pol
    only_I              = cpar%only_I
    output_cg_eigenvals = cpar%output_cg_eigenvals
    outdir              = cpar%outdir
    precond_type        = cpar%cg_precond
    info               => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)
    if (self%pol) then
       self%output_EB     = cpar%cs_output_EB(id_abs)
    else
       self%output_EB     = .false.
    end if
    operation          = cpar%operation


    ! Diffuse preconditioner variables
    npre      = npre+1
    lmax_pre  = max(lmax_pre,  self%lmax_amp)
    nside_pre = min(nside_pre, self%nside)
    nmaps_pre = max(nmaps_pre, self%nmaps)

    ! CMB low-l preconditioner
    if (trim(self%type) == 'cmb' .and. trim(precond_type) == 'diagonal' &
         & .and. self%lmax_pre_lowl > -1) then
       self%lmax_pre_lowl = cpar%cg_lmax_precond
       ntot               = (self%lmax_pre_lowl+1)**2
       nloc               = count(info%lm(1,:) <= self%lmax_pre_lowl)
       if (nloc > 0) allocate(self%invM_lowl(0:ntot-1,0:nloc-1))
    else
       self%lmax_pre_lowl = -1
    end if

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id_abs)) == 'zero' .or. trim(cpar%cs_input_amp(id_abs)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       self%x => comm_map(info, trim(cpar%cs_input_amp(id_abs)))
       do i = 1, self%x%info%nmaps
          self%x%map(:,i) = self%x%map(:,i) / (self%RJ2unit_(i)*self%cg_scale(i))
       end do
       call self%x%YtW

       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l < self%lmin_amp) self%x%alm(i,:) = 0.d0
       end do
    end if
    self%ncr = size(self%x%alm)

    ! Initialize output beam
    self%B_out => comm_B_bl(cpar, self%x%info, 0, 0, fwhm=cpar%cs_fwhm(id_abs), nside=self%nside,&
         & init_realspace=.false.)

    ! Read component mask
    if (trim(cpar%cs_mask(id_abs)) /= 'fullsky' .and. self%latmask < 0.d0) then
       self%mask => comm_map(self%x%info, trim(cpar%cs_mask(id_abs)), &
            & udgrade=.true.)
    end if

    ! Read processing mask
    if (trim(cpar%ds_procmask) /= 'none') then
       self%procmask => comm_map(self%x%info, trim(cpar%ds_procmask), &
            & udgrade=.true.)
    end if

    ! Read index sampling mask; downgrade to each channel; re-use for equal Nsides; skip for dense matrices
    if (trim(cpar%cs_indmask(id_abs)) /= 'fullsky') then
       indmask => comm_map(self%x%info, trim(cpar%cs_indmask(id_abs)), &
            & udgrade=.true.)
       allocate(self%indmask(numband))
       do i = 1, numband
          if (trim(data(i)%N%type) == 'QUmat') cycle
          exist = .false.
          do j = 1, i-1
             if (data(i)%info%nside == data(j)%info%nside .and. &
                  & data(i)%info%nmaps == data(j)%info%nmaps) then
                self%indmask(i)%p => self%indmask(j)%p
                exist = .true.
                exit
             end if
          end do
          if (.not. exist) then
             info_ud => comm_mapinfo(cpar%comm_chain, indmask%info%nside, 0, data(i)%info%nmaps, data(i)%info%nmaps==3)
             mask_ud => comm_map(info_ud)
             do j = 1, min(data(i)%info%nmaps, indmask%info%nmaps)
                mask_ud%map(:,j) = indmask%map(:,j)
             end do
             self%indmask(i)%p => comm_map(data(i)%info)
             call mask_ud%udgrade(self%indmask(i)%p)
             where (self%indmask(i)%p%map > 0.5d0)
                self%indmask(i)%p%map = 1.d0
             elsewhere
                self%indmask(i)%p%map = 0.d0
             end where
             call mask_ud%dealloc(); deallocate(mask_ud)
          end if
       end do
       call indmask%dealloc(); deallocate(indmask)
    end if

    ! Read deflation mask
    if (trim(cpar%cs_defmask(id_abs)) /= 'fullsky') then
       self%lmax_def      = 512
       self%nside_def     = 32
       self%fwhm_def      = 90.d0
       info_def => comm_mapinfo(self%comm, self%nside_def, self%lmax_def, 1, .false.)
       self%defmask => comm_map(info_def, trim(cpar%cs_defmask(id_abs)), udgrade=.true.)
       where (self%defmask%map < 0.5)
          self%defmask%map = 0.d0
       elsewhere
          self%defmask%map = 1.d0
       end where
    else
       self%lmax_def      = -1
       self%nside_def     = 0
       self%fwhm_def      = 0.d0
    end if

    ! Initialize prior mean
    if (trim(cpar%cs_prior_amp(id_abs)) /= 'none') then
       self%mu => comm_map(info, trim(cpar%cs_prior_amp(id_abs)))
       call self%mu%YtW
    end if

    ! Allocate mixing matrix
    call update_status(status, "init_pre_mix")
    self%ndet = maxval(data%ndet)
    allocate(self%F(numband,0:self%ndet), self%F_mean(numband,0:self%ndet,self%nmaps), self%F_null(numband,0:self%ndet))
    self%F_mean = 0.d0
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & self%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
       do j = 0, data(i)%ndet
          if (j<=1) then
            self%F(i,j)%p    => comm_map(info)
            self%F_null(i,j) =  .false.
          else
            do k=1, j
               if (all(data(i)%bp(k)%p%tau0==data(i)%bp(j)%p%tau0)) then
                  self%F(i,j)%p => self%F(i,k)%p
                  self%F_null(i,j) =  .false.
                  exit
               else if (k==j-1) then
                  self%F(i,j)%p    => comm_map(info)
                  self%F_null(i,j) =  .false.
               end if
            end do
          end if
          if (data(i)%bp(j)%p%nu_c < self%nu_min .or. data(i)%bp(j)%p%nu_c > self%nu_max) self%F_null(i,j) =  .true.
       end do
    end do
    call update_status(status, "init_postmix")


    ! Initialize power spectrum
    self%Cl => comm_Cl(cpar, self%x%info, id, id_abs)

    ! Initialize pointers for non-linear search
    if (.not. allocated(res_lowres)) allocate(res_lowres(numband))
    if (.not. allocated(res_smooth)) allocate(res_smooth(numband))
    if (.not. allocated(rms_smooth)) allocate(rms_smooth(numband))

    if (.not. allocated(dust_lowres)) allocate(dust_lowres(numband))
    if (.not. allocated(hotpah_lowres)) allocate(hotpah_lowres(numband))

    ! Set up monopole prior
    self%mono_prior_type = get_token(cpar%cs_mono_prior(id_abs), ":", 1)
    if (trim(self%mono_prior_type) /= 'none') then
       self%cg_samp_group_md = cpar%cg_samp_group_md
       temp_filename = get_token(cpar%cs_mono_prior(id_abs), ":", 2)
       if(temp2(1:1) /= '/') then
          if(trim(self%mono_prior_type) /= 'bandmono') then
            temp_filename = trim(cpar%datadir)// '/' // trim(temp_filename)
          end if
       end if
       call get_tokens(temp_filename, ",", tokens, i)
 
       if (trim(self%mono_prior_type) == 'lower_value_prior') then
          if (i < 4) then
             call report_error('monopole lower_value_prior needs filename,mean,rms,fwhm as input. Not enough inputs found')
          end if
          filename = trim(tokens(1))
          read(tokens(2),*) self%mono_prior_gaussian_mean
          if (self%mono_prior_gaussian_mean < 0.d0) call report_error('Lower value monopole prior requires a non-negative prior mean. component '//trim(self%label))
          read(tokens(3),*) self%mono_prior_gaussian_rms
          if (self%mono_prior_gaussian_rms < 0.d0) call report_error('Lower value monopole prior requires a non-negative prior RMS. component '//trim(self%label))
          read(tokens(4),*) self%mono_prior_fwhm
          if (self%mono_prior_fwhm <= 0.d0) call report_error('Lower value monopole prior requires a positive FWHM. component '//trim(self%label))
          self%mono_prior_map => comm_map(self%x%info, trim(filename))

          !scale prior mean and rms by cg_scale to make sure units match during correction
          self%mono_prior_gaussian_mean = self%mono_prior_gaussian_mean/self%cg_scale(1)
          self%mono_prior_gaussian_rms = self%mono_prior_gaussian_rms/self%cg_scale(1)

          !init smoothing beam for the component amplitude, to be used for monopole prior correction
          self%B_mono_prior => comm_B_bl(cpar, self%x%info, 0, 0, fwhm=self%mono_prior_fwhm, nside=self%nside,&
               & init_realspace=.false.)
       else if (trim(self%mono_prior_type) == 'crosscorr') then
          if (i < 4) then
             call report_error('monopole crosscorr prior needs filename,nside,fwhm,threshold(s) as input. Not enough inputs found. component '//trim(self%label))
          end if
          filename = trim(tokens(1))
          read(tokens(2),*) self%mono_prior_nside
          read(tokens(3),*) self%mono_prior_fwhm
          self%mono_prior_Nthresh=i-3
          allocate(self%mono_prior_threshold(self%mono_prior_Nthresh))
          k = 0
          do j = 4,i
             k=k+1
             read(tokens(j),*) param_dp
             self%mono_prior_threshold(k)=param_dp
          end do

          !read crosscorr map
          info      => comm_mapinfo(cpar%comm_chain, self%mono_prior_nside, &
               & -1, self%x%info%nmaps, self%x%info%pol)
          self%mono_prior_map => comm_map(info, trim(filename))
          !init smoothing beam for the component amplitude, to be used for monopole prior correction
          self%B_mono_prior => comm_B_bl(cpar, self%x%info, 0, 0, fwhm=self%mono_prior_fwhm, nside=self%nside,&
               & init_realspace=.false.)

       else if (trim(self%mono_prior_type) == 'bandmono') then
          filename = get_token(temp_filename, ",", 1)
          self%mono_prior_band=trim(filename)
       else if (trim(self%mono_prior_type) == 'monopole+tempfit') then
          filename = get_token(temp_filename, ",", 1)
          info_tempfit => comm_mapinfo(cpar%comm_chain, self%nside, 0, 2, .false.)
          self%mono_prior_map => comm_map(info_tempfit, trim(filename))
       else          
          filename = get_token(temp_filename, ",", 1)
          self%mono_prior_map => comm_map(self%x%info, trim(filename))
       end if
    end if

    call update_status(status, "init_diffuse_end")
    
  end subroutine initDiffuse

  module subroutine initLmaxSpecind(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, nmaps, p, ierr

    call update_status(status, "lmax_specind_start")

    !check nmaps of the component
    nmaps = 1
    if (cpar%cs_polarization(id_abs)) nmaps=3

    if (self%npar==0) return !do not go further, lmax_ind is set in initDiffuse 

    allocate(self%lmax_ind_pol(3,self%npar))  ! {integer}: lmax per. polarization (poltype index) per spec. ind.

    !self%lmax_ind      = cpar%cs_lmax_ind(id_abs) !not to be used anymore

    self%lmax_ind = -1 !default
    !self%lmax_ind_pol = -1 !default
    do i = 1,self%npar
       do p = 1, self%poltype(i)
          l = cpar%cs_lmax_ind_pol(p,i,id_abs)

          !assign lmax per spec ind per polarization (poltype)
          if (self%poltype(i)==1) then !all polarizations have the same lmax
             self%lmax_ind_pol(:,i) = l 
          else if (self%poltype(i)==2) then 
             if (p == 1) then 
                if (cpar%only_pol) cycle
                self%lmax_ind_pol(1,i) = l
             else
                !if only_pol, we don't want a possible Stokes I to slow down a potential lmax_ind == 0 speed up in SHT
                if (cpar%only_pol) self%lmax_ind_pol(1,i) = l 
                if (p > nmaps) then
                   self%lmax_ind_pol(2:3,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(2:3,i) = l
                end if
             end if
          else if (self%poltype(i)==3) then 
             if (p == 1) then 
                if (cpar%only_pol) cycle
                self%lmax_ind_pol(1,i) = l
             else if (p == 2) then 
                !if only_pol, we don't want a possible Stokes I to slow down a potential lmax_ind == 0 speed up in SHT
                if (cpar%only_pol) self%lmax_ind_pol(1,i) = l 
                if (p > nmaps) then
                   self%lmax_ind_pol(2,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(2,i) = l
                end if
             else
                if (p > nmaps) then
                   self%lmax_ind_pol(3,i) = self%lmax_ind_pol(nmaps,i)
                   cycle
                else
                   self%lmax_ind_pol(3,i) = l
                end if
             end if
          else
             write(*,fmt='(a,i1,a)') 'Incompatible poltype ',p,' for '//trim(self%label)//' '//trim(self%indlabel(i))
             stop
          end if
          if (l > self%lmax_ind) self%lmax_ind = l
       end do
    end do


    allocate(self%lmax_ind_mix(3,self%npar))  ! {integer}: lmax per. polarization (poltype index) per spec. ind.
                                              ! Set to zero if lmax = -1 and fullsky pixregions are used

    self%lmax_ind_mix = self%lmax_ind_pol
    ! Comment: We are defining lmax per poltype index (or polarization) so that only the valid polarizations contribute
    !          to the overall lmax_ind. It also makes sure that if all valid lmax_ind parameters per poltype index and 
    !          parameter is 0, then all non-active are set to zero as well, i.e. we save some time on SHT.
    !          
    ! Note: This subroutine must be called for all diffuse components with npar > 0, before initializing diffuse comp,
    !       i.e. initDiffuse(), but after the number of spec. ind. parameters (npar) and poltype has been set.
    call update_status(status, "lmax_specind_end")

  end subroutine initLmaxSpecind

  module subroutine initPixregSampling(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    integer(i4b) :: i, j, k, l, m, n, p, pr, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    character(len=16),         dimension(1000) :: pixreg_label
    character(len=16),         dimension(1000) :: pixreg_prior
    character(len=512),        dimension(1000) :: tokens
    integer(i4b), allocatable, dimension(:)    :: sum_pix
    real(dp),     allocatable, dimension(:)    :: sum_theta, sum_proplen, sum_nprop
    real(dp),     allocatable, dimension(:,:)  :: m_in, m_out, buffer
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    logical(lgt) :: all_fixed
    class(comm_mapinfo), pointer :: info2 => null(), info_ud => null(), info3 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 

    ! A routine that initializes all parameters needed to perform spectral index sampling for the pixel based
    ! parameters, such as synchrotron beta,dust temperature ect. parameters not included are CO line ratios and
    ! others which use different sampling routines
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! id: integer
    !       Integer ID of the diffuse component with respect to the active components
    !
    ! id_abs: integer
    !       Integer ID of the diffuse component with respect to all components defined in the parameter file
    !       (and also the id in the 'cpar' parameter)
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    !


 
    call update_status(status, "initPixreg_specind_start")
    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_ind, &
         & self%nmaps, self%pol)

    ! Set up smoothing scale information
    allocate(self%smooth_scale(self%npar))
    self%smooth_scale = cpar%cs_smooth_scale(id_abs,1:self%npar)


!!! (local) sampling specific parameters!!!

    allocate(self%pol_lnLtype(3,self%npar))        ! {chisq, ridge, marginal}: evaluation type for lnL
    allocate(self%pol_sample_nprop(3,self%npar))   ! {.true., .false.}: sample nprop on first iteration
    allocate(self%pol_sample_proplen(3,self%npar)) ! {.true., .false.}: sample proplen on first iteration

    !allocating and setting up default correlation limits for sampled spectral parameters, local sampling
    allocate(self%spec_corr_convergence(self%npar),self%spec_corr_limit(self%npar))
    self%spec_corr_convergence(:)=.false. !do not push back (add extra samples) during local sampling by default
    self%spec_corr_limit(:)=0.1d0 !assign a default correlation limit if non is defined 

    allocate(self%pol_pixreg_type(3,self%npar))    ! {0=ignore, 1=fullsky, 2=single_pix, 3=pixel_regions}
    allocate(self%nprop_uni(2,self%npar))          ! {integer}: upper and lower limits on nprop
    allocate(self%npixreg(3,self%npar))            ! {integer}: number of pixel regions per poltye per spec ind
    allocate(self%first_ind_sample(3,self%npar)) !used for pixelregion sampling
    self%first_ind_sample=.true.

    call update_status(status, "initPixreg_specind_pixreg_type")
    
    self%npixreg = 0
    self%pol_pixreg_type = 0
    self%priorsamp_local=.false.
    do i = 1,self%npar       
       do j = 1,self%poltype(i)
          if (j > self%nmaps) cycle
          if (cpar%only_pol .and. j == 1 .and. self%poltype(i) > 1) cycle
          if (self%lmax_ind_pol(j,i) < 0) then
             if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='fullsky') then
                self%pol_pixreg_type(j,i) = 1
                self%npixreg(j,i) = 1
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='single_pix') then
                self%pol_pixreg_type(j,i) = 2
                if (cpar%nside_smooth(smooth_scale) < self%theta(i)%p%info%nside) then
                   self%npixreg(j,i) = 12*(cpar%nside_smooth(self%smooth_scale(i))**2)
                else
                   self%npixreg(j,i) = 12*(self%theta(i)%p%info%nside**2)
                end if
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='pixreg') then
                self%pol_pixreg_type(j,i) = 3
                self%npixreg(j,i) = cpar%cs_spec_npixreg(j,i,id_abs) 
             else
                write(*,*) 'Unspecified pixel region type (', cpar%cs_spec_pixreg(j,i,id_abs), ') for poltype',j,'of spectral index',i,'in component',id_abs
                stop
             end if
          else
             !check if alm using defined pixel regions to suggest theta (then to smooth and get alms from map)
             !if so, pixreg values needs to be added. 
             if (cpar%almsamp_pixreg) then
                self%pol_pixreg_type(j,i) = 3
                self%npixreg(j,i) = cpar%cs_spec_npixreg(j,i,id_abs) 
             end if
          end if
       end do
    end do

!    do i = 1,self%npar       
!       do j = 1,self%poltype(i)
!          if (self%myid == 0) write(*,*) 'i j pol_pixreg_type', i, j, self%pol_pixreg_type(j,i)          
!       end do
!    end do
    
    call update_status(status, "initPixreg_specind_allocate")

    if (any(self%pol_pixreg_type(:,:) > 0)) then
       !find highest number of pixel regions and init some key parameters
       k = 0
       m = 0
       do i = 1,self%npar
          do j = 1,self%poltype(i)
             if (j > self%nmaps) cycle
             if (cpar%only_pol .and. j == 1 .and. self%poltype(i) > 1) cycle

             if (self%pol_pixreg_type(j,i) > 0) then
                if (self%npixreg(j,i) > k) k = self%npixreg(j,i)
             end if
             if (self%pol_pixreg_type(j,i) == 3) then
                if (self%npixreg(j,i) > m) m = self%npixreg(j,i)
             end if
             if (self%lmax_ind_pol(j,i) >= 0) cycle
             self%pol_lnLtype(j,i)  = cpar%cs_spec_lnLtype(j,i,id_abs)
             if (trim(self%pol_lnLtype(j,i)) == 'prior') then
                self%priorsamp_local=.true.
                self%pol_sample_nprop(j,i) = .false.
                self%pol_sample_proplen(j,i) = .false.
             else
                self%pol_sample_nprop(j,i) = cpar%cs_spec_samp_nprop(j,i,id_abs)
                self%pol_sample_proplen(j,i) = cpar%cs_spec_samp_proplen(j,i,id_abs)
             end if
          end do
          if (all(self%lmax_ind_pol(:min(self%nmaps,self%poltype(i)),i) >= 0)) cycle
          self%nprop_uni(:,i)=cpar%cs_spec_uni_nprop(:,i,id_abs)
          self%spec_corr_convergence(i)=cpar%cs_spec_corr_convergence(i,id_abs)
          if (self%spec_corr_convergence(i)) self%spec_corr_limit(i)=abs(cpar%cs_spec_corr_limit(i,id_abs))
       end do

       if (self%myid == 0 .and. cpar%verbosity > 1) write(*,*) '|    Number of pixel regions', k

       allocate(self%nprop_pixreg(k,3,self%npar))
       allocate(self%npix_pixreg(k,3,self%npar))
       allocate(self%proplen_pixreg(k,3,self%npar))
       allocate(self%B_pp_fr(self%npar))
       allocate(self%B_smooth_amp(self%npar))
       allocate(self%B_smooth_specpar(self%npar))
       allocate(self%theta_pixreg(0:k,3,self%npar)) ! extra 0 for fullsky maybe; always nmaps=3...
       allocate(self%theta_pixreg_buff(0:k,3,self%npar))
       allocate(self%prior_pixreg(k,3,self%npar))
       self%theta_pixreg = 1.d0 !just some default values, is set later in the code
       self%theta_pixreg_buff = 1.d0 !just some default values, is set later in the code
       self%nprop_pixreg = 0    ! default values, is set later in the code
       self%proplen_pixreg = 1.d0 ! default values, is set later in the code
!obs read theta from correct init

       if (any(self%pol_pixreg_type(:,:) == 3)) then
          allocate(self%fix_pixreg(m,3,self%npar))
          allocate(self%pixreg_priors(MAXVAL(self%npixreg(:,:)),3,self%npar))
          self%fix_pixreg(:,:,:) = .false.
          do i = 1,self%npar
             self%pixreg_priors(:,:,i) = self%p_gauss(1,i)
             do j = 1,self%poltype(i)
                if (j > self%nmaps) cycle
                if (self%pol_pixreg_type(j,i) == 3) then

                   ! Loop over priors on regions
                   if (.not. trim(cpar%cs_spec_pixreg_priors(j,i,id_abs)) == 'none') then
                      if (self%npixreg(j,i) > 32) write(*,*) "Max pixregs is 20 for this, you're trying",  self%npixreg(j,i)
                      call get_tokens(cpar%cs_spec_pixreg_priors(j,i,id_abs), ",", pixreg_prior, n)
                      if (n == self%npixreg(j,i)) then
                         do pr = 1,self%npixreg(j,i)
                            read(pixreg_prior(pr),*) self%pixreg_priors(pr,j,i)
                         end do
                      else
                         write(*,*) "Must have same number of priors as there are number of regions for par", i, " and poltype", j
                      end if
                   end if
                   !write(*,*) "pixreg priors", self%pixreg_priors(:,j,i)
                   ! Loop over frozen regions
                   if (.not. trim(cpar%cs_spec_fix_pixreg(j,i,id_abs)) == 'none') then
                      do pr = 1,self%npixreg(j,i)
                         call get_tokens(cpar%cs_spec_fix_pixreg(j,i,id_abs), ",", pixreg_label, n)
                         do m = 1, n
                            call str2int(trim(pixreg_label(m)),k,ierr)
                            if (ierr == 0) then
                               if (pr == k) then
                                  self%fix_pixreg(pr,j,i) = .true.
                                  exit
                               end if
                            end if
                         end do
                      end do
                   end if

                end if
             end do
             all_fixed=.true. !assume all pixregs are fixed
             !This should check if all pixregs for all active poltypes are frozen
             ! I.e. if any poltype doen't have pixreg we always allow fixing a full poltype
             ! only if all poltypes have pixreg and everything is frozen we exit
             do j = 1,self%poltype(i)
                if (j > self%nmaps) cycle
                if (self%pol_pixreg_type(j,i) == 3) then
                   if (any(self%fix_pixreg(:self%npixreg(j,i),j,i) .eqv. .false.)) &
                        & all_fixed=.false. !there is an un-frozen pixel region
                else
                   all_fixed=.false. !not all poltypes are pixelregions
                end if
             end do
             if (all_fixed .eqv. .true.) then
                write(*,fmt='(a,a)') 'Component "'//trim(self%label)//'", spec. ind "'&
                     & //trim(self%indlabel(i))//'", all poltypes have pixel region sampling '//&
                     & 'and all regions have been fixed. This only the prior RMS should do. Exiting'
                stop
             end if
          end do !npar
       end if

    end if

    ! Set up spectral index sampling masks, proposal length maps and nprop maps
    allocate(self%pol_ind_mask(self%npar)) ! masks per spectral index (all poltypes)
    allocate(self%pol_nprop(self%npar))    ! nprop map per spectral index (all poltypes
    allocate(self%pol_proplen(self%npar))  ! proplen map per spectral index (all poltypes)
    allocate(self%ind_pixreg_map(self%npar))   ! pixel region map per spectral index (all poltypes)
    allocate(self%spec_mono_combined(self%npar)) !logical array, combined monopole and spectral parameter sampling
    allocate(self%spec_mono_mask(self%npar)) !map pointer array, monopole sampling mask
    allocate(self%spec_mono_type(self%npar)) !map pointer array, monopole sampling mask
    allocate(self%spec_mono_freeze(self%npar)) !map pointer array, monopole sampling mask

    info => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_ind, &
         & self%nmaps, self%pol)

    if (any(self%pol_pixreg_type > 0)) then
       if (self%nmaps > 1 .and. any(self%poltype > 1)) then
          allocate(self%ind_pixreg_arr(0:info%np-1,3,self%npar)) 
       else
          allocate(self%ind_pixreg_arr(0:info%np-1,1,self%npar)) 
       end if
       self%ind_pixreg_arr = 0 !all pixels assigned to pixelregion 0 (not to be sampled), will read in pixreg later
    end if

    if (self%priorsamp_local) allocate(self%theta_prior(2,3,self%npar))
    
    do i = 1,self%npar
       call update_status(status, "initPixreg_spec_monopole_mask")
       self%spec_mono_combined(i)=cpar%cs_spec_mono_combined(id_abs,i)
       if (self%spec_mono_combined(i)) then
          self%spec_mono_type(i)=trim(cpar%cs_spec_mono_type(id_abs,i))
          self%spec_mono_freeze(i)=trim(cpar%cs_spec_mono_freeze(id_abs,i))

          if (self%spec_mono_type(i) == 'monopole' .or. self%spec_mono_type(i) == 'monopole+dipole' .or. &
               & self%spec_mono_type(i) == 'monopole-dipole') then

             info2 => comm_mapinfo(cpar%comm_chain, self%nside, -1, &
                  & self%nmaps, self%pol)
          else
             call report_error('Invalid monopole correction type for the combined monopole sampling of '&
                  & //trim(self%label)//' '//trim(self%indlabel(i))//': '//trim(self%spec_mono_type(i)))
          end if

          if (trim(cpar%cs_spec_mono_mask(id_abs,i)) == 'fullsky') then
             self%spec_mono_mask(i)%p => comm_map(info2)
             self%spec_mono_mask(i)%p%map = 1.d0
          else
             self%spec_mono_mask(i)%p => comm_map(info2, trim(cpar%cs_spec_mono_mask(id_abs,i)))

             if (min(self%poltype(i),self%nmaps) > &
                  & self%spec_mono_mask(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' monopole mask has fewer maps (', & 
                     & self%spec_mono_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                     & ') for component nr. ',id_abs
                call mpi_finalize(ierr)
             else if (self%nside /= self%spec_mono_mask(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//&
                     & ' monopole mask (combined sampler) has different nside (', & 
                     & self%spec_mono_mask(i)%p%info%nside,') than the spectral parameter (', &
                     & info2%nside, ') for component nr. ',id_abs
                call mpi_finalize(ierr)
             end if

          end if
          
          where (self%spec_mono_mask(i)%p%map > 0.5d0)
             self%spec_mono_mask(i)%p%map = 1.d0
          elsewhere
             self%spec_mono_mask(i)%p%map = 0.d0
          end where
       end if
       
       if (any(self%lmax_ind_pol(:min(self%nmaps,self%poltype(i)),i) < 0 .and. &
            & self%pol_pixreg_type(:min(self%nmaps,self%poltype(i)),i) > 0)) then
          call update_status(status, "initPixreg_specind_mask")
          ! spec. ind. mask
          if (trim(cpar%cs_spec_mask(i,id_abs)) == 'fullsky') then
             self%pol_ind_mask(i)%p => comm_map(info)
             self%pol_ind_mask(i)%p%map = 1.d0
          else
             ! Read map from FITS file
             self%pol_ind_mask(i)%p => comm_map(info, trim(cpar%cs_spec_mask(i,id_abs)))

             if (min(self%poltype(i),self%nmaps) > &
                  & self%pol_ind_mask(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' mask has fewer maps (', & 
                     & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                     & ') for component nr. ',id_abs
                stop
             else if (self%theta(i)%p%info%nside /= self%pol_ind_mask(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' mask has different nside (', & 
                     & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                     & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
                stop
             end if
          end if

          call update_status(status, "initPixreg_specind_proplen")

          ! prop. len. map
          if (trim(cpar%cs_spec_proplen(i,id_abs)) == 'fullsky') then
             self%pol_proplen(i)%p => comm_map(info)
             self%pol_proplen(i)%p%map = 1.d0
          else
             ! Read map from FITS file
             self%pol_proplen(i)%p => comm_map(info, trim(cpar%cs_spec_proplen(i,id_abs)))
             if (min(self%poltype(i),self%nmaps) > &
                  & self%pol_proplen(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' proplen map has fewer maps (', & 
                     & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                     & ') for component nr. ',id_abs
                stop
             else if (self%theta(i)%p%info%nside /= self%pol_proplen(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' proplen map has different nside (', & 
                     & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                     & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
                stop
             end if
          end if

          ! replace proplen of input map for given poltype with user specified value, if given
          do j = 1,self%poltype(i)
             if (j > self%nmaps .or. trim(self%pol_lnLtype(j,i)) == 'prior') cycle
             if (cpar%cs_spec_proplen_init(j,i,id_abs) > 0.d0) then
                self%pol_proplen(i)%p%map(:,j)=cpar%cs_spec_proplen_init(j,i,id_abs)
             end if
          end do

          call update_status(status, "initPixreg_specind_nprop")

          ! nprop map
          if (trim(cpar%cs_spec_nprop(i,id_abs)) == 'fullsky') then
             self%pol_nprop(i)%p => comm_map(info)
             self%pol_nprop(i)%p%map = 1.d0
          else
             ! Read map from FITS file
             self%pol_nprop(i)%p => comm_map(info, trim(cpar%cs_spec_nprop(i,id_abs)))
             if (min(self%poltype(i),self%nmaps) > &
                  & self%pol_nprop(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' nprop map has fewer maps (', & 
                     & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                     & ') for component nr. ',id_abs
                stop
             else if (self%theta(i)%p%info%nside /= self%pol_nprop(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' nprop map has different nside (', & 
                     & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                     & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
                stop
             end if
          end if
          ! replace nprop of input map for given poltype with user specified value, if given       
          do j = 1,self%poltype(i)
             if (j > self%nmaps .or. trim(self%pol_lnLtype(j,i)) == 'prior') cycle
             if (cpar%cs_spec_nprop_init(j,i,id_abs) > 0) then
                self%pol_nprop(i)%p%map(:,j)=cpar%cs_spec_nprop_init(j,i,id_abs)*1.d0
             end if
          end do

          ! limit nprop
          self%pol_nprop(i)%p%map = min(max(self%pol_nprop(i)%p%map,self%nprop_uni(1,i)*1.d0), &
               & self%nprop_uni(2,i)*1.d0)

       end if !any lmax_ind_pol < 0

       call update_status(status, "initPixreg_specind_pixel_regions")

       ! initialize pixel regions if relevant 
       if (any(self%pol_pixreg_type(1:self%poltype(i),i) /= 0)) then

          info2  => comm_mapinfo(self%theta(i)%p%info%comm, self%nside, &
               & 3*self%nside, 1, .false.) 

          smooth_scale = self%smooth_scale(i)
          if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
             !create beam for pre-sampling smoothing
             info3  => comm_mapinfo(self%theta(i)%p%info%comm, self%nside, &
                  & self%lmax_amp, self%nmaps, self%pol) 

             self%B_smooth_amp(i)%p => comm_B_bl(cpar, info3, 1, 1, &
                  & fwhm=cpar%fwhm_smooth(smooth_scale), &
                  & nside=self%nside, &
                  & init_realspace=.false.)

             self%B_smooth_specpar(i)%p => comm_B_bl(cpar, info2, 1, 1, &
                  & fwhm=cpar%fwhm_smooth(smooth_scale), &
                  & nside=self%nside, &
                  & init_realspace=.false.)

             !create beam for post processing smoothing
             self%B_pp_fr(i)%p => comm_B_bl(cpar, info2, 1, 1, &
                  & fwhm=cpar%fwhm_postproc_smooth(smooth_scale), &
                  & nside=self%nside, &
                  & init_realspace=.false.)
          end if

          call update_status(status, "initPixreg_specind_pixreg_map")
          if (any(self%pol_pixreg_type(1:self%poltype(i),i) == 3)) then !pixreg map defined from file ! Trygve edited to any
             if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'fullsky') then
                self%ind_pixreg_map(i)%p => comm_map(info)
                self%ind_pixreg_map(i)%p%map = 1.d0
             else if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'none') then
                self%ind_pixreg_map(i)%p => comm_map(info)
                self%ind_pixreg_map(i)%p%map = 0.d0
             else
                ! Read map from FITS file
                self%ind_pixreg_map(i)%p => comm_map(info, trim(cpar%cs_spec_pixreg_map(i,id_abs)))
                if (min(self%poltype(i),self%nmaps) > &
                     & self%ind_pixreg_map(i)%p%info%nmaps) then
                   write(*,fmt='(a,i2,a,i2,a,i2)') trim(self%indlabel(i))//' pixreg map has fewer maps (', & 
                        & self%pol_ind_mask(i)%p%info%nmaps,') than poltype (',self%poltype(i), &
                        & ') for component nr. ',id_abs
                   stop
                else if (self%theta(i)%p%info%nside /= self%ind_pixreg_map(i)%p%info%nside) then
                   write(*,fmt='(a,i4,a,i4,a,i2)') trim(self%indlabel(i))//' pixreg map has different nside (', & 
                        & self%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                        & self%theta(i)%p%info%nside, ') for component nr. ',id_abs
                   stop
                end if
             end if
          else !pixreg map assigned later.
             self%ind_pixreg_map(i)%p => comm_map(info)
             self%ind_pixreg_map(i)%p%map = 1.d0
          end if !any pixregtype == 3

          ! Check if pixel region type = 'fullsky' og 'single_pix' for given poltype index  
          do j = 1,self%poltype(i)
             if (j > self%nmaps) cycle
             !fullsky (or fullsky gauusian prior)
             if (self%pol_pixreg_type(j,i) == 1) then
                self%ind_pixreg_map(i)%p%map(:,j) = 1.d0
             else if (self%pol_pixreg_type(j,i) == 2) then 
                !single pix at smoothing scale pix size
                !make a map at smooth scale nside, give pixels value from 1 to npix
                !udgrade to full resolution
                if (cpar%nside_smooth(smooth_scale) > self%nside) then
                   !only allow nside up to full resolution
                   info_ud  => comm_mapinfo(self%theta(i)%p%info%comm, &
                        & self%theta(i)%p%info%nside, -1, 1, .false.) 
                else
                   info_ud  => comm_mapinfo(self%theta(i)%p%info%comm, &
                        & cpar%nside_smooth(smooth_scale), -1, 1, .false.) 
                end if

                if (cpar%nside_smooth(smooth_scale) < self%theta(i)%p%info%nside) then
                   tp => comm_map(info_ud)
                   tp%map(:,1) = tp%info%pix*1.d0 + 1.d0 !set pixel value equal to pixel number+1 !obs
                   allocate(m_in(0:info_ud%npix-1,1))
                   allocate(m_out(0:self%ind_pixreg_map(i)%p%info%npix-1, 1))
                   allocate(buffer(0:self%ind_pixreg_map(i)%p%info%npix-1, 1))
                   m_in(tp%info%pix,1) = tp%map(:,1)
                   call udgrade_ring(m_in, info_ud%nside, m_out, &
                        & self%ind_pixreg_map(i)%p%info%nside)
                   call mpi_allreduce(m_out, buffer, size(m_out), MPI_DOUBLE_PRECISION, &
                        & MPI_SUM, self%ind_pixreg_map(i)%p%info%comm, ierr)

                   self%ind_pixreg_map(i)%p%map(:,j) = buffer(self%ind_pixreg_map(i)%p%info%pix,1)
                   deallocate(m_in, m_out, buffer)
                   call tp%dealloc(); deallocate(tp)
                   tp => null()
                else
                   self%ind_pixreg_map(i)%p%map(:,j) = self%ind_pixreg_map(i)%p%info%pix*1.d0 +1.d0
                end if
             end if
          end do

          call update_status(status, "initPixreg_specind_pixreg_values")

          ! check if there is a non-smoothed init map for theta of pixelregions
          if (cpar%cs_pixreg_init_theta(i,id_abs) == 'none') then
             tp => comm_map(self%theta(i)%p%info)
             tp%map = self%theta(i)%p%map !take average from existing theta map   !!!what if no map
          else
             !read map from init map (non-smoothed theta map)
             tp => comm_map(self%theta(i)%p%info, trim(cpar%cs_pixreg_init_theta(i,id_abs)))
          end if

          !compute the average theta in each pixel region for the poltype indices that sample theta using pixel regions
          do j = 1,self%poltype(i)
!             if (self%myid == 0) write(*,*) 'A0 j self%poltype(i)          ', j, self%poltype(i) 
             if (j > self%nmaps) cycle
             
!             self%theta_pixreg(:,j,i)=self%p_gauss(1,i) !why is this always set to prior mean and not default!!
             !             self%theta_pixreg(:,j,i)= mean(tp%map(:,j)) ! putting to init/default value instead of prior
!             if (self%myid == 0) write(*,*) 'A1 j self%pol_pixreg_type(j,i)', j, self%pol_pixreg_type(j,i)
             if (self%pol_pixreg_type(j,i) < 1) cycle
!             if (self%myid == 0) write(*,*) 'A2 j                          ', j

             n=self%npixreg(j,i)
             allocate(sum_pix(1:n),sum_theta(0:n),sum_nprop(n),sum_proplen(n))
             sum_theta=0.d0
             sum_pix=0
             self%nprop_pixreg(:,j,i) = 0
             self%proplen_pixreg(:,j,i) = 1.d0
             sum_nprop=0.d0
             sum_proplen=0.d0

             do k = 0,self%theta(i)%p%info%np-1
                do m = 1,n
                   if ( self%ind_pixreg_map(i)%p%map(k,j) > (m-0.5d0) .and. &
                        & self%ind_pixreg_map(i)%p%map(k,j) < (m+0.5d0) ) then
                      !sum_theta(m)=sum_theta(m)+self%theta(i)%p%map(k,j)
                      sum_theta(m)=sum_theta(m)+tp%map(k,j)
                      sum_pix(m)=sum_pix(m)+1 
                      self%ind_pixreg_arr(k,j,i)=m !assign pixel region index 
                      if (self%lmax_ind_pol(j,i) < 0) then
                         sum_proplen(m)=sum_proplen(m)+self%pol_proplen(i)%p%map(k,j)
                         sum_nprop(m)=sum_nprop(m)+self%pol_nprop(i)%p%map(k,j)
                      end if
                      exit
                   end if
                end do
                !region=0 fullsky case
                sum_theta(0) = sum_theta(0)+tp%map(k,j)
             end do
             !allreduce
             call mpi_allreduce(MPI_IN_PLACE, sum_theta, n+1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_pix, n, MPI_INTEGER, MPI_SUM, info%comm, ierr)
             if (self%lmax_ind_pol(j,i) < 0) then
                call mpi_allreduce(MPI_IN_PLACE, sum_proplen, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
                call mpi_allreduce(MPI_IN_PLACE, sum_nprop, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             end if

             self%theta_pixreg(0,j,i)=sum_theta(0)/tp%info%npix ! init at fullsky mean value
             do k = 1,n
                !if (cpar%myid == cpar%root) write(*,*) 'pixreg',k,'  -- number of pixels',sum_pix(k)
                if (sum_pix(k) > 0) then
                   self%theta_pixreg(k,j,i) = sum_theta(k)/(1.d0*sum_pix(k))
                   if (self%lmax_ind_pol(j,i) < 0) then
                      self%nprop_pixreg(k,j,i) = IDINT(sum_nprop(k)/(1.d0*sum_pix(k)))
                      self%proplen_pixreg(k,j,i) = sum_proplen(k)/(1.d0*sum_pix(k))
                   end if
                else
                    self%theta_pixreg(k,j,i) = self%theta_pixreg(0,j,i) ! put to fullsky value
!                   self%theta_pixreg(k,j,i) = self%p_gauss(1,i) ! the prior as theta
                   if (self%lmax_ind_pol(j,i) < 0) then
                      self%nprop_pixreg(k,j,i) = 0
                      self%proplen_pixreg(k,j,i) = 1.d0
                   end if
                end if
                self%npix_pixreg(k,j,i)=sum_pix(k)
             end do
!             self%theta_pixreg(0,j,i)=self%p_gauss(1,i) !all pixels in region 0 has the prior as theta



             

!             if (self%myid == 0) write(*,*) 'B0 npix', sum_pix(0), tp%info%npix
!             if (self%myid == 0) write(*,*) 'C0', i, j, self%theta_pixreg(0,j,i)
!             if (self%myid == 0) write(*,*) 'C1', i, j, self%theta_pixreg(1,j,i)
!             !if (self%myid == 0) write(*,*) 'C2', i, j, self%theta_pixreg(2,j,i)
!             if (self%myid == 0) write(*,*) 'CC self%theta_pixreg(0,1,1)', self%theta_pixreg(0,1,1)


             
             deallocate(sum_pix,sum_theta,sum_proplen,sum_nprop)

             if (trim(self%pol_lnLtype(j,i)) == 'prior') then
                self%theta_prior(:,j,i) = cpar%cs_theta_prior(:,j,i,id_abs)
                self%nprop_pixreg(:,j,i) = 1
                self%pol_nprop(i)%p%map(:,j) = 1
                self%proplen_pixreg(:,j,i) = 1.d0
                self%pol_proplen(i)%p%map(:,j) = 1.d0
             end if

          end do !poltype
          
          call tp%dealloc(); deallocate(tp)

          call update_status(status, "initPixreg_specind_precalc_sampled_theta")

          do j = 1,self%poltype(i)
! if (self%myid == 0) write(*,*) 'i, j poltype nmaps', i, j, self%poltype(i), self%nmaps
             if (j > self%nmaps) cycle
             !Should also assign (and smooth) theta map if RMS /= 0 and we're not initializing from HDF
             if (self%p_gauss(2,i) /= 0.d0 .and. (trim(cpar%init_chain_prefix) == 'none' .or. &
                  & (trim(cpar%init_chain_prefix) /= 'none' .and. trim(self%init_from_HDF) == 'none'))) then 
                !We are not initializing from HDF
                if (self%poltype(i) == 1) then
                   p_min = 1; p_max = self%nmaps
                   if (cpar%only_pol) p_min = 2
                else if (self%poltype(i) == 2) then
                   if (j == 1) then
                      p_min = 1; p_max = 1
                   else
                      p_min = 2; p_max = self%nmaps
                   end if
                else if (self%poltype(i) == 3) then
                   p_min = j
                   p_max = self%nmaps
                else
                   write(*,*) 'Unsupported polarization type'
                   stop
                end if
                if (p_min > p_max) cycle !just a guaranty that we dont smooth for nothing

!                if (self%myid == 0) write(*,*) 'd1', self%theta(i)%p%map(0,1:self%nmaps)
                
                smooth_scale = self%smooth_scale(i)
                if (cpar%num_smooth_scales > 0 .and. smooth_scale >= 0) then

                   !ind. map with 1 map (will be smoothed like zero spin map using the existing code)
                   tp => comm_map(info2)

                   do k = 0,info2%np-1
                      tp%map(k,1) = self%theta_pixreg(self%ind_pixreg_arr(k,j,i),j,i) !theta_pixreg!!!
                   end do
                   if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                      !smooth single map as intensity (i.e. zero spin)
                      call smooth_map(info2, .false., &
                           & self%B_pp_fr(i)%p%b_l*0.d0+1.d0, tp, &  
                           & self%B_pp_fr(i)%p%b_l, tp_smooth)

                      tp%map=tp_smooth%map
                      call tp_smooth%dealloc(); deallocate(tp_smooth)
                   end if

!                   if (self%myid == 0) write(*,*) 'd2', self%theta(i)%p%map(0,1:self%nmaps)
!                   if (self%myid == 0) write(*,*) 'd2- tp%map(0,1)', tp%map(0,1), self%ind_pixreg_arr(0,j,i)

                   do k = p_min,p_max
                      self%theta(i)%p%map(:,k) = tp%map(:,1)
                   end do
!                   if (self%myid == 0) write(*,*) 'd3', self%theta(i)%p%map(0,1:self%nmaps)

!                   if (self%myid == 0) write(*,*) 'd4 pmin pmak', p_min, p_max
!                   if (self%myid == 0) write(*,*) 'd5 self%theta_pixreg(0,1,1)', self%theta_pixreg(0,1,1)
!                   if (self%myid == 0) write(*,*) 'd6'


                   
                   if (self%lmax_ind_pol(j,i) >= 0) then
                      info3 => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_ind_pol(j,i), &
                           & 1, .false.)
                      tp_smooth => comm_map(info3)
                      tp_smooth%map = tp%map
                      call tp_smooth%YtW_scalar()
                      do k = p_min,p_max
                         self%theta(i)%p%alm(0:info3%nalm-1,k) = tp_smooth%alm(0:info3%nalm-1,1)
                      end do
                      call tp_smooth%dealloc(); deallocate(tp_smooth)
                   end if
                   call tp%dealloc(); deallocate(tp)
                else
                   if (cpar%num_smooth_scales <= 0) then
                      write(*,*) 'need to define smoothing scales'
                      stop
                   else if (smooth_scale < 0) then
                      write(*,*) 'need to define smoothing scale for component '//&
                           & trim(self%label)//', parameter '//trim(self%indlabel(i))
                      stop
                   end if
                end if !num smooth scales > 0
             end if
          end do !poltype

       end if !any pixreg_type > 0

       call update_status(status, "initPixreg_specind_lmax_ind_mix")

       !final check to see if lmax < 0 and pixel region is fullsky
       do j = 1,self%poltype(i)
          if (j > self%nmaps) then
             self%lmax_ind_mix(j:,i) = self%lmax_ind_mix(self%nmaps,i)
             cycle
          end if

          if (self%poltype(i) == 1) then
             p_min = 1; p_max = 3
          else if (self%poltype(i) == 2) then
             if (j == 1) then
                if (cpar%only_pol) cycle
                p_min = 1; p_max = 1
             else
                p_min = 2; p_max = 3
                if (cpar%only_pol) p_min = 1 !add the first index to this case
             end if
          else if (self%poltype(i) == 3) then
             if (cpar%only_pol .and. j == 1) cycle
             p_min = j
             p_max = j
             if (cpar%only_pol .and. j == 2) p_min = 1 !add the first index to this case
          else
             write(*,*) 'Unsupported polarization type'
             stop
          end if

          if (self%lmax_ind_pol(j,i) >= 0) then
             self%lmax_ind_mix(p_min:p_max,i) = self%lmax_ind_pol(j,i) !in case only_pol and poltype = 2 has lmax > 0
          else if (self%pol_pixreg_type(j,i)==1) then !pixel region is defined fullsky
             self%lmax_ind_mix(p_min:p_max,i) = 0
          else
             self%lmax_ind_mix(p_min:p_max,i) = self%lmax_ind_pol(j,i)
          end if
       end do

    end do !npar
    
    if (self%theta(1)%p%info%myid == 0 .and. .false.) then !debug

       write(*,*) ''
       write(*,*) 'lmax_ind_pol'
       do i = 1,self%npar
          write(*,*) i, self%lmax_ind_pol(:,i)
       end do
       write(*,*) ''
       write(*,*) 'lmax_ind_mix'
       do i = 1,self%npar
          write(*,*) i, self%lmax_ind_mix(:,i)
       end do
       write(*,*) ''
       if (all(self%lmax_ind_mix(:,:) == 0)) write(*,*) '  all lmax_ind_mix zero'

    end if

    call update_status(status, "initPixreg_specind_end")

  end subroutine initPixregSampling

  module subroutine initSpecindProp(self,cpar, id, id_abs)
    !
    !  Subroutine that initializes the spectral index proposal matrix
    !
    !  Arguments:
    !  ----------
    !  self: comm_diffuse_component
    !       Object that contains all of the diffuse component properties
    !  cpar: comm_params
    !       Object that contains parameters input to Commander from parameter
    !       file.
    !  id: int
    !       index of the components, out of only the ones being used
    !  id_abs: int
    !       index of the component, out of all, whether they are used or not. 
    !  Returns:
    !  --------
    !  None, but modifies self by changing the proposal matrices.
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    character(len=512) :: filename
    integer(i4b) :: i, j, l, m, ntot, nloc, p, q
    real(dp) :: fwhm_prior, sigma_prior
    logical(lgt) :: exist
    real(dp),     allocatable, dimension(:)   :: L_arr
    integer(i4b), allocatable, dimension(:)   :: corrlen_arr     
    
    ! Init alm sampling params (Trygve)
    allocate(self%corrlen(self%npar, self%nmaps))
    self%corrlen    = 0     ! Init correlation length
    allocate(corrlen_arr(self%nmaps))
    corrlen_arr     = 0     ! Init correlation length
    
    ! Init bool for L-flags
    allocate(self%L_read(self%npar), self%L_calculated(self%npar))
    self%L_read = .false.
    self%L_calculated = .false.

    ! save minimum chisq per iteration
    allocate(self%chisq_min(self%npar, self%nmaps))

    self%nalm_tot = (self%lmax_ind + 1)**2

    ! Init smooth prior
    allocate(self%sigma_priors(0:self%nalm_tot-1,self%npar)) !a_00 is given by different one

    fwhm_prior = cpar%almsamp_prior_fwhm   !600.d0 ! 1200.d0
    do j = 1, self%npar
       self%sigma_priors(0,j) = self%p_gauss(2,j)
       if (self%nalm_tot > 1) then
          ! Saving and smoothing priors
          i = 1
          do l = 1, self%lmax_ind ! Skip a_00 - m-major ordering (0,0)(1,0)(2,0)(1,1)(1,-1)(2,1)(2,-1)(2,2)(2,-2)
             self%sigma_priors(i,j) = 0.01*self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
             i = i + 1
          end do
          do l = 1, self%lmax_ind
             do m = 1, l
                self%sigma_priors(i,j) = self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                i = i + 1
                self%sigma_priors(i,j) = self%sigma_priors(0,j)*exp(-0.5d0*l*(l+1)*(fwhm_prior * pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
                i = i + 1
             end do
          end do
       end if
    end do

    ! Initialize cholesky matrix
    if (cpar%almsamp_pixreg)  then
       allocate(self%L(0:maxval(self%npixreg), 0:maxval(self%npixreg), self%nmaps, self%npar)) ! Set arbitrary number of max regions
       allocate(L_arr(0:maxval(self%npixreg))) ! Set arbitrary number of max regions
    else
       allocate(self%L(0:self%nalm_tot-1, 0:self%nalm_tot-1, self%nmaps, self%npar))
       allocate(L_arr(0:self%nalm_tot-1))
    end if

    allocate(self%steplen(self%nmaps, self%npar)) !a_00 is given by different one              
    self%L = 0.d0
    self%steplen = 1 !0.3d0

    ! Filename formatting
    do j = 1, self%npar
       if (all(self%lmax_ind_pol(:min(self%nmaps,self%poltype(j)),j)<0)) cycle 
       if (trim(cpar%cs_almsamp_init(j,id_abs)) == 'none') then ! If present cholesky file
          if (cpar%almsamp_pixreg) then
             do p = 1, maxval(self%npixreg)
                self%L(p,p,:,j) = self%sigma_priors(0,j)
             end do
          else
             do p = 0, self%nalm_tot-1
                self%L(p,p,:,j) = self%sigma_priors(p,j)
             end do
          end if
       else
          self%L_read(j) = .true.
          if ( self%myid == 0 ) write(*,*) "|    Initializing alm tuning from ", trim(cpar%cs_almsamp_init(j,id_abs))
          open(unit=11, file=trim(cpar%cs_almsamp_init(j,id_abs)), recl=10000)
          read(11,*) corrlen_arr
          self%corrlen(j,:) = corrlen_arr

          do p = 1, self%nmaps
             read(11,*)
             do q = 0, size(self%L(:,1,p,j))-1
                read(11,*) L_arr
                self%L(q,:,p,j) = L_arr
             end do
          end do
          close(11)
       end if
    end do

    !if ( self%myid == 0 ) write(*,*) "finished Initializing alm tuning"
    
    deallocate(corrlen_arr, L_arr)
    
  end subroutine initSpecindProp


  module subroutine initDiffPrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    select case (trim(precond_type))
    case ("diagonal")
       call initDiffPrecond_diagonal(comm)
    case ("pseudoinv")
       call initDiffPrecond_pseudoinv(comm)
    case default
       call report_error("Preconditioner type not supported: "//trim(precond_type))
    end select

  end subroutine initDiffPrecond

  module subroutine initDiffPrecond_diagonal(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: p1 => null(), p2 => null()
    real(dp),     allocatable, dimension(:,:) :: mat

    call update_status(status, "init_diffpre1")
    if (npre == 0) return
    if (allocated(P_cr%invM_diff)) return
    
    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%nextComp()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if
    call update_status(status, "init_diffpre2")
    
    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(P_cr%invM_diff(0:info_pre%nalm-1,info_pre%nmaps))
    !!$OMP PARALLEL PRIVATE(mat, ind, j, i1, l, m, q, i2, k1, p1, k2, n)
    allocate(mat(npre,npre), ind(npre))
    do j = 1, info_pre%nmaps
       !!$OMP DO SCHEDULE(guided)
       do i1 = 0, info_pre%nalm-1
          call info_pre%i2lm(i1, l, m)
          mat = 0.d0
          do q = 1, numband
             call data(q)%info%lm2i(l,m,i2)
             if (i2 == -1) cycle
             if (j > data(q)%info%nmaps) cycle
             do k1 = 1, npre
                p1 => diffComps(k1)%p
                if (l > p1%lmax_amp) cycle
                if (j > p1%nmaps) cycle
                do k2 = 1, npre
                   p2 => diffComps(k2)%p
                   if (l > p2%lmax_amp) cycle
                   if (j > p2%nmaps) cycle
                   mat(k1,k2) = mat(k1,k2) + &
                        & data(q)%N%invN_diag%alm(i2,j) * &  ! invN_{lm,lm}
                        & data(q)%B(0)%p%b_l(l,j)**2 * &          ! b_l^2
                        & p1%F_mean(q,0,j) * p2%F_mean(q,0,j)    ! F(c1)*F(c2)
                end do
             end do
          end do

          n = 0
          allocate(P_cr%invM_diff(i1,j)%comp2ind(npre))
          P_cr%invM_diff(i1,j)%comp2ind = -1
          do k1 = 1, npre
             if (mat(k1,k1) > 0.d0) then
                n = n+1
                ind(n) = k1
                P_cr%invM_diff(i1,j)%comp2ind(k1) = n
             end if
          end do
          P_cr%invM_diff(i1,j)%n = n
          allocate(P_cr%invM_diff(i1,j)%ind(n))
          allocate(P_cr%invM_diff(i1,j)%M0(n,n), P_cr%invM_diff(i1,j)%M(n,n))
          P_cr%invM_diff(i1,j)%ind = ind(1:n)
          P_cr%invM_diff(i1,j)%M0   = mat(ind(1:n),ind(1:n))
       end do
       !!$OMP END DO
    end do
    deallocate(ind, mat)
    !!$OMP END PARALLEL
    call wall_time(t2)
    call update_status(status, "init_diffpre3")

  end subroutine initDiffPrecond_diagonal


  module subroutine initDiffPrecond_pseudoinv(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, k1, k2, q, l, m, n
    real(dp)     :: t1, t2
    integer(i4b), allocatable, dimension(:) :: ind
    class(comm_comp),         pointer :: c => null()
    class(comm_diffuse_comp), pointer :: p1 => null(), p2 => null()
    real(dp),     allocatable, dimension(:,:) :: mat

    if (npre == 0) return
    if (allocated(P_cr%invM_diff)) return
    
    if (.not. allocated(diffComps)) then
       ! Set up an array of all the diffuse components
       allocate(diffComps(npre))
       c => compList
       i =  1
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             diffComps(i)%p => c
             i              =  i+1
          end select
          c => c%nextComp()
       end do
       info_pre => comm_mapinfo(comm, nside_pre, lmax_pre, nmaps_pre, nmaps_pre==3)
    end if
    
    ! Allocate space for pseudo-inverse of U
    call wall_time(t1)
    allocate(P_cr%invM_diff(0:lmax_pre,info_pre%nmaps))
    do j = 1, info_pre%nmaps
       do l = 0, lmax_pre
          allocate(P_cr%invM_diff(l,j)%M(npre,numband+npre))
       end do
    end do

  end subroutine initDiffPrecond_pseudoinv


  module subroutine updateDiffPrecond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    select case (trim(precond_type))
    case ("diagonal")
       call updateDiffPrecond_diagonal(samp_group, force_update)
    case ("pseudoinv")
       call updateDiffPrecond_pseudoinv(samp_group, force_update)
    case default
       call report_error("Preconditioner type not supported: "//trim(precond_type))
    end select

  end subroutine updateDiffPrecond

  module subroutine updateDiffPrecond_diagonal(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, k1, k2, l, m, n, ierr, unit, p, q
    real(dp)     :: W_ref, t1, t2
    logical(lgt), save :: first_call = .true.
    integer(i4b), allocatable, dimension(:) :: ind
    real(dp),     allocatable, dimension(:) :: W
    real(dp),     allocatable, dimension(:) :: cond, W_min
    real(dp),     allocatable, dimension(:,:) :: alm

    if (npre == 0) return
    if (.not. recompute_diffuse_precond .and. .not. force_update) return
    
    ! Initialize current preconditioner to F^t * B^t * invN * B * F
    !self%invM    = self%invM0
    do j = 1, info_pre%nmaps
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n > 0) P_cr%invM_diff(i,j)%M = P_cr%invM_diff(i,j)%M0
       end do
    end do

    ! Right-multiply with sqrt(Cl)
    call wall_time(t1)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       if (.not. diffComps(k1)%p%active_samp_group(samp_group)) cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          if (.not. diffComps(k2)%p%active_samp_group(samp_group)) cycle
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = P_cr%invM_diff(i,j)%M(q,p)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
          !if (info_pre%myid == 0) write(*,*) 'a', k1, k2, alm(4,1)
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre, diag=.true.)
          !if (info_pre%myid == 0) write(*,*) 'b', k1, k2, alm(4,1)
          !call mpi_finalize(j)
          !stop
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) P_cr%invM_diff(i,j)%M(q,p) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do

    ! Left-multiply with sqrt(Cl)
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       if (.not. diffComps(k1)%p%active_samp_group(samp_group)) cycle
       !$OMP PARALLEL PRIVATE(alm, k2, j, i, p, q)
       allocate(alm(0:info_pre%nalm-1,info_pre%nmaps))
       !$OMP DO SCHEDULE(guided)
       do k2 = 1, npre
          if (.not. diffComps(k2)%p%active_samp_group(samp_group)) cycle
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) then
                   alm(i,j) = P_cr%invM_diff(i,j)%M(p,q)
                else
                   alm(i,j) = 0.d0
                end if
             end do
          end do
!          if (info_pre%myid == 0) write(*,*) 'c', k1, k2, alm(4,1)
          call diffComps(k1)%p%Cl%sqrtS(alm=alm, info=info_pre, diag=.true.)
!          if (info_pre%myid == 0) write(*,*) 'd', k1, k2, alm(4,1)
          do j = 1, info_pre%nmaps
             do i = 0, info_pre%nalm-1
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                q = P_cr%invM_diff(i,j)%comp2ind(k2)
                if (p /= -1 .and. q /= -1) P_cr%invM_diff(i,j)%M(p,q) = alm(i,j)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(alm)
       !$OMP END PARALLEL
    end do
    !call wall_time(t2)
    !write(*,*) 'sqrtS = ', t2-t1

    ! Nullify temperature block if only polarization
    if (only_pol) then
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,1)%n == 0) cycle                
          P_cr%invM_diff(i,1)%M = 0.d0
       end do
    end if


    ! Add unity 
    do k1 = 1, npre
       if (trim(diffComps(k1)%p%cltype) == 'none') cycle
       if (.not. diffComps(k1)%p%active_samp_group(samp_group)) cycle
       !!$OMP PARALLEL PRIVATE(i,l,m,j,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, info_pre%nalm-1
          call info_pre%i2lm(i, l, m)
          !if (info_pre%myid == 1) then
          !write(*,*) info_pre%myid, i, l, m
          !end if
          if (l <= diffComps(k1)%p%lmax_amp) then
             do j = 1, info_pre%nmaps
                if (P_cr%invM_diff(i,j)%n == 0) cycle                
                p = P_cr%invM_diff(i,j)%comp2ind(k1)
                if (p > 0) P_cr%invM_diff(i,j)%M(p,p) = P_cr%invM_diff(i,j)%M(p,p) + 1.d0
             end do
          end if
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
    end do


    ! Nullify elements that are not involved in current sample group
    do j = 1, info_pre%nmaps
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n == 0) cycle                
          do k1 = 1, npre
             if (diffComps(k1)%p%active_samp_group(samp_group)) cycle
             p = P_cr%invM_diff(i,j)%comp2ind(k1)
             if (p == -1) cycle
             P_cr%invM_diff(i,j)%M(p,:) = 0.d0
             P_cr%invM_diff(i,j)%M(:,p) = 0.d0
          end do
       end do
    end do

    ! Print out worst condition number in first call
    call wall_time(t1)
    if (first_call .and. output_cg_eigenvals) then
       allocate(cond(0:lmax_pre), W_min(0:lmax_pre))
       cond  = 0.d0
       W_min = 1.d30
       do j = 1, nmaps_pre
          do i = 0, info_pre%nalm-1
             call info_pre%i2lm(i, l, m)
             n = P_cr%invM_diff(i,j)%n
             if (n > 0) then
                allocate(W(n), ind(n))
                call get_eigenvalues(P_cr%invM_diff(i,j)%M, W)
                W_ref    = minval(abs(W))
                cond(l)  = max(cond(l), maxval(abs(W/W_ref)))
                W_min(l) = min(W_min(l), minval(W))
                deallocate(W, ind)
             end if
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, cond,  lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MAX, info_pre%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, W_min, lmax_pre+1, MPI_DOUBLE_PRECISION, MPI_MIN, info_pre%comm, ierr)
       if (info_pre%myid == 0) then
          write(*,*) 'Precond -- largest condition number = ', real(maxval(cond),sp)
          write(*,*) 'Precond -- smallest eigenvalue      = ', real(minval(W_min),sp)
          unit = getlun()
          open(unit, file=trim(outdir)//'/precond_eigenvals.dat', recl=1024)
          do l = 0, lmax_pre
             !write(unit,fmt='(i8,2e16.8)') l, cond(l), W_min(l)
             write(unit,*) l, cond(l), W_min(l)
          end do
          close(unit)
       end if

       first_call = .false.
       deallocate(cond, W_min)
    end if
    call wall_time(t2)
    !write(*,*) 'eigen = ', t2-t1

    ! Invert matrix
    call wall_time(t1)
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n > 0) then
             if (any(P_cr%invM_diff(i,j)%M /= 0.d0)) call invert_matrix_with_mask(P_cr%invM_diff(i,j)%M)
          end if
       end do
    end do
    call wall_time(t2)
    !write(*,*) 'invert = ', t2-t1

    ! Disable preconditioner update
    recompute_diffuse_precond = .false.
       
  end subroutine updateDiffPrecond_diagonal


  module subroutine updateDiffPrecond_pseudoinv(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

    integer(i4b) :: i, j, k, l, n, ierr, p, q
    real(dp)     :: t1, t2, Cl
    real(dp),     allocatable, dimension(:,:) :: mat

    if (npre == 0) return
    if (.not. recompute_diffuse_precond .and. .not. force_update) return

    ! Build pinv_U
    call wall_time(t1)
    allocate(mat(numband+npre,npre))
    do j = 1, info_pre%nmaps
       do l = 0, lmax_pre

          ! Data section
          mat = 0.d0
          do q = 1, numband
             if (l > data(q)%info%lmax) cycle
             if (j > data(q)%info%nmaps) cycle
             do k = 1, npre
                if (l > diffComps(k)%p%lmax_amp) cycle
                if (j > diffComps(k)%p%nmaps) cycle
                if (.not. diffComps(k)%p%active_samp_group(samp_group)) cycle
                mat(q,k) = data(q)%N%alpha_nu(j) * &
                          & data(q)%B(0)%p%b_l(l,j) * &
                          & diffComps(k)%p%F_mean(q,0,j)

                if (trim(diffComps(k)%p%Cl%type) /= 'none') then
                   mat(q,k) = mat(q,k)*sqrt(diffComps(k)%p%Cl%getCl(l,j))
                end if
                if (mat(q,k) /= mat(q,k)) then
                   if (trim(diffComps(k)%p%Cl%type) /= 'none') then
                      write(*,*) q, j, l, k, data(q)%N%alpha_nu(j), data(q)%B(0)%p%b_l(l,j), &
                           & diffComps(k)%p%F_mean(q,0,j),diffComps(k)%p%Cl%getCl(l,j)
                   else
                      write(*,*) q, j, l, k, data(q)%N%alpha_nu(j), data(q)%B(0)%p%b_l(l,j),&
                           & diffComps(k)%p%F_mean(q,0,j)
                   end if

                end if
             end do
          end do

          ! Prior section
          do k = 1, npre
             if (trim(diffComps(k)%p%Cl%type) == 'none' .or. l > diffComps(k)%p%lmax_amp .or. &
                  & .not. diffComps(k)%p%active_samp_group(samp_group)) cycle
             mat(numband+k,k) = 1.d0
          end do

          ! Store pseudo-inverse of U
          call compute_pseudo_inverse(mat, P_cr%invM_diff(l,j)%M)
          !P_cr%invM_diff(l,j)%M = transpose(mat)
       end do
    end do
    deallocate(mat)
    call wall_time(t2)

    ! Set up array of active components
    if (allocated(ind_pre)) deallocate(ind_pre)
    n = 0
    do k = 1, npre
       if (diffComps(k)%p%active_samp_group(samp_group)) n = n+1
    end do
    if (n > 0) then
       allocate(ind_pre(n))
       i = 0
       do k = 1, npre
          if (diffComps(k)%p%active_samp_group(samp_group)) then
             i = i+1
             ind_pre(i) = k
          end if
       end do
    end if

    ! Disable preconditioner update
    recompute_diffuse_precond = .false.

  end subroutine updateDiffPrecond_pseudoinv
    
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  module subroutine updateDiffuseMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta
    real(dp),  dimension(:,:,:),               intent(in),    optional :: beta  ! Not used here
    integer(i4b),                              intent(in),    optional :: band
    class(map_ptr), dimension(:),              intent(inout), optional :: df    ! Derivative of mixmat with respect to parameter par; for Jeffreys prior
    integer(i4b),                              intent(in),    optional :: par   ! Parameter ID for derivative

    integer(i4b) :: i, j, k, l, n, p, p_min, p_max, nmaps, ierr
    real(dp)     :: lat, lon, t1, t2
    logical(lgt) :: precomp, mixmatnull, bad ! NEW
    character(len=2) :: ctext
    real(dp),        allocatable, dimension(:,:,:) :: theta_p
    real(dp),        allocatable, dimension(:)     :: nu, s, buffer, buff2
    class(comm_mapinfo),          pointer          :: info => null(), info_tp => null()
    class(comm_map),              pointer          :: t => null(), tp => null(), td => null()
    class(map_ptr),  allocatable, dimension(:)     :: theta_prev


    if (trim(self%type) == 'md') return

    call update_status(status, "mixupdate1 " // trim(self%label))
    
    ! Copy over alms from input structure, and compute pixel-space parameter maps
    if (present(theta)) then
       do i = 1, self%npar
          self%theta(i)%p%alm = theta(i)%alm          
       end do
    end if
    
    ! Compute mixing matrix
    do i = 1, numband
       
       ! Only update requested band if present
       if (present(band)) then
          if (i /= band) cycle
       end if

       ! Compute spectral parameters at the correct resolution for this channel
       if (self%npar > 0) then
          nmaps = min(data(i)%info%nmaps, self%theta(1)%p%info%nmaps)
          allocate(theta_p(0:data(i)%info%np-1,nmaps,self%npar))
          !if (self%myid==0) write(*,*) 'ee1 npix nmaps npar', data(i)%info%np-1,nmaps,self%npar
          
          do j = 1,self%npar
             info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                  & self%theta(j)%p%info%lmax, nmaps, data(i)%info%pol)
             td => comm_map(info)
             
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) >= 0)) then
                t => comm_map(info)
                t%alm(:,1:nmaps) = self%theta(j)%p%alm(:,1:nmaps)
                call t%Y_scalar()
                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) < 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do
                call t%dealloc(); deallocate(t)
             end if

             ! if any polarization is local sampled. Only set theta using polarizations with local sampling
             if (any(self%lmax_ind_pol(1:self%poltype(j),j) < 0)) then
                info_tp => comm_mapinfo(self%theta(j)%p%info%comm, self%theta(j)%p%info%nside, &
                     & 2*self%theta(j)%p%info%nside, nmaps, .false.)
                tp    => comm_map(info_tp)
                tp%map(:,1:nmaps) = self%theta(j)%p%map(:,1:nmaps)
                call tp%YtW_scalar()
                info => comm_mapinfo(data(i)%info%comm, data(i)%info%nside, &
                     & data(i)%info%lmax, nmaps, .false.)
                t    => comm_map(info)
                call tp%alm_equal(t)
                call t%Y_scalar()
                where (t%map < self%p_uni(1,j))
                   t%map = self%p_uni(1,j)
                end where
                where (t%map > self%p_uni(2,j))
                   t%map = self%p_uni(2,j)
                end where
                call tp%dealloc(); deallocate(tp)

                
                call wall_time(t2)
                do p = 1,self%poltype(j)
                   if (self%lmax_ind_pol(p,j) >= 0) cycle
                   if (self%poltype(j) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(j)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(j)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(j) 
                      stop
                   end if

                   do k = p_min,p_max
                      td%map(:,k) = t%map(:,k)
                   end do
                end do

                call t%dealloc(); deallocate(t)

                !if (info%myid == 0) write(*,*) 'udgrade = ', t2-t1
             end if
             theta_p(:,:,j) = td%map
             !if (info%myid == 0) write(*,*) 'q1, j=',j,  minval(theta_p(:,:,j)), maxval(theta_p(:,:,j))
             call td%dealloc(); deallocate(td)
          end do
       end if

       do l = 0, data(i)%ndet
          
          ! Don't update null mixing matrices
          if (self%F_null(i,l)) then
             if (present(df)) df(i)%p%map = 0.d0
             cycle
          end if
          
          
          ! Loop over all pixels, computing mixing matrix for each
          !allocate(theta_p(self%npar,self%nmaps))
          call wall_time(t1)
          do j = 0, self%F(i,l)%p%info%np-1
             if (self%latmask >= 0.d0) then
                call pix2ang_ring(data(i)%info%nside, data(i)%info%pix(j+1), lat, lon)
                lat = 0.5d0*pi-lat
                if (abs(lat) < self%latmask) then
                   self%F(i,l)%p%map(j,:) = 0.d0
                   cycle
                end if
             end if

!          if (all(self%lmax_ind_mix(1:min(self%nmaps,data(i)%info%nmaps)) == 0)) then  !if (self%lmax_ind == 0) then
!             cycle
!          end if

             ! NEW ! Check band sensitivity before mixing matrix update
             ! Possible labels are "broadband", "cmb", "synch", "dust", "co10", "co21", "co32", "ff", "ame"
             if (data(i)%comp_sens == "broadband") then
                ! If broadband, calculate mixing matrix
                mixmatnull = .false.
             else
                ! If component sensitivity, only calculate mixmat on that component.
                mixmatnull = .true.
                If (data(i)%comp_sens == self%label) then
                   mixmatnull = .false.
                end if
             end if
             
             ! Temperature
             if (.not. only_pol) then
                if (self%npar > 0) then
                   if (mixmatnull) then
                      self%F(i,l)%p%map(j,1) = 0.0
                   else
                      if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                         write(*,*) 'dust beta and T can not be null, crashing'
                         write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                         stop !debug, replace by proper stop and error message
                      end if
                      self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval(theta_p(j,1,:)) * data(i)%gain * self%cg_scale(1)
                   end if
                else
                   if (mixmatnull) then 
                      self%F(i,l)%p%map(j,1) = 0.0
                   else
                      self%F(i,l)%p%map(j,1) = self%F_int(1,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(1)
                   end if
                end if
             end if
             
             ! Polarization
             if (self%nmaps == 3 .and. data(i)%info%nmaps == 3) then
                ! Stokes Q
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else if (all(self%poltype < 2)) then
                   self%F(i,l)%p%map(j,2) = self%F(i,l)%p%map(j,1) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,2,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval(theta_p(j,2,:)) * data(i)%gain * self%cg_scale(2)
                   else
                      self%F(i,l)%p%map(j,2) = self%F_int(2,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(2)
                   end if
                end if
                
                ! Stokes U
                if (self%npar == 0) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else if (all(self%poltype < 3)) then
                   self%F(i,l)%p%map(j,3) = self%F(i,l)%p%map(j,2) 
                else
                   if (trim(self%label) == 'dust' .and. any(theta_p(j,1,:)==0.d0)) then
                      write(*,*) i, l, j, real(theta_p(j,1,:),sp)
                   end if
                   if (self%npar > 0) then
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval(theta_p(j,3,:)) * data(i)%gain * self%cg_scale(3)
                   else
                      self%F(i,l)%p%map(j,3) = self%F_int(3,i,l)%p%eval([0.d0]) * data(i)%gain * self%cg_scale(3)
                   end if
                end if
             end if
          
             ! Compute derivative if requested
             if (present(df) .and. l == 0) then
                if (self%npar > 0) then
                   do k = 1, nmaps
                      if (k <= self%poltype(par)) then
                         df(i)%p%map(j,k) = self%F_int(k,i,l)%p%eval_deriv(theta_p(j,k,:),par) * data(i)%gain * self%cg_scale(k)
                      end if
                   end do
                else
                   df(i)%p%map = 0.d0
                end if
             end if

          end do

          call wall_time(t2)
          !if (self%x%info%myid == 0) write(*,*) 'eval = ', t2-t1
                
          ! Compute mixing matrix average; for preconditioning only
          allocate(buffer(self%nmaps), buff2(self%nmaps))
          do j = 1, min(self%nmaps, data(i)%info%nmaps)
             self%F_mean(i,l,j) = sum(self%F(i,l)%p%map(:,j))
          end do
          buff2 = self%F_mean(i,l,:)
          !call mpi_barrier(mpi_comm_world, ierr)
          call mpi_allreduce(buff2, buffer, self%nmaps, &
               & MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
          self%F_mean(i,l,:) = buffer / self%F(i,l)%p%info%npix
          deallocate(buffer,buff2)
    
       end do
       if (allocated(theta_p)) deallocate(theta_p)
    end do

    call update_status(status, "mixupdate2 " // trim(self%label))

    ! Request preconditioner update
    recompute_diffuse_precond = .true.

  end subroutine updateDiffuseMixmat



  module function evalDiffuseBand(self, band, amp_in, pix, alm_out, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: res

    integer(i4b) :: i, j, np, nmaps, lmax, nmaps_comp, d
    logical(lgt) :: alm_out_
    real(dp)     :: t1, t2
    class(comm_mapinfo), pointer :: info 
    class(comm_map),     pointer :: m

    alm_out_ = .false.; if (present(alm_out)) alm_out_ = alm_out
    d = 0; if (present(det)) d = det

    if (self%F_null(band,0)) then
       if (alm_out_) then
          if (.not. allocated(res)) allocate(res(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       else
          if (.not. allocated(res)) allocate(res(0:data(band)%info%np-1,data(band)%info%nmaps))
       end if
       res = 0.d0
       return
    end if

    ! Initialize amplitude map
    nmaps =  min(data(band)%info%nmaps, self%nmaps)
    info  => comm_mapinfo(data(band)%info%comm, data(band)%info%nside, data(band)%info%lmax, nmaps, nmaps==3)
    m     => comm_map(info)
    if (present(amp_in)) then
       m%alm(:,1:nmaps) = amp_in(:,1:nmaps)
       !m%alm(:,1:nmaps) = amp_in
    else
       call self%x%alm_equal(m)
       !m%alm(:,1:nmaps) = self%x%alm(:,1:nmaps)
    end if

    if (apply_mixmat) then
       ! Scale to correct frequency through multiplication with mixing matrix

       if (all(self%lmax_ind_mix(1:nmaps,:) == 0) .and. self%latmask < 0.d0) then
          do i = 1, m%info%nmaps
             m%alm(:,i) = m%alm(:,i) * self%F_mean(band,d,i)
          end do
       else
          call m%Y()
          m%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band,d)%p%map(:,1:nmaps)
          call m%YtW()
       end if
    end if

    ! Convolve with band-specific beam
    call data(band)%B(d)%p%conv(trans=.false., map=m)
       
    ! Return correct data product
    if (alm_out_) then
       if (.not. data(band)%B(d)%p%almFromConv) call m%YtW()
       if (.not. allocated(res)) allocate(res(0:data(band)%info%nalm-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) res = 0.d0
       res(:,1:nmaps) = m%alm(:,1:nmaps)
    else
       if (data(band)%B(d)%p%almFromConv) call m%Y()
       if (.not. allocated(res)) allocate(res(0:data(band)%info%np-1,data(band)%info%nmaps))
       if (nmaps /= data(band)%info%nmaps) res = 0.d0
       res(:,1:nmaps) = m%map(:,1:nmaps)
    end if

    ! Clean up
    call m%dealloc(); deallocate(m)
    nullify(info)

  end function evalDiffuseBand

  ! Return component projected from map
  module function projectDiffuseBand(self, band, map, alm_in, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: res

    integer(i4b) :: i, nmaps, d
    logical(lgt) :: alm_in_
    class(comm_mapinfo), pointer :: info_in => null(), info_out => null()
    class(comm_map),     pointer :: m       => null(), m_out    => null()

    if (self%F_null(band,0)) then
       if (.not. allocated(res)) allocate(res(0:self%x%info%nalm-1,self%x%info%nmaps))
       res = 0.d0
       return
    end if

    d = 0; if (present(det)) d = det

    nmaps     =  min(self%x%info%nmaps, map%info%nmaps)
    info_in   => comm_mapinfo(self%x%info%comm, map%info%nside, map%info%lmax, nmaps, nmaps==3)
    info_out  => comm_mapinfo(self%x%info%comm, map%info%nside, self%lmax_amp, nmaps, nmaps==3)
    m         => comm_map(info_in)
    m_out     => comm_map(info_out)
    alm_in_ = .false.; if (present(alm_in)) alm_in_ = alm_in

    ! Convolve with band-specific beam
    if (alm_in) then
       m%alm = map%alm(:,1:nmaps)
    else
       m%map = map%map(:,1:nmaps)
       call m%Yt()
    end if
    call data(band)%B(d)%p%conv(trans=.true., map=m)
    
    if (all(self%lmax_ind_mix(1:nmaps,:) == 0) .and. self%latmask < 0.d0) then
       do i = 1, nmaps
          m%alm(:,i) = m%alm(:,i) * self%F_mean(band,d,i)
       end do
    else
       if (data(band)%B(d)%p%almFromConv) call m%Y()
       m%map(:,1:nmaps) = m%map(:,1:nmaps) * self%F(band,d)%p%map(:,1:nmaps)
       call m%YtW()
    end if
    call m%alm_equal(m_out)

    if (.not. allocated(res)) allocate(res(0:self%x%info%nalm-1,self%x%info%nmaps))
    res = m_out%alm

    call m%dealloc(); deallocate(m)
    call m_out%dealloc(); deallocate(m_out)

  end function projectDiffuseBand


  module subroutine applyDiffPrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    select case (trim(precond_type))
    case ("diagonal")
       call applyDiffPrecond_diagonal(x)
    case ("pseudoinv")
       call applyDiffPrecond_pseudoinv(x)
    case default
       call report_error("Preconditioner type not supported: "//trim(precond_type))
    end select

  end subroutine applyDiffPrecond


  module subroutine applyDiffPrecond_diagonal(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y

    if (npre == 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
    allocate(y(npre,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       call cr_extract_comp(diffComps(i)%p%id, x, alm)
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Multiply with preconditioner
    do j = 1, nmaps_pre
       do i = 0, info_pre%nalm-1
          if (P_cr%invM_diff(i,j)%n == 0) cycle
          y(P_cr%invM_diff(i,j)%ind,i,j) = &
               & matmul(P_cr%invM_diff(i,j)%M, y(P_cr%invM_diff(i,j)%ind,i,j))
       end do
    end do

    ! Reformat y(npre,nalm,nmaps) structure into linear array
    do i = 1, npre
       nmaps = diffComps(i)%p%x%info%nmaps
       allocate(alm(0:diffComps(i)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(i)%p%x%info%nalm-1
          call diffComps(i)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = y(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(i)%p%id, .false., alm, x)
       deallocate(alm)
    end do
    
    deallocate(y)

  end subroutine applyDiffPrecond_diagonal


  module subroutine applyDiffPrecond_pseudoinv(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, ii, j, k, l, m, p, q, qq, nmaps, npre_int
    real(dp)                  :: t1, t2
    real(dp), allocatable, dimension(:)     :: w, w2
    real(dp), allocatable, dimension(:,:)   :: alm
    real(dp), allocatable, dimension(:,:,:) :: y, z
    class(comm_map), pointer                :: invN_x => null()

    if (.not. allocated(ind_pre)) return

    npre_int = size(ind_pre)
    if (npre_int <= 0) return
    
    ! Reformat linear array into y(npre,nalm,nmaps) structure
!    call update_status(status, "pseudo1")
    allocate(y(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    allocate(z(npre_int,0:info_pre%nalm-1,info_pre%nmaps))
    y = 0.d0
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       call cr_extract_comp(diffComps(ii)%p%id, x, alm)
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          y(i,k,1:nmaps) = alm(j,1:nmaps)
       end do
       deallocate(alm)
    end do

    ! Frequency-dependent terms
    z = 0.d0
    do k = 1, numband
!       call update_status(status, "pseudo2")
       invN_x => comm_map(data(k)%info)
       nmaps  =  data(k)%info%nmaps
       
       ! Sum over (U^plus)^t
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                invN_x%alm(i,p) = invN_x%alm(i,p) + P_cr%invM_diff(l,p)%M(qq,k) * y(q,j,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
!       call update_status(status, "pseudo3")

       ! Multiply by T
       call invN_x%WY
!       call update_status(status, "pseudo4")
       call data(k)%N%N(invN_x)
!       call update_status(status, "pseudo5")
       call invN_x%YtW
!       call update_status(status, "pseudo6")
       do i = 1, nmaps
          invN_x%alm(:,i) = invN_x%alm(:,i) * data(k)%N%alpha_nu(i)**2
       end do

       ! Sum over U^plus
       !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,m,j,q,qq,p)
       !!$OMP DO SCHEDULE(guided)
       do i = 0, data(k)%info%nalm-1
          call data(k)%info%i2lm(i, l, m)
          if (l > info_pre%lmax) cycle
          call info_pre%lm2i(l,m,j)
          do q = 1, npre_int
             qq = ind_pre(q)
             do p = 1, nmaps
                z(q,j,p) = z(q,j,p) + P_cr%invM_diff(l,p)%M(qq,k) * invN_x%alm(i,p)
             end do
          end do
       end do
       !!$OMP END DO
       !!$OMP END PARALLEL
!       call update_status(status, "pseudo7")

       call invN_x%dealloc(); deallocate(invN_x)
    end do

    ! Prior terms
!    call update_status(status, "pseudo7.1")
    call wall_time(t1)
    !!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,l,k,j,m,w,w2,p)
    allocate(w(npre_int), w2(npre_int))
    !!$OMP DO SCHEDULE(guided)
    do i = 0, info_pre%nalm-1
       do p = 1, info_pre%nmaps
          !call info_pre%i2lm(i, l, m)
          l = info_pre%lm(1,i)
          m = info_pre%lm(2,i)
          w        = y(:,i,p)
          w2       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w2(j) = w2(j) + P_cr%invM_diff(l,p)%M(ind_pre(k),numband+ind_pre(j))*w(k)
             end do
          end do
          w       = 0.d0
          do j = 1, npre_int
             do k = 1, npre_int
                w(j) = w(j) + P_cr%invM_diff(l,p)%M(ind_pre(j),numband+ind_pre(k))*w2(k)
             end do
          end do
          z(:,i,p) = z(:,i,p) + w
       end do
    end do
    deallocate(w, w2)
    !!$OMP END DO
    !!$OMP END PARALLEL
    call wall_time(t2)
    !if (info_pre%myid == 0 .or. info_pre%myid == 25) write(*,*) info_pre%myid, ', nalm = ', info_pre%nalm, real(t2-t1,sp)
!    call update_status(status, "pseudo8")

    ! Reformat z(npre,nalm,nmaps) structure into linear array
    do i = 1, npre_int
       ii = ind_pre(i)
       nmaps = diffComps(ii)%p%x%info%nmaps
       allocate(alm(0:diffComps(ii)%p%x%info%nalm-1,nmaps))
       alm = 0.d0
       do j = 0, diffComps(ii)%p%x%info%nalm-1
          call diffComps(ii)%p%x%info%i2lm(j, l, m)
          call info_pre%lm2i(l, m, k)
          alm(j,1:nmaps) = z(i,k,1:nmaps)
       end do
       call cr_insert_comp(diffComps(ii)%p%id, .false., alm, x)
       deallocate(alm)
    end do
!    call update_status(status, "pseudo9")
    
    deallocate(y, z)

  end subroutine applyDiffPrecond_pseudoinv


  
  ! Dump current sample to HEALPix FITS file
  module subroutine dumpDiffuseToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    !
    ! Routine that writes a diffuce component to FITS (and HDF) files. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! iter: integer
    !       Sample number in the Gibb's chain.
    !
    ! chainfile: hdf_file
    !       HDF file to write the component to
    !
    ! output_hdf: logical
    !       Logical parameter to tell whether or not to write the component to the specified HDF file
    !
    ! postfix: string
    !       A string label to be added to the end of FITS-files.
    !       (default format: cXXXX_kYYYYYY; XXXX = chain number, YYYYYY = sample number)
    !
    ! dir: string
    !       Output directory to which output is written
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, j, k, m, ierr, unit
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: vals(10)
    logical(lgt)       :: exist, first_call = .true.
    character(len=6)   :: itext
    character(len=512) :: filename, path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer :: map => null(), tp => null()
    real(dp), allocatable, dimension(:,:) :: sigma_l
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg

    if (.not. self%output) return

    if (trim(self%type) == 'md') then
       filename = trim(dir)//'/md_' // trim(postfix) // '.dat'
       vals = 0.d0
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l == 0) then                 ! Monopole
             vals(1)  = 1.d0/sqrt(4.d0*pi) * self%x%alm(i,1)             * self%RJ2unit_(1)
             vals(5)  = 1.d0/sqrt(4.d0*pi) * self%mu%alm(i,1)            * self%RJ2unit_(1)
             vals(9)  = 1.d0/sqrt(4.d0*pi) * sqrt(self%Cl%Dl(0,1))       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m == -1) then   ! Y dipole
             vals(3)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(7)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m ==  0) then   ! Z dipole
             vals(4)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(8)  = 1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
          end if
          if (l == 1 .and. m ==  1) then   ! X dipole
             vals(2)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%x%alm(i,1)        * self%RJ2unit_(1)
             vals(6)  = -1.d0/sqrt(4.d0*pi/3.d0) * self%mu%alm(i,1)       * self%RJ2unit_(1)
             vals(10) =  1.d0/sqrt(4.d0*pi/3.d0) * sqrt(self%Cl%Dl(1,1))  * self%RJ2unit_(1)
          end if
       end do
       call mpi_allreduce(MPI_IN_PLACE, vals, 10, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       if (self%x%info%myid == 0) then
          if (output_hdf) then
             call int2string(iter, itext)
             path = trim(adjustl(itext))//'/md/'//trim(adjustl(self%label))
             call write_hdf(chainfile, trim(adjustl(path)), vals(1:4))
          end if
          
          inquire(file=trim(filename), exist=exist)
          unit = getlun()
          if (first_call) then
             open(unit, file=trim(filename), recl=1024)
             write(unit,*) '# Band          M_def       X_def       Y_def'//&
                  & '      Z_def     M_mean      X_mean      Y_mean      Z_mean'// &
                  & '      M_rms       D_rms'
             first_call = .false.
          else
             open(unit, file=trim(filename), recl=1024, position='append')
          end if
          write(unit,'(a12,10e12.3)') trim(self%label), vals
          close(unit)
       end if
    else

       call update_status(status, "writeFITS_1")

       ! Write amplitude
       map => comm_map(self%x)
       do i = 1, map%info%nmaps
          map%alm(:,i) = map%alm(:,i) * self%RJ2unit_(i) * self%cg_scale(i)  ! Output in requested units
       end do

       !call update_status(status, "writeFITS_2")

       if (output_hdf) then
          call int2string(iter, itext)
          path = '/'//trim(adjustl(itext))//'/'//trim(adjustl(self%label))
          !write(*,*) 'path', trim(path)
          if (self%x%info%myid == 0) call create_hdf_group(chainfile, trim(adjustl(path)))
       end if

       !call update_status(status, "writeFITS_3")

       filename = trim(self%label) // '_' // trim(postfix) // '.fits'
       call self%B_out%conv(trans=.false., map=map)
       call map%Y
       do i = 1, map%info%nmaps
          map%alm(:,i) = self%x%alm(:,i) * self%RJ2unit_(i) * self%cg_scale(i)  ! Replace convolved with original alms
       end do
       !call update_status(status, "writeFITS_4")

       !call self%apply_proc_mask(map)

       if (output_hdf) then
          !write(*,*) 'path2', trim(path)//'/amp_'
          call map%writeFITS(trim(dir)//'/'//trim(filename), &
               & hdffile=chainfile, hdfpath=trim(path)//'/amp_', output_hdf_map=.false.)
          !if we have set the overall scale parameter
          if (self%x_scale /= 1.d0 .and. self%myid == 0) then
            call write_hdf(chainfile, trim(path)//'/x_scale', self%x_scale)
          end if
       else
          call map%writeFITS(trim(dir)//'/'//trim(filename))
       end if
       call map%dealloc(); deallocate(map)
       call update_status(status, "writeFITS_5")

       if (self%output_EB) then
          map => comm_map(self%x)

          do i = 1, map%info%nmaps
             map%alm(:,i) = map%alm(:,i) * self%RJ2unit_(i) * self%cg_scale(i)  ! Output in requested units
          end do
          
          filename = trim(self%label) // '_' // trim(postfix) // '_TEB.fits'
          call self%B_out%conv(trans=.false., map=map)
          call map%Y_EB
          !call self%apply_proc_mask(map)
          call map%writeFITS(trim(dir)//'/'//trim(filename))
          call map%dealloc(); deallocate(map)
       end if
       !call update_status(status, "writeFITS_6")
       
       allocate(sigma_l(0:self%x%info%lmax,self%x%info%nspec))
       call self%x%getSigmaL(sigma_l)
       k = 1
       do i = 1, self%x%info%nmaps
          do j = i, self%x%info%nmaps
             sigma_l(:,k) = sigma_l(:,k) * self%RJ2unit(j)
             k = k+1
          end do
       end do
       if (self%x%info%myid == 0) then
          filename = trim(dir)//'/sigma_l_'//trim(self%label) // '_' // trim(postfix) // '.dat'
          call write_sigma_l(filename, sigma_l)
          if (output_hdf) call write_hdf(chainfile, trim(adjustl(path))//'/sigma_l', sigma_l)             
       end if
       deallocate(sigma_l)
       !call update_status(status, "writeFITS_7")

       ! Write spectral index maps
       do i = 1, self%npar
          filename = trim(self%label) // '_' // trim(self%indlabel(i)) // '_' // &
               & trim(postfix) // '.fits'

          ! need to only take alms from polarizations with lmax > 0
          if (self%lmax_ind >= 0) then
             info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                  & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
             
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (all(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                call self%theta(i)%p%Y_scalar()
             else if (any(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                tp => comm_map(info)
                tp%alm = self%theta(i)%p%alm
                call tp%Y_scalar
                do p = 1,self%poltype(i)
                   if (self%lmax_ind_pol(p,i) < 0) cycle
                   if (self%poltype(i) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(i)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(i)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                      stop
                   end if

                   do k = p_min,p_max
                      self%theta(i)%p%map(:,k) = tp%map(:,k)
                   end do
                end do
                call tp%dealloc(); deallocate(tp)
             end if
          end if


          if (output_hdf) then
             call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename), &
                  & hdffile=chainfile, hdfpath=trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
                  & '_')
          else
             call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
          end if
          
          !write proposal length and number of proposals maps if local sampling was used
          if (self%output_localsamp_maps .and. any(self%lmax_ind_pol(:min(self%nmaps,self%poltype(i)),i) < 0 .and. &
               & self%pol_pixreg_type(:min(self%nmaps,self%poltype(i)),i) > 0)) then
             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_proplen_'  // trim(postfix) // '.fits'
             call self%pol_proplen(i)%p%writeFITS(trim(dir)//'/'//trim(filename))

             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_nprop_'  // trim(postfix) // '.fits'
             call self%pol_nprop(i)%p%writeFITS(trim(dir)//'/'//trim(filename))

          end if

          !if pixelregions, create map without smoothed thetas (for input in new runs)
          if (self%output_localsamp_maps .and. any(self%pol_pixreg_type(1:min(self%nmaps,self%poltype(i)),i) > 0)) then
             
             info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                  & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
             tp => comm_map(info)
             tp%map = self%theta(i)%p%map
             do p = 1,self%poltype(i)
                if (self%pol_pixreg_type(p,i) /=3) cycle
                if (self%poltype(i) == 1) then
                   p_min=1
                   p_max=info%nmaps
                   if (only_pol) p_min = 2
                else if (self%poltype(i)==2) then
                   if (p == 1) then
                      p_min = 1
                      p_max = 1
                   else
                      p_min = 2
                      p_max = info%nmaps
                   end if
                else if (self%poltype(i)==3) then
                   p_min = p
                   p_max = p
                else
                   write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                   stop
                end if

                do j = 0,info%np-1
                   tp%map(j,p_min:p_max) = self%theta_pixreg(self%ind_pixreg_arr(j,p,i),p,i)
                end do
             end do
             filename = trim(self%label) // '_' // trim(self%indlabel(i)) // &
                  & '_noSmooth_'  // trim(postfix) // '.fits'
             call tp%writeFITS(trim(dir)//'/'//trim(filename))
             call tp%dealloc(); deallocate(tp)

          end if

          !output theta, proposal length and number of proposals per pixel region to HDF
          if (output_hdf) then
             npol=min(self%nmaps,self%poltype(i))!only concerned about the maps/poltypes in use
             if (any(self%pol_pixreg_type(:npol,i) > 0)) then
                npr=0
                do j = 1,npol
                   if (self%npixreg(j,i)>npr) npr = self%npixreg(j,i)
                end do
                if (npr == 0) then !no pixelregions, theta = prior
                   if (self%theta(i)%p%info%myid == 0) write(*,*) 'No defined pixel regions for ',trim(self%label)//'_'//&
                        & trim(self%indlabel(i))//' in writing to HDF'
                else
                   allocate(dp_pixreg(npr,npol),int_pixreg(npr,npol))
                   !pixel region values for theta
                   dp_pixreg=self%theta_pixreg(1:npr,1:npol,i)
                   !if (self%theta(i)%p%info%myid == 0) write(*,*) 'zzz1', npr, npol, i
                   !if (self%theta(i)%p%info%myid == 0) write(*,*) 'zzz2', minval(dp_pixreg(:,1)), maxval(dp_pixreg(:,1))
                   !if (self%theta(i)%p%info%myid == 0) write(*,*) 'zzz3', minval(dp_pixreg(:,2)), maxval(dp_pixreg(:,2))
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_val', real(dp_pixreg,sp))
                   !pixel region values for proposal length
                   dp_pixreg=self%proplen_pixreg(1:npr,1:npol,i)
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_proplen', real(dp_pixreg,sp))
                   !pixel region values for number of proposals
                   int_pixreg=self%nprop_pixreg(1:npr,1:npol,i)
                   if (self%theta(i)%p%info%myid == 0) call write_hdf(chainfile, trim(path)//'/'//&
                        & trim(adjustl(self%indlabel(i)))//'_pixreg_nprop', int_pixreg)

                   deallocate(dp_pixreg,int_pixreg)
                end if
             end if

          end if




       end do

       ! Output Sampled SED's
       if (output_hdf .and. allocated(self%SEDtab) .and. self%x%info%myid == 0) then
         call write_hdf(chainfile, trim(path)//'/SED', self%SEDtab)
       end if
       
       ! Write mixing matrices
       if (self%output_mixmat) then
          do i = 1, numband
             if (self%F_null(i,0)) cycle
             filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
                  & trim(postfix) // '.fits'
             call self%F(i,0)%p%writeFITS(trim(dir)//'/'//trim(filename))
          end do
       end if
       call update_status(status, "writeFITS_9")
    end if

  end subroutine dumpDiffuseToFITS

  ! Dump current sample to HEALPix FITS file
  module subroutine initDiffuseHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_diffuse_comp),  intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath

    integer(i4b)       :: i, j, l, m, ierr
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: md(4)
    character(len=512) :: path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer     :: tp => null()
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg
    
    if (trim(self%type) == 'md') then
       call read_hdf(hdffile, trim(adjustl(hdfpath))//'md/'//trim(adjustl(self%label)), md)
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l == 0) then                 
             self%x%alm(i,1) =  md(1) * sqrt(4.d0*pi)      / self%RJ2unit_(1)  ! Monopole
          else if (l == 1 .and. m == -1) then   
             self%x%alm(i,1) =  md(3) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! Y dipole
          else if (l == 1 .and. m ==  0) then   
             self%x%alm(i,1) =  md(4) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! Z dipole
          else if (l == 1 .and. m ==  1) then   
             self%x%alm(i,1) = -md(2) * sqrt(4.d0*pi/3.d0) / self%RJ2unit_(1)  ! X dipole
          end if
       end do
    else
       path = trim(adjustl(hdfpath))//trim(adjustl(self%label))
       call self%Cl%initHDF(hdffile, path)
       call self%x%readHDF(hdffile, trim(adjustl(path))//'/amp_alm', .false.)
       !call self%x%readHDF(hdffile, trim(adjustl(path))//'/amp_map', .true.)    ! Read amplitudes
       do i = 1, self%x%info%nmaps
         self%x%alm(:,i) = self%x%alm(:,i) / (self%RJ2unit_(i) * self%cg_scale(i))
       end do
       do i = 0, self%x%info%nalm-1
          call self%x%info%i2lm(i,l,m)
          if (l < self%lmin_amp) self%x%alm(i,:) = 0.d0
       end do

       if (trim(self%type) == 'MBBtab') then
         call read_hdf(hdffile, trim(adjustl(path))//'/SED', self%SEDtab)
       end if

       do i = 1, self%npar
          call self%theta(i)%p%readHDF(hdffile, trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
               & '_map', .true.)
          if (self%lmax_ind >= 0) then
             call self%theta(i)%p%readHDF(hdffile, trim(path)//'/'//trim(adjustl(self%indlabel(i)))//&
                  & '_alm', .false.)
             ! need to only take alms from polarizations with lmax > 0
             ! if any polarization is alm sampled. Only use alms to set polarizations with alm sampling
             if (all(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                call self%theta(i)%p%Y_scalar
             else if (any(self%lmax_ind_pol(1:self%poltype(i),i) >= 0)) then
                info => comm_mapinfo(self%theta(i)%p%info%comm, self%theta(i)%p%info%nside, &
                     & self%theta(i)%p%info%lmax, self%theta(i)%p%info%nmaps, self%theta(i)%p%info%pol)
                tp => comm_map(info)
                tp%alm = self%theta(i)%p%alm
                call tp%Y_scalar
                do p = 1,self%poltype(i)
                   if (self%lmax_ind_pol(p,i) < 0) cycle
                   if (self%poltype(i) == 1) then
                      p_min=1
                      p_max=info%nmaps
                      if (only_pol) p_min = 2
                   else if (self%poltype(i)==2) then
                      if (p == 1) then
                         p_min = 1
                         p_max = 1
                      else
                         p_min = 2
                         p_max = info%nmaps
                      end if
                   else if (self%poltype(i)==3) then
                      p_min = p
                      p_max = p
                   else
                      write(*,*) '  Unknown poltype in component ',self%label,', parameter ',self%indlabel(i) 
                      stop
                   end if

                   do j = p_min,p_max
                      self%theta(i)%p%map(:,j) = tp%map(:,j)
                   end do
                end do
                call tp%dealloc(); deallocate(tp)
             end if
          end if

          !Need to initialize pixelregions and local sampler from chain as well (where relevant)
          npol=min(self%nmaps,self%poltype(i))!only concerned about the maps/poltypes in use
          if (any(self%pol_pixreg_type(:npol,i) > 0) .and. cpar%sample_specind) then
             npr=0
             do j = 1,npol
                if (self%npixreg(j,i)>npr) npr = self%npixreg(j,i)
             end do
             if (npr == 0) then !no pixelregions, theta = prior
                if (self%theta(i)%p%info%myid == 0) write(*,*) 'No defined pixel regions for ',trim(self%label)//'_'//&
                     & trim(self%indlabel(i))
             else
                allocate(dp_pixreg(npr,npol),int_pixreg(npr,npol))
                !pixel region values for theta
                if (self%theta(i)%p%info%myid == 0) call read_hdf_dp_2d_buffer(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_val', dp_pixreg)
                call mpi_bcast(dp_pixreg, size(dp_pixreg),  MPI_DOUBLE_PRECISION, 0, self%theta(i)%p%info%comm, ierr)
                self%theta_pixreg(1:npr,1:npol,i)=dp_pixreg
                !pixel region values for proposal length
                if (self%theta(i)%p%info%myid == 0) call read_hdf_dp_2d_buffer(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_proplen', dp_pixreg)
                call mpi_bcast(dp_pixreg, size(dp_pixreg),  MPI_DOUBLE_PRECISION, 0, self%theta(i)%p%info%comm, ierr)
                self%proplen_pixreg(1:npr,1:npol,i)=dp_pixreg
                !pixel region values for number of proposals
                if (self%theta(i)%p%info%myid == 0) call read_hdf_int_2d_buffer(hdffile, trim(path)//'/'//&
                     & trim(adjustl(self%indlabel(i)))//'_pixreg_nprop', int_pixreg)
                call mpi_bcast(int_pixreg, size(int_pixreg),  MPI_INTEGER, 0, self%theta(i)%p%info%comm, ierr)
                self%nprop_pixreg(1:npr,1:npol,i)=int_pixreg

                deallocate(dp_pixreg,int_pixreg)
             end if
          end if
       end do !i = 1,npar
    end if

    call self%updateMixmat

  end subroutine initDiffuseHDF
  
  module subroutine add_to_npre(n, nside, lmax, nmaps)
    implicit none
    integer(i4b), intent(in) :: n, nside, lmax, nmaps
    npre      = npre + n
    nside_pre = min(nside_pre, nside)
    lmax_pre  = max(lmax_pre, lmax)
    nmaps_pre = max(nmaps_pre, nmaps)
  end subroutine add_to_npre

  ! Sample spectral parameters. This subroutine is obsolete, it is defined and used in comm_nonlin_mod instead
  module subroutine sampleDiffuseSpecInd(self, cpar, handle, id, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter

    integer(i4b) :: i, j, k, l, n, p, q, pix, ierr, ind(1), counter, n_ok
    integer(i4b) :: i_min, i_max, status, n_gibbs, n_pix, n_pix_tot, flag, npar, np, nmaps
    real(dp)     :: a, b, a_tot, b_tot, s, t1, t2, x_min, x_max, delta_lnL_threshold
    real(dp)     :: mu, sigma, par, w, mu_p, sigma_p, a_old, chisq, chisq_old, chisq_tot, unitconv
    real(dp)     :: x(1), theta_min, theta_max
    logical(lgt) :: ok
    logical(lgt), save :: first_call = .true.
    class(comm_comp), pointer :: c => null()
    real(dp),     allocatable, dimension(:)   :: lnL, P_tot, F, theta, a_curr
    real(dp),     allocatable, dimension(:,:) :: amp, buffer, alm_old
    !Following is for the local sampler
    real(dp)     :: old_theta, new_theta, mixing_old, mixing_new, lnL_new, lnL_old, res_lnL, delta_lnL, accept_rate
    integer(i4b) :: i_s, p_min, p_max, pixreg_nprop
    integer(i4b) :: n_spec_prop, n_accept, n_corr_prop, n_prop_limit, n_corr_limit, corr_len
    logical(lgt) :: first_sample
    class(comm_mapinfo),            pointer :: info => null()

    return 

  end subroutine sampleDiffuseSpecInd
  
  module subroutine print_precond_mat
    implicit none

    integer(i4b) :: l, m, i, j, p
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: mat

    if (info_pre%myid /= 0) return
    p = 2

    open(58,file=trim(outdir)//'/precond_W.dat',recl=1024)
    if (trim(precond_type) == 'diagonal') then
       do l = 0, info_pre%lmax
          call info_pre%lm2i(l, 0, i)
          write(*,*) 
          write(*,*) l 
          do j = 1, size(P_cr%invM_diff(i,p)%M(j,:),1)
             write(*,*) real(P_cr%invM_diff(i,p)%M(j,:),sp)
          end do
          allocate(W(P_cr%invM_diff(i,p)%n))
          call get_eigenvalues(P_cr%invM_diff(i,p)%M, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
    else
       allocate(mat(npre,npre))
       do l = 0, info_pre%lmax
          mat = matmul(P_cr%invM_diff(l,p)%M, transpose(P_cr%invM_diff(l,p)%M))
          write(*,*) 
          write(*,*) l 
          do j = 1, npre
             write(*,*) real(mat(j,:),sp)
          end do
          allocate(W(npre))
          call get_eigenvalues(mat, W)
          write(58,*) l, real(W,sp)
          deallocate(W)
       end do
       deallocate(mat)
    end if
    close(58)


  end subroutine print_precond_mat

  module subroutine updateDiffuseFInt(self, band)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: band

    integer(i4b) :: i, j, k

    if (present(band)) then
       do i = 1, data(band)%info%nmaps
          do j = 0, data(band)%ndet
             call self%F_int(i,band,j)%p%update(pol=i)
          end do
       end do
    else
       do k = 1, numband
          do i = 1, data(k)%info%nmaps
             do j = 0, data(k)%ndet
                call self%F_int(i,k,j)%p%update(pol=i)
             end do
          end do
       end do
    end if

  end subroutine updateDiffuseFInt

  module subroutine updateLowlPrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, lp, mp, myid, nalm, ierr, nmaps
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer

    if (self%lmax_pre_lowl < 0) return

    call wall_time(t1)

    ! Initialize arrays
    nalm = (self%lmax_pre_lowl+1)**2
    allocate(invM(0:nalm-1,0:nalm-1), buffer(0:nalm-1,0:nalm-1))
    invM = 0.d0

    info => comm_mapinfo(self%comm, 2, 2*self%lmax_pre_lowl, self%nmaps, self%nmaps==3)
    map  => comm_map(info)
    tot  => comm_map(info)
    do l = 0, self%lmax_pre_lowl
       !if (info%myid == 0) write(*,*) '  Low-ell init l = ', l, self%lmax_pre_lowl, trim(self%label)
       do m = -l, l
          ! Set up unit vector
          call info%lm2i(l,m,j)
          map%alm = 0.d0
          if (j >= 0) map%alm(j,1) = 1.d0

          call self%Cl%sqrtS(alm=map%alm, info=map%info) ! Multiply with sqrt(Cl)

          ! Add frequency dependent terms
          tot%alm = 0.d0
          do i = 1, numband

             ! Compute component-summed map, ie., column-wise matrix elements
             nmaps    = min(self%nmaps, data(i)%info%nmaps)
             info     => comm_mapinfo(self%comm, data(i)%N%nside_chisq_lowres, 2*self%lmax_pre_lowl, &
                  & nmaps, nmaps==3)
             map2     => comm_map(info)
             do k = 1, nmaps
                map2%alm(:,k) = map%alm(:,k) * self%F_mean(i,0,k)
             end do
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'a', l, m, i
             end if
             call data(i)%B(0)%p%conv(trans=.false., map=map2)
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'b', l, m, i
             end if
             call map2%Y()
             if (any(map2%map /= map2%map)) then
                write(*,*) 'c', l, m, i
             end if

             ! Multiply with invN
             call data(i)%N%InvN_lowres(map2)
             if (any(map2%map /= map2%map)) then
                write(*,*) 'd', l, m, i
             end if

             ! Project summed map into components, ie., row-wise matrix elements
             call map2%Yt()
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'e', l, m, i
             end if
             call data(i)%B(0)%p%conv(trans=.true., map=map2)
             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'f', l, m, i
             end if

             map2%alm(:,1) = map2%alm(:,1) * self%F_mean(i,0,1)

             if (any(map2%alm /= map2%alm)) then
                write(*,*) 'g/', l, m, i
             end if


             ! Add up alms
             tot%alm(:,1) = tot%alm(:,1) + map2%alm(:,1)

             call map2%dealloc(); deallocate(map2)
          end do

             if (any(tot%alm /= tot%alm)) then
                write(*,*) 'h', l, m
             end if

          ! Add prior term and multiply with sqrt(S) for relevant components
          call self%Cl%sqrtS(alm=tot%alm, info=tot%info)

             if (any(tot%alm /= tot%alm)) then
                write(*,*) 'i', l, m
             end if


          ! Add (unity) prior term, and store column
          i = l**2 + l + m
          do lp = 0, self%lmax_pre_lowl
             do mp = -lp, lp
                call tot%info%lm2i(lp,mp,k)
                if (k >= 0) then
                   j = lp**2 + lp + mp
                   invM(i,j) = tot%alm(k,1)
                   if (i == j) invM(i,j) = invM(i,j) + 1.d0
                end if
             end do
          end do

       end do
    end do

    call mpi_allreduce(invM, buffer, size(invM), MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
    invM = buffer

    if (self%x%info%myid == 0) then
       do i = 0, 3
          write(*,*) real(invM(i,0:3),sp)
       end do
    end if

    call invert_matrix(invM, cholesky=.true.)

    if (self%x%info%myid == 0) then
       do i = 0, 3
!          write(*,*) real(invM(i,0:3),sp)
       end do
    end if

    ! Store matrix rows
    q = 0
    do l = 0, self%lmax_pre_lowl
       do m = -l, l
          call self%x%info%lm2i(l,m,j)
          if (j >= 0) then
             k = l**2 + l + m
             do i = 0, nalm-1
                self%invM_lowl(i,q) = invM(i,k)
             end do
             q = q+1
          end if
       end do
    end do

    deallocate(invM, buffer)
    call map%dealloc(); deallocate(map)
    call tot%dealloc(); deallocate(tot)

    call wall_time(t2)
    if (info%myid == 0) write(*,*) '  Low-ell init = ', t2-t1


  end subroutine updateLowlPrecond

  module subroutine applyLowlPrecond(self, alm, alm0)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(inout)  :: alm
    real(dp),                 dimension(0:,1:), intent(in)     :: alm0

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, q, myid, nalm, ntot, ierr
    real(dp), allocatable, dimension(:) :: y, yloc, buffer

    ntot = (self%lmax_pre_lowl+1)**2
    allocate(y(0:ntot-1))
    if (allocated(self%invM_lowl)) then
       nalm = size(self%invM_lowl,2)
       allocate(yloc(0:nalm-1))
       q   = 0
       do l = 0, self%lmax_pre_lowl
          do m = -l, l
             call self%x%info%lm2i(l,m,j)
             if (j >= 0) then
                yloc(q) = alm(j,1)
                q       = q+1
             end if
          end do
       end do
       y    = matmul(self%invM_lowl,yloc)
       deallocate(yloc)
    else
       nalm = 0
       y    = 0.d0
    end if
    allocate(buffer(0:ntot-1))
    call mpi_allreduce(y, buffer, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, &
         & self%x%info%comm, ierr)
    y = buffer
    deallocate(buffer)

    alm = alm0
    do l = 0, self%lmax_pre_lowl
       do m = -l, l
          call self%x%info%lm2i(l,m,j)
          if (j >= 0) then
             k = l**2 + l + m
             alm(j,1) = y(k)
          end if
       end do
    end do


    deallocate(y)

  end subroutine applyLowlPrecond


  module subroutine updateDeflatePrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self

    real(dp)                  :: t1, t2, norm
    integer(i4b)              :: i, j, k, l, m, p, q, lp, mp, myid, nalm, ierr, nmaps, nmax, ndef, nside
    logical(lgt)              :: add
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()
    real(dp),        allocatable, dimension(:,:) :: invM, buffer, Z
    real(dp),        allocatable, dimension(:) :: W

    if (self%lmax_def < 0) return

    call wall_time(t1)

    if (.not. allocated(self%Z_def)) then
       nside     = min(self%nside,128)
       self%ndef = 0
       q         = (nside/self%nside_def)**2
       info => comm_mapinfo(self%comm, nside, self%lmax_def, 1, .false.)
       map  => comm_map(info)

       nalm      = info%nalm
       nmax      = min(10000, info%npix)
       allocate(Z(0:nalm-1,nmax))

       call setup_needlets(info, self%nside_def, self%defmask, Z, self%ndef) 
       if (self%myid == 0) write(*,*) 'nz = ', self%ndef

       allocate(self%Z_def(0:nalm-1,self%ndef))
       self%Z_def = Z(:,1:self%ndef)

       if (self%myid == 0) allocate(self%invM_def(self%ndef,self%ndef))

       call map%dealloc(); deallocate(map)
       deallocate(Z)
    end if

    ndef = size(self%Z_def,2)
    info => comm_mapinfo(self%comm, 2, self%lmax_def, self%nmaps, self%nmaps==3)
    map  => comm_map(info)
    tot  => comm_map(info)
    allocate(invM(ndef,ndef), buffer(ndef,ndef))
    do l = 1, ndef

       if (info%myid == 0 .and. mod(l,100) == 1) write(*,*) '  Precomputing deflate mode = ', l, ndef, trim(self%label)
       ! Set up basis vector
       map%alm = 0.d0
       map%alm(:,1) = self%Z_def(:,l)

       call self%Cl%sqrtS(alm=map%alm, info=map%info) ! Multiply with sqrt(Cl)

       ! Add frequency dependent terms
       tot%alm = 0.d0
       do i = 1, numband

          ! Compute component-summed map, ie., column-wise matrix elements
          nmaps  = data(i)%info%nmaps
          info     => comm_mapinfo(self%comm, data(i)%N%nside_chisq_lowres, self%lmax_def, &
               & nmaps, nmaps==3)
          map2     => comm_map(info)
          do k = 1, min(self%nmaps, nmaps)
             map2%alm(:,k) = map%alm(:,k) * self%F_mean(i,0,k)
          end do
          call data(i)%B(0)%p%conv(trans=.false., map=map2)
          call map2%Y()

          ! Multiply with invN
          call data(i)%N%InvN_lowres(map2)

             ! Project summed map into components, ie., row-wise matrix elements
          call map2%Yt()
          call data(i)%B(0)%p%conv(trans=.true., map=map2)
          map2%alm(:,1) = map2%alm(:,1) * self%F_mean(i,0,1)

          ! Add up alms
          do k = 1, min(self%nmaps, nmaps)
             tot%alm(:,k) = tot%alm(:,k) + map2%alm(:,k)
          end do

          call map2%dealloc(); deallocate(map2)
       end do

       ! Add prior term and multiply with sqrt(S) for relevant components
       call self%Cl%sqrtS(alm=tot%alm, info=tot%info)

       ! Add (unity) prior term, and store column
       tot%alm(:,1) = tot%alm(:,1) + self%Z_def(:,l)

       if (info%nalm > 0) then
          do j = 1, ndef
             invM(l,j) = sum(self%Z_def(:,j)*tot%alm(:,1))
          end do
!          write(*,*) self%myid, l, sum(invM(l,:)), sum(tot%alm(:,1))
       else
          invM(l,:) = 0.d0
       end if
          
    end do

    call mpi_reduce(invM, buffer, size(invM), MPI_DOUBLE_PRECISION, MPI_SUM, 0, info%comm, ierr)

    if (self%myid == 0) then
       self%invM_def = buffer
       allocate(W(ndef))
       call get_eigenvalues(buffer, W)
       write(*,*) 'W =', real(W,sp)
       deallocate(W)
       call invert_matrix(self%invM_def, cholesky=.true.)
    end if

    deallocate(invM, buffer)
    call map%dealloc(); deallocate(map)
    call tot%dealloc(); deallocate(tot)

    call wall_time(t2)
    if (info%myid == 0) write(*,*) '  Deflate init = ', t2-t1

  end subroutine updateDeflatePrecond



  module subroutine applyDeflatePrecond(self, alm, Qalm)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(in)     :: alm
    real(dp),                 dimension(0:,1:), intent(out)    :: Qalm

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, k, l, m, myid, nalm, ntot, ierr, ndef
    real(dp), allocatable, dimension(:) :: y, ytot
    class(comm_map),     pointer :: map => null(), map2 => null(), tot => null()
    class(comm_mapinfo), pointer :: info => null()

    Qalm= 0.d0
    if (self%lmax_def < 0) return

    ndef = size(self%Z_def,2)
    allocate(y(ndef), ytot(ndef))

    ! Multiply with Z
    info => comm_mapinfo(self%comm, 2, self%lmax_def, 1, .false.)
    map  => comm_map(info)
    map%alm = 0.d0
    do i = 0, info%nalm-1
       call info%i2lm(i,l,m)
       call self%x%info%lm2i(l,m,j)
       if (j > -1) map%alm(i,1) = alm(j,1)
    end do
    if (info%nalm > 0) then
       y = matmul(transpose(self%Z_def), map%alm(:,1))
    else
       y = 0.d0
    end if
    call mpi_reduce(y, ytot, ndef, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
         & self%x%info%comm, ierr)

    ! Multiply with invE
    if (self%myid == 0) then
       y = matmul(self%invM_def, ytot)
    end if
    call mpi_bcast(y, ndef, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)

    ! Multiply with Z^t
    map%alm(:,1) = matmul(self%Z_def, y)
    
    do i = 0, info%nalm-1
       call info%i2lm(i,l,m)
       call self%x%info%lm2i(l,m,j)
       if (j > -1) Qalm(j,1) = map%alm(i,1)
    end do
    
    call map%dealloc(); deallocate(map)
    deallocate(y, ytot)

  end subroutine applyDeflatePrecond

  module subroutine setup_needlets(info, nside_def, defmask, Z, ndef)
    implicit none
    class(comm_mapinfo), pointer,          intent(in)  :: info
    integer(i4b),                          intent(in)  :: nside_def
    class(comm_map),                       intent(in)  :: defmask
    real(dp),            dimension(0:,1:), intent(out) :: Z
    integer(i4b),                          intent(out) :: ndef

    integer(i4b) :: i, j, l, p, q, nside, m, res, ierr
    real(dp)     :: B, bl, Bj
    character(len=4) :: ntext
    logical(lgt) :: add
    class(comm_map),     pointer :: map => null()
    real(dp), allocatable, dimension(:) :: x, t, f, psi, phi, b0
    type(spline_type) :: sb, spsi, sphi

    ! Set up needlet kernel
    m    = 10000
    B    = 2.72d0
    allocate(x(m), t(m), f(m), psi(m), phi(m), b0(m))
    do i = 1, m
       x(i) = -1.d0 + 2.d0*(i-1.d0)/(m-1.d0)
    end do
    f = exp(-1.d0/(1.d0-x**2))

    psi(1) = 0.d0
    do i = 2, m
       psi(i) = psi(i-1) + 0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
    end do
    psi = psi / psi(m)
    call spline(spsi, x, psi)

    do i = 1, m
       t(i) = (i-1.d0)/(m-1.d0)
       if (t(i) < 1.d0/B) then
          phi(i) = 1.d0
       else
          phi(i) = splint(spsi, 1.d0-2.d0*B/(B-1.d0)*(t(i)-1.d0/B))
       end if
    end do
    call spline(sphi, t, phi)

    do i = 1, m
       t(i)  = B*(i-1.d0)/(m-1.d0)
       if (t(i) < 1.d0) then
          b0(i) = sqrt(splint(sphi, t(i)/B) - splint(sphi, t(i)))
          !b0(i) = (splint(sphi, t(i)/B) - splint(sphi, t(i)))
       else
          b0(i) = sqrt(splint(sphi, t(i)/B))
          !b0(i) = (splint(sphi, t(i)/B))
       end if
       if (b0(i) < 1d-10 .or. b0(i) /= b0(i)) b0(i) = 0.d0
    end do
    call spline(sb, t, b0)

    map   => comm_map(info)
    ndef  = 0.d0
    nside = 1
    Bj    = B
    do while(nside <= nside_def)
       do i = 0, 12*nside**2-1
          p = i*(defmask%info%nside/nside)**2
          call nest2ring(defmask%info%nside, p, j)
          add = .false.
          do p = 1, defmask%info%np
             if (defmask%info%pix(p) == j) add = (defmask%map(p-1,1) == 0.d0)  
          end do
          call mpi_allreduce(MPI_IN_PLACE, add, 1, MPI_LOGICAL, MPI_LOR, info%comm, ierr)
          if (.not. add) cycle
          
          ndef    = ndef+1
          if (ndef > size(Z,2)) then
             write(*,*) ' ERROR: Too many deflated basis functions'
             stop
          end if

          map%map = 0.d0
          p       = i*(info%nside/nside)**2
          call nest2ring(info%nside, p, j)
          do p = 1, info%np
             if (info%pix(p) == j) map%map(p,1) = 1.d0
          end do
          call map%YtW

          do j = 0, info%nalm-1
             l = info%lm(1,j)
             if (l > Bj) then
                bl = 0.d0
             else
                bl = splint(sb, l/Bj)
             end if
             map%alm(j,1) = bl * map%alm(j,1) 
          end do

          if (mod(ndef,10) == 0) then
             call int2string(ndef, ntext)
             call map%Y
             call map%writeFITS('need'//ntext//'.fits')
          end if

          Z(:,ndef) = map%alm(:,1)
       end do
       if (info%myid == 0) write(*,*) nside, ndef

       nside = 2*nside
       Bj    = B*Bj
    end do


    call free_spline(sb)
    call free_spline(spsi)
    call free_spline(sphi)
    call map%dealloc(); deallocate(map)
    deallocate(x, t, f, psi, phi, b0)

  end subroutine setup_needlets

  module subroutine applyMonoDipolePrior(self, handle)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    type(planck_rng),         intent(inout)          :: handle

    integer(i4b) :: i, j, k, l, m, ierr, pix
    real(dp)     :: mu(0:3), a, b, Amat(0:3,0:3), bmat(0:3), v(0:3), corr_res(3)
    class(comm_map), pointer :: map, lr_map 
    class(comm_mapinfo), pointer :: info => null()
    real(dp), dimension(:), allocatable :: mask_list, corr_list, amp_list, intersect, all_thetas
    real(dp)     :: mean_intersect, std_intersect 
    real(dp)     :: diff_mono, diff_comp, mono_mix, comp_mix, prior_vals(2)
    class(comm_comp), pointer :: c => null()

    if (trim(self%mono_prior_type) == 'none') then ! No active monopole prior
       return
    end if

    ! Compute monopole offset given the specified prior
    mu = 0.d0

    ! Generate real-space component map
    map => comm_map(self%x)
    if (trim(self%mono_prior_type) == 'crosscorr' .or. trim(self%mono_prior_type) == 'lower_value_prior') then
       call self%B_mono_prior%conv(trans=.false., map=map) !smooth to prior specific FWHM
    else
       call self%B_out%conv(trans=.false., map=map) !smooth to output FWHM of component
    end if
    call map%Y

    if (trim(self%mono_prior_type) == 'monopole') then        ! Set monopole to zero outside user-specified mask

       ! Compute mean outside mask; no noise weighting for now at least
       a = sum(map%map(:,1) * self%mono_prior_map%map(:,1))
       b = sum(self%mono_prior_map%map(:,1))
       call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       mu(0) = a / b
       mu(1:3) = 0.d0 !in order to not subtract any dipole in alm space!

       ! Subtract mean in real space
       self%x%map(:,1) = self%x%map(:,1) - mu(0)

       if (self%x%info%myid == 0) then
             write(*,fmt='(a)') ' |  Monopole prior correction for component: '//trim(self%label)
             write(*,fmt='(a,f11.3)') ' |   Monopole: ',mu(0)*self%cg_scale(1)
             write(*,fmt='(a)') ' | '
       end if

    else if (trim(self%mono_prior_type) == 'monopole-dipole' .or. trim(self%mono_prior_type) == 'monopole+dipole') then        ! Set monopole or monopole+dipole to zero outside user-specified mask. In both cases the dipole is computed, but only in the monopole+dipole case it is removed
       ! Generate real-space component map

       ! Compute mean outside mask; no noise weighting for now at least
       Amat = 0.d0; bmat = 0.d0
       do i = 0, self%x%info%np-1
          if (self%mono_prior_map%map(i,1) < 0.5d0) cycle
          v(0) = 1.d0
          call pix2vec_ring(self%x%info%nside, self%x%info%pix(i+1), v(1:3))
          do j = 0, 3
             do k = 0, 3
                Amat(j,k) = Amat(j,k) + v(j)*v(k)
             end do
             bmat(j) = bmat(j) + v(j)*map%map(i,1)
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, Amat, 16, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, bmat,  4, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call solve_system_real(Amat, mu, bmat)

       ! Subtract mean (and dipole) in real space
       if (trim(self%mono_prior_type) == 'monopole-dipole') then
          ! Subtract mean in real space 
          self%x%map(:,1) = self%x%map(:,1) - mu(0)
          if (self%x%info%myid == 0) then
             write(*,fmt='(a)') ' |  Monopole prior correction (with dipole estimate) for component: '//trim(self%label)
             write(*,fmt='(a,f10.3,a,3f10.3,a)') ' |    Monopole (dipole):', &
                  & mu(0)*self%cg_scale(1),'  ( ',mu(1:3)*self%cg_scale(1), ' )'
             write(*,fmt='(a)') ' | '

          end if
          mu(1:3)=0.d0 !in order to not subtract the dipole in alm space!
       else
          ! Subtract mean and dipole in real space
          do i = 0, self%x%info%np-1
             v(0) = 1.d0
             call pix2vec_ring(self%x%info%nside, self%x%info%pix(i+1), v(1:3))
             self%x%map(i,1) = self%x%map(i,1) - sum(v*mu)
          end do
          if (self%x%info%myid == 0) then
             write(*,fmt='(a)') ' Monopole+dipole prior correction for component: '//trim(self%label)
             write(*,fmt='(a)') '      Monopole   Dipole_x   Dipole_y   Dipole_z'
             write(*,fmt='(a,4f11.3)') '   ',mu*self%cg_scale(1)
             write(*,fmt='(a)') ' | '
          end if
       end if

    else if (trim(self%mono_prior_type) == 'monopole+tempfit') then 

       ! Compute mean outside mask; no noise weighting for now at least
       Amat = 0.d0; bmat = 0.d0
       do i = 0, self%x%info%np-1
          if (self%mono_prior_map%map(i,2) < 0.5d0) cycle
          Amat(1,1) = Amat(1,1) + 1.d0
          Amat(2,1) = Amat(2,1) + self%mono_prior_map%map(i,1)
          Amat(1,2) = Amat(1,2) + self%mono_prior_map%map(i,1)
          Amat(2,2) = Amat(2,2) + self%mono_prior_map%map(i,1)**2 
          bmat(1)   = bmat(1)   + map%map(i,1)
          bmat(2)   = bmat(2)   + map%map(i,1)*self%mono_prior_map%map(i,1)
       end do
       call mpi_allreduce(MPI_IN_PLACE, Amat(1:2,1:2), 4, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, bmat(1:2),     2, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call solve_system_real(Amat(1:2,1:2), mu(1:2), bmat(1:2))

    else if (trim(self%mono_prior_type) == 'crosscorr') then ! Enforce zero intercept in correlation with specified map
       !ud-grade amplitude map if necessary
       info     => comm_mapinfo(self%comm, self%mono_prior_nside, -1, self%nmaps, self%pol)
       lr_map => comm_map(info)

       call map%udgrade(lr_map)

       allocate(mask_list(0:lr_map%info%npix-1))
       allocate(corr_list(0:lr_map%info%npix-1))
       allocate(amp_list(0:lr_map%info%npix-1))

       !1. for each threshold:
       allocate(intersect(self%mono_prior_Nthresh))
       m = 0
       do i = 1,self%mono_prior_Nthresh
          !resetting the lists
          mask_list = 0.d0
          corr_list = 0.d0
          amp_list = 0.d0

          !creating mask
          do j = 0,lr_map%info%np-1
             !    1: mask pixels above threshold
             if (self%mono_prior_map%map(j,1) > self%mono_prior_threshold(i)) then
                mask_list(self%mono_prior_map%info%pix(j+1)) = 0.d0
             else
                mask_list(self%mono_prior_map%info%pix(j+1)) = 1.d0
             end if
          end do
          !    2: set up a pixel list for the corr-map and comp. amplitude (like ud-grade) in order to gather all values
          corr_list(self%mono_prior_map%info%pix) = self%mono_prior_map%map(:,1)
          amp_list(lr_map%info%pix) = lr_map%map(:,1)
          call mpi_allreduce(MPI_IN_PLACE, mask_list, lr_map%info%npix, MPI_DOUBLE_PRECISION, MPI_SUM, lr_map%info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, corr_list, lr_map%info%npix, MPI_DOUBLE_PRECISION, MPI_SUM, lr_map%info%comm, ierr)
          call mpi_allreduce(MPI_IN_PLACE, amp_list, lr_map%info%npix, MPI_DOUBLE_PRECISION, MPI_SUM, lr_map%info%comm, ierr)

          !    4: root processor computes the intersection (and slope and correlation coefficient) and stores value to a list
          !restructure so that we only have unmasked pixels, only 1 proc needs to do this
          if (lr_map%info%myid == 0) then
             k=-1

             !open(58, file='corr.dat') !debug
             ! set up lists with only the pixels with cross-correlation data under the threshold
             do j = 0,lr_map%info%npix-1
                if (mask_list(j) > 0.5d0) then
                   k = k + 1
                   corr_list(k) = corr_list(j) !external data
                   amp_list(k) = amp_list(j)   !component amplitude
                   !write(58,*) corr_list(j), amp_list(j)  !debug
                end if
             end do
             !close(58) !debug

             !calculate intersection 
             corr_res = calc_linear_regression(corr_list(0:k), amp_list(0:k))
             !write(*,*) 'threshold',self%mono_prior_threshold(i),'corr_res', corr_res

             if (corr_res(3) <= 1.0d0 .and. corr_res(3) >= -1.0d0) then
                m = m + 1
                intersect(m) = corr_res(1)
             end if
          end if
       end do !i = 1,self%mono_prior_Nthresh

       ! calculate mean and std of intersect values from regression 
       if (lr_map%info%myid == 0) then
          if (m <= 0) call report_error("No intersection value found from crosscorrelation, monopole prior component "//trim(self%label))
          if (m == 1) then
             mu(0) = intersect(1)
             mean_intersect = intersect(1)
             std_intersect = 0.d0
          else
             mean_intersect = sum(intersect(1:m))/m
             std_intersect = sqrt(sum((intersect(1:m)-mean_intersect)**2)/m)

             if (std_intersect > 0.d0) then
                a = rand_gauss(handle)     
                mu(0) = mean_intersect + std_intersect * a
             else
                mu(0) = mean_intersect
             end if
          end if
       end if


       ! bcast drawn monopole value/intersection value
       call mpi_bcast(mu(0), 1, MPI_DOUBLE_PRECISION, 0, lr_map%info%comm, ierr)
       call mpi_bcast(mean_intersect, 1, MPI_DOUBLE_PRECISION, 0, lr_map%info%comm, ierr)
       call mpi_bcast(std_intersect, 1, MPI_DOUBLE_PRECISION, 0, lr_map%info%comm, ierr)
       
       !mu(0) = mu(0) - 50
       mu(1:3) = 0.d0

       ! Subtract mean in real space 
       self%x%map(:,1) = self%x%map(:,1) - mu(0)
       if (self%x%info%myid == 0) then
          write(*,fmt='(a)') ' |  Cross-correlation prior correction for component: '//trim(self%label)
          write(*,fmt='(a,i2)') ' |    Number of linear fits (thresholds): ',&
               & self%mono_prior_Nthresh
          write(*,fmt='(a,f14.3,f14.3)') ' |    Drawing intersect to subtract from prior (mu,RMS)  ', &
               & mean_intersect*self%cg_scale(1), &
               & std_intersect*self%cg_scale(1) 
          write(*,fmt='(a,f14.3,f14.3)') ' |    New value             ', &
               & (mean_intersect-mu(0))*self%cg_scale(1)
          write(*,fmt='(a,f14.3,f14.3)') ' |    Old value             ', &
               & mean_intersect*self%cg_scale(1)
          write(*,fmt='(a,f14.3,f14.3)') ' |    Difference            ', &
               & -mu(0)*self%cg_scale(1)
          write(*,fmt='(a)') ' |  '

       end if

       deallocate(intersect)
       deallocate(mask_list)
       deallocate(amp_list)
       deallocate(corr_list)
       call lr_map%dealloc(); deallocate(lr_map)

    else if (trim(self%mono_prior_type) == 'lower_value_prior') then
       
       !allocate map lists for the mask and smoothed amplitude map
       allocate(mask_list(0:map%info%npix-1))
       allocate(amp_list(0:map%info%npix-1))
       mask_list = 0.d0
       amp_list = 0.d0

       do j = 0,map%info%np-1
          !    1: mask pixels above threshold
          if (self%mono_prior_map%map(j,1) > 0.5d0) then
             mask_list(self%mono_prior_map%info%pix(j+1)) = 1.d0
          else
             mask_list(self%mono_prior_map%info%pix(j+1)) = 0.d0
          end if
       end do
       amp_list(map%info%pix) = map%map(:,1)
       call mpi_allreduce(MPI_IN_PLACE, mask_list, map%info%npix, MPI_DOUBLE_PRECISION, MPI_SUM, map%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, amp_list, map%info%npix, MPI_DOUBLE_PRECISION, MPI_SUM, map%info%comm, ierr)

       !    4: root processor calculates the lowest value inside the mask, draws a new value from the input prior and finds the monopole needed to make the lowest value equal to the drawn value

       if (map%info%myid == 0) then
          k=-1
          do j = 0,map%info%npix-1
             if (mask_list(j) > 0.5d0) then
                if (k < 0) then
                   k = j
                else
                   if (amp_list(j) < amp_list(k)) k = j
                end if
             end if
          end do
          !amp_list(k) is now the lowest value in the smoothed map (inside the mask)

          if (k < 0) call report_error("No lower value found from lowest value prior on amplitude monopole prior correction. All pixels masked. Component "//trim(self%label))
          !Draw new lower value
          if (self%mono_prior_gaussian_rms > 0.d0 ) then
             mean_intersect = -1.d0
             do while (mean_intersect < 0.d0)
                a = rand_gauss(handle)     
                mean_intersect = self%mono_prior_gaussian_mean + &
                     & self%mono_prior_gaussian_rms * a
             end do
          else
             mean_intersect = self%mono_prior_gaussian_mean
          end if
          mu(0) = amp_list(k) - mean_intersect 
          !when subrtracting this mu(0), the lowest value in the map becomes mean_intersect 
       end if


       ! bcast drawn monopole value/intersection value
       call mpi_bcast(mu(0), 1, MPI_DOUBLE_PRECISION, 0, map%info%comm, ierr)
       call mpi_bcast(mean_intersect, 1, MPI_DOUBLE_PRECISION, 0, map%info%comm, ierr)
       
       mu(1:3) = 0.d0

       ! Subtract mean in real space 
       self%x%map(:,1) = self%x%map(:,1) - mu(0)
       if (self%x%info%myid == 0) then
          write(*,fmt='(a)') ' Lowest value prior correction for component: '//trim(self%label)
          write(*,fmt='(a,f14.3,f14.3)') '   Prior value (mu,RMS)  ', &
               & self%mono_prior_gaussian_mean*self%cg_scale(1), &
               & self%mono_prior_gaussian_rms*self%cg_scale(1) 
          write(*,fmt='(a,f14.3,f14.3)') '   New value             ', &
               & mean_intersect*self%cg_scale(1)
          write(*,fmt='(a,f14.3,f14.3)') '   Old value             ', amp_list(k)*self%cg_scale(1)
          write(*,fmt='(a,f14.3,f14.3)') '   Difference            ', -mu(0)*self%cg_scale(1)
          write(*,fmt='(a)') ' | '
       end if

       deallocate(mask_list)
       deallocate(amp_list)

    else if (trim(self%mono_prior_type) == 'bandmono') then
       !get monopole value of frequency-band monopole and draw a new value
       ! real(dp) :: mu(0:3), a, b, Amat(0:3,0:3), bmat(0:3), v(0:3), corr_res(3)
       c => compList
       a=0.d0
       b=0.d0
       prior_vals=0.d0
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             if (trim(c%label)==trim(self%mono_prior_band)) then
                mono_mix=c%RJ2unit_(1)
                do i = 0, c%x%info%nalm-1
                   call c%x%info%i2lm(i,l,m)
                   if (l == 0) then ! monopole
                      a = c%x%alm(i,1)
                      ! draw new monopole value 
                      !mono [alm_uKRJ] = mu_in_alm_uK_RJ + rms_in_alm_uKRJ * rand_gauss
                      b = c%mu%alm(i,1) + sqrt(c%Cl%Dl(0,1))*rand_gauss(handle) 
                      c%x%alm(i,1) = b ! set new monopole value
                      prior_vals(1)=c%mu%alm(i,1)/sqrt(4.d0*pi)
                      prior_vals(2)=sqrt(c%Cl%Dl(0,1))/sqrt(4.d0*pi)
                   end if
                end do
             end if
          end select
          c => c%nextComp()
       end do

       ! MPI reduce existing and new monopole
       call mpi_allreduce(MPI_IN_PLACE, a, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
       call mpi_allreduce(MPI_IN_PLACE, prior_vals, 2, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)



       diff_mono = (b - a)/sqrt(4.d0*pi) !to get it in pixel space units

       !need to get mixing matrix for data band for the given component
       allocate(all_thetas(self%npar))
       comp_mix=0.d0
       do j = 1,numband
          if (trim(self%mono_prior_band)==trim(data(j)%label)) then
             do i = 1, self%npar
                mean_intersect = 0.d0
                do pix = 0,self%x%info%np-1
                   mean_intersect = mean_intersect + self%theta(i)%p%map(pix,1) 
                end do
                call mpi_allreduce(MPI_IN_PLACE, mean_intersect, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
                all_thetas(i) = mean_intersect/self%x%info%npix
             end do

             ! get conversion factor from amplitude to data 
             ! (i.e. mixing matrix element for prior band)         
             comp_mix = self%F_int(1,j,0)%p%eval(all_thetas) * &
                        & data(j)%gain * self%cg_scale(1)
             exit
          end if
       end do
       
       if (comp_mix > 0.d0) then
          diff_comp = diff_mono * mono_mix/comp_mix 
          !diff[monoRJ]*[monoRJ2data]/[compRJ2data] = diff[compRJ]
       else
          diff_comp = 0.d0
          if (self%x%info%myid == 0) write(*,fmt='(a)') &
               &'   Could not compute mixing matrix to prior frequency band "'&
               &//trim(self%mono_prior_band)//'" for component "'//trim(self%label)//&
               & '". Reversing the sampled monopole.'

          c => compList
          do while (associated(c))
             select type (c)
             class is (comm_diffuse_comp)
                if (trim(c%label)==trim(self%mono_prior_band)) then
                   mono_mix=c%RJ2unit_(1)
                   do i = 0, c%x%info%nalm-1
                      call c%x%info%i2lm(i,l,m)
                      if (l == 0) then ! monopole
                         c%x%alm(i,1) = a ! go back to original value
                      end if
                   end do
                end if
             end select
             c => c%nextComp()
          end do
       end if

       mu = 0.d0
       mu(0) = diff_comp
       ! a positive shift in band monopole -> negative shift in amplitude 
       ! -> subtract the monopole

       ! Subtract mean in real space 
       self%x%map(:,1) = self%x%map(:,1) - mu(0)
       if (self%x%info%myid == 0) then 
          write(*,fmt='(a)') ' |  Band monopole prior correction for -- comp: '//trim(self%label)//' -- prior band: '//trim(self%mono_prior_band)
          write(*,fmt='(a,f14.3,f14.3)') ' |  Band monopole prior (mu,RMS) ', prior_vals(1),prior_vals(2)
          write(*,fmt='(a,f14.3)') ' |  New band monopole            ',b*mono_mix/sqrt(4.d0*pi)
          write(*,fmt='(a,f14.3)') ' |  Old band monopole            ',a*mono_mix/sqrt(4.d0*pi)
          write(*,fmt='(a,f14.3)') ' |  Change to band monopole      ',diff_mono*mono_mix
          write(*,fmt='(a,f14.3)') ' |  Change to component monopole ',-diff_comp*self%cg_scale(1)
          write(*,fmt='(a)') ' | '
       end if

       deallocate(all_thetas)

    end if
    
    ! Prepare template corrected map in harmonic space
    if (trim(self%mono_prior_type) == 'monopole+tempfit') then
       map%map      = 0.d0
       map%map(:,1) = mu(1) + mu(2)*self%mono_prior_map%map(:,1)
       call map%YtW()
       self%x%alm(:,1) = self%x%alm(:,1) - map%alm(:,1)
    else
       do i = 0, self%x%info%nalm-1
          if (self%x%info%lm(1,i) == 0 .and. self%x%info%lm(2,i) == 0) then
             self%x%alm(i,1) = self%x%alm(i,1) - mu(0) * sqrt(4.d0*pi)
          end if
          if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == -1) then
             self%x%alm(i,1) = self%x%alm(i,1) - mu(2) * sqrt(4.d0*pi/3.d0)
          end if
          if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 0) then
             self%x%alm(i,1) = self%x%alm(i,1) - mu(3) * sqrt(4.d0*pi/3.d0)
          end if
          if (self%x%info%lm(1,i) == 1 .and. self%x%info%lm(2,i) == 1) then
             self%x%alm(i,1) = self%x%alm(i,1) + mu(1) * sqrt(4.d0*pi/3.d0)
          end if
       end do
    end if

    call map%dealloc(); deallocate(map)
    
  end subroutine applyMonoDipolePrior

  module subroutine nullify_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band

    integer(i4b) :: i
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             do i = 0, c%x%info%nalm-1
                if (c%x%info%lm(1,i) == 0 .and. c%x%info%lm(2,i) == 0) then
                   c%x%alm(i,1) = 0.d0
                   return
                end if
             end do
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine nullify_monopole_amp


  module function get_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band
    real(dp)                     :: get_monopole_amp

    integer(i4b) :: l, m, ierr
    real(dp)     :: mono
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             mono = -1.d100
             if (c%x%info%nalm > 0) then
                call c%x%info%i2lm(0,l,m)
                if (l == 0) then
                   mono = 1.d0/sqrt(4.d0*pi) * c%x%alm(0,1) * c%RJ2unit_(1)
                end if
             end if
             call mpi_allreduce(MPI_IN_PLACE, mono, 1, MPI_DOUBLE_PRECISION, MPI_MAX, c%x%info%comm, ierr)
             get_monopole_amp = mono
             return
          end if
       end select
       c => c%nextComp()
    end do

  end function get_monopole_amp

  module subroutine set_monopole_amp(band, mono)
    implicit none
    character(len=*), intent(in) :: band
    real(dp),         intent(in) :: mono

    integer(i4b) :: l, m
    class(comm_comp), pointer :: c => null()

    c => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (trim(c%label) == trim(band)) then
             if (c%x%info%nalm > 0) then
                call c%x%info%i2lm(0,l,m)
                if (l == 0) then
                   c%x%alm(0,1) = sqrt(4.d0*pi) * mono /c%RJ2unit_(1)
                end if
             end if
             return
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine set_monopole_amp

end submodule comm_diffuse_comp_smod
