module comm_spindust_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  use spline_1D_mod
  implicit none

  private
  public comm_spindust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_spindust_comp
     real(dp)          :: nu_p0, nu_min, nu_max
     type(spline_type) :: SED_spline
   contains
     procedure :: S    => evalSED
  end type comm_spindust_comp

  interface comm_spindust_comp
     procedure constructor
  end interface comm_spindust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_spindust_comp), pointer   :: constructor

    integer(i4b) :: i, j, k, l, m, n, ierr, ind(1)
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta 
    character(len=512) :: temptxt, partxt
    real(dp), allocatable, dimension(:,:) :: SED
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: theta_single_pol => null() 
    class(comm_map),     pointer :: theta_single_pol_smooth => null() 

    ! General parameters
    allocate(constructor)
    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    constructor%npar         = 1
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%poltype(1), constructor%indlabel(1))
    allocate(constructor%nu_min_ind(1), constructor%nu_max_ind(1))
    constructor%poltype(1)   = cpar%cs_poltype(1,id_abs)
    constructor%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    constructor%indlabel(1)  = 'nu_p'
    constructor%nu_min_ind(1) = cpar%cs_nu_min(id_abs,1)
    constructor%nu_max_ind(1) = cpar%cs_nu_max(id_abs,1)

    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    ! Component specific parameters for 2 parameter model
    !constructor%npar         = 2
    !allocate(constructor%theta_def(2), constructor%p_gauss(2,2), constructor%p_uni(2,2))
    !allocate(constructor%poltype(2), constructor%indlabel(2))
    !allocate(constructor%nu_min_ind(2), constructor%nu_max_ind(2))
    !do i = 1, 2
    !   constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    !   constructor%theta_def(i) = cpar%cs_theta_def(i,id_abs)
    !   constructor%p_uni(:,i)   = cpar%cs_p_uni(id_abs,:,i)
    !   constructor%p_gauss(:,i) = cpar%cs_p_gauss(id_abs,:,i)
    !   constructor%nu_min_ind(i) = cpar%cs_nu_min(id_abs,i)
    !   constructor%nu_max_ind(i) = cpar%cs_nu_max(id_abs,i)
    !end do
    !constructor%indlabel  = ['nu_p','alpha']
    
    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(1))
    
    if (trim(cpar%cs_input_ind(1,id_abs)) == 'default') then
       constructor%theta(1)%p => comm_map(info)
       constructor%theta(1)%p%map = constructor%theta_def(1)
    else
       ! Read map from FITS file, and convert to alms
       constructor%theta(1)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(1,id_abs)))
    end if
    if (constructor%lmax_ind >= 0) call constructor%theta(1)%p%YtW_scalar

    ! Initialize spectral index map for 2 parameter model
    !allocate(constructor%theta(constructor%npar))
    !do i = 1, constructor%npar
    !   if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
    !      constructor%theta(i)%p => comm_map(info)
    !      constructor%theta(i)%p%map = constructor%theta_def(i)
    !   else
    !      ! Read map from FITS file, and convert to alms
    !      constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_$
    !   end if
    !   if (constructor%lmax_ind >= 0) call constructor%theta(1)%p%YtW_scalar
    !end do

    ! Initialize spectral template !CHANGE??
    call read_spectrum(trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(1,id_abs)), SED)
    ind                = maxloc(SED(:,2))
    constructor%nu_p0  = SED(ind(1),1)
    constructor%nu_min = minval(SED(:,1))
    constructor%nu_max = maxval(SED(:,1))
    SED                = log(SED)
    call spline(constructor%SED_spline, SED(:,1), SED(:,2))
    deallocate(SED)

    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                   constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             constructor%F_int(k,i,j)%p => comm_F_int_1D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    ! Initialize mixing matrix
    call constructor%updateMixmat

    ! Set up smoothing scale information
    allocate(constructor%smooth_scale(constructor%npar))
    constructor%smooth_scale = cpar%cs_smooth_scale(id_abs,1:constructor%npar)

    !!! (local) sampling specific parameters!!!

    allocate(constructor%pol_lnLtype(3,constructor%npar))        ! {chisq, ridge, marginal}: evaluation type for lnL
    allocate(constructor%pol_sample_nprop(3,constructor%npar))   ! {.true., .false.}: sample nprop on first iteration
    allocate(constructor%pol_sample_proplen(3,constructor%npar)) ! {.true., .false.}: sample proplen on first iteration
    allocate(constructor%pol_pixreg_type(3,constructor%npar))    ! {1=fullsky, 2=single_pix, 3=pixel_regions}
    allocate(constructor%nprop_uni(2,constructor%npar))          ! {integer}: upper and lower limits on nprop
    allocate(constructor%lmax_ind_pol(3,constructor%npar))       ! {integer}: lmax per. poltype sample per spec. index
    allocate(constructor%npixreg(3,constructor%npar))            ! {integer}: number of pixel regions per poltye per spec ind

    do i = 1,constructor%npar
       constructor%nprop_uni(:,i)=cpar%cs_spec_uni_nprop(:,i,id_abs)
       do j = 1,constructor%poltype(i)
          !assign lmax per spec ind per polarization sample type (poltype)
          constructor%lmax_ind_pol(j,i) = cpar%cs_lmax_ind_pol(j,i,id_abs)
          !find the highest lmax of the component in total
          constructor%lmax_ind = max(constructor%lmax_ind,constructor%lmax_ind_pol(j,i)) 
          constructor%pol_lnLtype(j,i)  = cpar%cs_spec_lnLtype(j,i,id_abs)
          constructor%pol_sample_nprop(j,i) = cpar%cs_spec_samp_nprop(j,i,id_abs)
          constructor%pol_sample_proplen(j,i) = cpar%cs_spec_samp_proplen(j,i,id_abs)

          if (constructor%lmax_ind_pol(j,i) < 0) then
             if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='fullsky') then
                constructor%pol_pixreg_type(j,i) = 1
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='single_pix') then
                constructor%pol_pixreg_type(j,i) = 2
             else if (trim(cpar%cs_spec_pixreg(j,i,id_abs))=='pixreg') then
                constructor%pol_pixreg_type(j,i) = 3
                constructor%npixreg(j,i) = cpar%cs_spec_npixreg(j,i,id_abs) 
             else
                write(*,*) 'Unspecified pixel region type for poltype',j,'of spectral index',i,'in component',id_abs
                stop
             end if
          else
             constructor%pol_pixreg_type(j,i) = 0
          end if
       end do
    end do

    ! Set up spectral index sampling masks, proposal length maps and nprop maps
    allocate(constructor%pol_ind_mask(constructor%npar)) ! masks per spectral index (all poltypes)
    allocate(constructor%pol_nprop(constructor%npar))    ! nprop map per spectral index (all poltypes
    allocate(constructor%pol_proplen(constructor%npar))  ! proplen map per spectral index (all poltypes)
    allocate(constructor%ind_pixreg_map(constructor%npar))   ! pixel region map per spectral index (all poltypes)

    if (any(constructor%pol_pixreg_type(:,:) == 3)) then
       k=0
       do i = 1,constructor%npar
          do j = 1,constructor%poltype(i)
             if (constructor%pol_pixreg_type(j,i) == 3) then
                if (constructor%npixreg(j,i) > k) k = constructor%npixreg(j,i)
             end if
          end do
       end do
       allocate(constructor%theta_pixreg(0:k,3,constructor%npar))
    end if

    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & 3, constructor%pol)

    if (any(constructor%pol_pixreg_type == 3)) then
       if (any(constructor%poltype > 1)) then
          allocate(constructor%ind_pixreg_arr(0:info%np-1,3,constructor%npar)) 
       else
          allocate(constructor%ind_pixreg_arr(0:info%np-1,1,constructor%npar)) 
       end if
       constructor%ind_pixreg_arr = 0 !all pixels assigned to pixelregion 0 (not to be sampled), will read in pixreg later
    end if

    do i = 1,constructor%npar
       ! spec. ind. mask
       if (trim(cpar%cs_spec_mask(i,id_abs)) == 'fullsky') then
          constructor%pol_ind_mask(i)%p => comm_map(info)
          constructor%pol_ind_mask(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          constructor%pol_ind_mask(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_mask(i,id_abs)))

          if (constructor%poltype(i) > constructor%pol_ind_mask(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(constructor%indlabel(i))//' mask has fewer maps (', & 
                  & constructor%pol_ind_mask(i)%p%info%nmaps,') than poltype (',constructor%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (constructor%theta(i)%p%info%nside /= constructor%pol_ind_mask(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(constructor%indlabel(i))//' mask has different nside (', & 
                  & constructor%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & constructor%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if

       ! prop. len. map
       if (trim(cpar%cs_spec_proplen(i,id_abs)) == 'fullsky') then
          constructor%pol_proplen(i)%p => comm_map(info)
          constructor%pol_proplen(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          constructor%pol_proplen(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_proplen(i,id_abs)))
          if (constructor%poltype(i) > constructor%pol_proplen(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(constructor%indlabel(i))//' proplen map has fewer maps (', & 
                  & constructor%pol_ind_mask(i)%p%info%nmaps,') than poltype (',constructor%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (constructor%theta(i)%p%info%nside /= constructor%pol_proplen(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(constructor%indlabel(i))//' proplen map has different nside (', & 
                  & constructor%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & constructor%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if
       
       ! replace proplen of input map for given poltype with user specified value, if given
       do j = 1,constructor%poltype(i)
          if (.not. trim(cpar%cs_spec_proplen_init(j,i,id_abs))=='none') then
             partxt=trim(cpar%cs_spec_proplen_init(j,i,id_abs))
             read(partxt,*) par_dp
             constructor%pol_proplen(i)%p%map(:,j)=par_dp
          end if
       end do
    
       ! nprop map
       if (trim(cpar%cs_spec_nprop(i,id_abs)) == 'fullsky') then
          constructor%pol_nprop(i)%p => comm_map(info)
          constructor%pol_nprop(i)%p%map = 1.d0
       else
          ! Read map from FITS file
          constructor%pol_nprop(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_nprop(i,id_abs)))
          if (constructor%poltype(i) > constructor%pol_nprop(i)%p%info%nmaps) then
             write(*,fmt='(a,i2,a,i2,a,i2)') trim(constructor%indlabel(i))//' nprop map has fewer maps (', & 
                  & constructor%pol_ind_mask(i)%p%info%nmaps,') than poltype (',constructor%poltype(i), &
                  & ') for component nr. ',id_abs
             stop
          else if (constructor%theta(i)%p%info%nside /= constructor%pol_nprop(i)%p%info%nside) then
             write(*,fmt='(a,i4,a,i4,a,i2)') trim(constructor%indlabel(i))//' nprop map has different nside (', & 
                  & constructor%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                  & constructor%theta(i)%p%info%nside, ') for component nr. ',id_abs
             stop
          end if
       end if
       ! replace nprop of input map for given poltype with user specified value, if given       
       do j = 1,constructor%poltype(i)
          if (.not. trim(cpar%cs_spec_nprop_init(j,i,id_abs))=='none') then
             partxt=trim(cpar%cs_spec_nprop_init(j,i,id_abs))
             read(partxt,*) k
             constructor%pol_nprop(i)%p%map(:,j)=k*1.d0
          end if
       end do
       ! limit nprop
       constructor%pol_nprop(i)%p%map = min(max(constructor%pol_nprop(i)%p%map,constructor%nprop_uni(1,i)*1.d0), &
            & constructor%nprop_uni(2,i)*1.d0)

       ! initialize pixel regions if relevant 
       if (any(constructor%pol_pixreg_type(:,i) == 3)) then
          ! Read map from FITS file
          constructor%ind_pixreg_map(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_pixreg_map(i,id_abs)))

          !compute the average theta in each pixel region for the poltype indices that sample theta using pixel regions
          do j =1,constructor%poltype(i)
             constructor%theta_pixreg(:,j,i)=0.d0
             if (.not. constructor%pol_pixreg_type(j,i) == 3) cycle

             n=constructor%npixreg(j,i)
             allocate(sum_pix(n),sum_theta(n))
             sum_theta=0.d0
             sum_pix=0

             do k = 0,constructor%theta(i)%p%info%np-1
                do m = 1,n
                   if ( constructor%ind_pixreg_map(i)%p%map(k,j) > (m-0.5d0) .and. &
                        & constructor%ind_pixreg_map(i)%p%map(k,j) < (m+0.5d0) ) then
                      sum_theta(m)=sum_theta(m)+constructor%theta(i)%p%map(k,j)
                      sum_pix(m)=sum_pix(m)+1
                      constructor%ind_pixreg_arr(k,j,i)=m !assign pixel region index 
                      exit
                   end if
                end do
             end do

             !allreduce
             call mpi_allreduce(MPI_IN_PLACE, sum_theta, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_pix, n, MPI_INTEGER, MPI_SUM, info%comm, ierr)
            
             do k = 1,n
                if (sum_pix(n) > 0) then
                   constructor%theta_pixreg(k,j,i)=sum_theta(k)/(1.d0*sum_pix(k))
                else
                   constructor%theta_pixreg(k,j,i)=constructor%p_gauss(1,i) ! the prior as theta
                end if
             end do
             constructor%theta_pixreg(0,j,i)=constructor%p_gauss(1,i) !all pixels in region 0 has the prior as theta

             deallocate(sum_pix,sum_theta)

             !Should also smooth and assign to theta map.
             if (constructor%poltype(i) == 1) then
                p_min = 1; p_max = constructor%nmaps
                if (cpar%only_pol) p_min = 2
             else if (constructor%poltype(i) == 2) then
                if (j == 1) then
                   p_min = 1; p_max = 1
                else
                   p_min = 2; p_max = constructor%nmaps
                end if
             else if (constructor%poltype(i) == 3) then
                p_min = j
                p_max = j
             else
                write(*,*) 'Unsupported polarization type'
                stop
             end if
             smooth_scale = constructor%smooth_scale(i)
             if (cpar%num_smooth_scales > 0) then
                if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                   if (p_min<=p_max) then !just a guaranty that we dont smooth for nothing
                   
                      info2  => comm_mapinfo(constructor%theta(i)%p%info%comm, cpar%nside_smooth(smooth_scale), &
                           & cpar%lmax_smooth(smooth_scale), 1, constructor%theta(i)%p%info%pol) !only want 1 map

                      !spec. ind. map with 1 map (will be smoothed like zero spin map using the existing code)
                      theta_single_pol => comm_map(info2)

                      do k = 0,info2%np-1
                         theta_single_pol%map(k,1) = constructor%theta_pixreg(constructor%ind_pixreg_arr(k,j,i),j,i)
                      end do
                      !smooth single map as intensity (i.e. zero spin)
                      call smooth_map(info2, .false., &
                           & data(1)%B_postproc(smooth_scale)%p%b_l*0.d0+1.d0, theta_single_pol, &  
                           & data(1)%B_postproc(smooth_scale)%p%b_l, theta_single_pol_smooth)

                      do k = p_min,p_max
                         constructor%theta(i)%p%map(:,k) = theta_single_pol_smooth%map(:,1)
                      end do
                      call theta_single_pol_smooth%dealloc()
                      theta_single_pol_smooth => null()
                      call theta_single_pol%dealloc()
                      theta_single_pol => null()

                   end if
                end if
             end if !num smooth scales > 0

          end do
          
       end if

    end do

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_spindust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band    
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    real(dp) :: scale, nu_p

    nu_p    = theta(1)
    !alpha   = theta(2)
    scale   = self%nu_p0 / (nu_p*1.d9) ! nu_p is in GHz

    if (scale*nu < self%nu_min .or. scale*nu > self%nu_max) then
       evalSED = 0.d0
    else
       evalSED = exp(splint(self%SED_spline, log(scale*nu))) / &
               & exp(splint(self%SED_spline, log(scale*self%nu_ref(pol)))) * (self%nu_ref(pol)/nu)**2
    !           & exp(splint(self%SED_spline, log(scale*self%nu_ref))) * (self%nu_ref/nu)**(2.d0-alpha)
    end if

  end function evalSED
  
end module comm_spindust_comp_mod
