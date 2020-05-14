module comm_physdust_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  use spline_2D_mod
  implicit none

  private
  public comm_physdust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_physdust_comp
     integer(i4b) :: num_nu, num_u, num_comp, num_u_int
     real(dp)     :: log_umax, gamma, alpha
     real(dp)     :: u_arr_min, u_arr_max, du_arr
     real(dp), allocatable, dimension(:)         :: log_dust_wav, dust_u, extcrv, extpol, scat, amps
     real(dp), allocatable, dimension(:,:,:,:,:) :: coeff_arr
   contains
     procedure :: S    => evalSED
  end type comm_physdust_comp

  interface comm_physdust_comp
     procedure constructor
  end interface comm_physdust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_physdust_comp), pointer   :: constructor

    real(dp)     :: dummy
    real(dp), allocatable, dimension(:,:,:) :: comp_mat
    integer(i4b) ::  unit, num_nu
    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix, sum_nprop
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen 
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 


    ! General parameters
    allocate(constructor)


    constructor%npar         = 2
    allocate(constructor%poltype(constructor%npar))
    allocate(constructor%lmax_ind_pol(3,constructor%npar))       ! {integer}: lmax per. poltype sample per spec. index
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
       do j = 1, constructor%poltype(i)
          !assign lmax per spec ind per polarization sample type (poltype)
          constructor%lmax_ind_pol(j,i) = cpar%cs_lmax_ind_pol(j,i,id_abs)
       end do
    end do

    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%indlabel(1))
    constructor%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    constructor%indlabel(1)  = 'Umin'

    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)


    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(1)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_input_ind(i,id_abs)))
       end if

       if (constructor%lmax_ind >= 0) then
          ! if any polarization is local sampled, only use alms to set polarizations with alm sampling
          if (any(constructor%lmax_ind_pol(1:constructor%poltype(i),i) < 0)) then
             tp => comm_map(info)
             tp%alm=constructor%theta(i)%p%alm
             call tp%Y_scalar
             do p = 1,constructor%poltype(i)
                if (constructor%lmax_ind_pol(p,i) < 0) cycle
                if (constructor%poltype(i) == 1) then
                   p_min=1
                   p_max=constructor%nmaps
                   if (cpar%only_pol) p_min = 2
                else if (constructor%poltype(i)==2) then
                   if (p == 1) then
                      p_min = 1
                      p_max = 1
                   else
                      p_min = 2
                      p_max = constructor%nmaps
                   end if
                else if (constructor%poltype(i)==3) then
                   p_min = p
                   p_max = p
                else
                   write(*,*) '  Unknown poltype in component ',id_abs,', parameter ',i 
                   stop
                end if

                do j = p_min,p_max
                   constructor%theta(i)%p%map(:,j) = tp%map(:,p)
                end do
             end do
             call tp%dealloc()
          else
             call constructor%theta(i)%p%YtW_scalar
          end if
       end if
    end do

    ! Initialize SED templates
    unit                  = getlun()
    num_nu                = 1000    ! Number of frequencies in dust files
    constructor%num_nu    = num_nu
    constructor%num_u     = 11      ! Number of radiation field strengths
    constructor%u_arr_min = -0.5d0
    constructor%u_arr_max =  0.5d0
    constructor%du_arr    = (constructor%u_arr_max - constructor%u_arr_min)/(constructor%num_u - 1.d0)
    constructor%num_comp  = 4 ! Number of dust components

    allocate(constructor%coeff_arr(4,4,num_nu,constructor%num_u,constructor%num_comp))
    allocate(constructor%log_dust_wav(num_nu), constructor%dust_u(constructor%num_u))
    allocate(constructor%extcrv(num_nu), constructor%extpol(num_nu), constructor%scat(num_nu))
    allocate(constructor%amps(constructor%num_comp))

    constructor%log_umax  = cpar%cs_auxpar(1,id_abs)
    constructor%gamma     = cpar%cs_auxpar(2,id_abs)
    constructor%alpha     = cpar%cs_auxpar(3,id_abs)
    constructor%amps      = cpar%cs_auxpar(4:3+constructor%num_comp,id_abs)

    ! Make array of log U
    do i = 1, constructor%num_u
       constructor%dust_u(i) = -0.5d0 + constructor%du_arr*(i-1)
    enddo

    allocate(comp_mat(num_nu,constructor%num_u,constructor%num_comp))    
    do k = 1, constructor%num_comp
       open(unit,file=trim(cpar%datadir)//'/'//trim(cpar%cs_SED_template(k,id_abs)),recl=4096)
       do i = 1, num_nu
          read(unit,fmt='(1P26E11.3)') constructor%log_dust_wav(i), constructor%extcrv(i),&
               & constructor%extpol(i), constructor%scat(i), (comp_mat(i,j,k),j=1,constructor%num_u), &
                (dummy,j=1,constructor%num_u)
          constructor%log_dust_wav(i) = log(constructor%log_dust_wav(i))
       enddo
       call splie2_full_precomp(constructor%log_dust_wav,constructor%dust_u,&
            & log(comp_mat(:,:,k)), constructor%coeff_arr(:,:,:,:,k))
        close(59)
     enddo
     deallocate(comp_mat)


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

    ! Set up smoothing scale information
    allocate(constructor%smooth_scale(constructor%npar))
    constructor%smooth_scale = cpar%cs_smooth_scale(id_abs,1:constructor%npar)


    !!! (local) sampling specific parameters!!!

    allocate(constructor%pol_lnLtype(3,constructor%npar))        ! {chisq, ridge, marginal}: evaluation type for lnL
    allocate(constructor%pol_sample_nprop(3,constructor%npar))   ! {.true., .false.}: sample nprop on first iteration
    allocate(constructor%pol_sample_proplen(3,constructor%npar)) ! {.true., .false.}: sample proplen on first iteration
    allocate(constructor%pol_pixreg_type(3,constructor%npar))    ! {1=fullsky, 2=single_pix, 3=pixel_regions}
    allocate(constructor%nprop_uni(2,constructor%npar))          ! {integer}: upper and lower limits on nprop
    allocate(constructor%npixreg(3,constructor%npar))            ! {integer}: number of pixel regions per poltye per spec ind
    allocate(constructor%first_ind_sample(3,constructor%npar)) !used for pixelregion sampling
    constructor%first_ind_sample=.true.

    constructor%npixreg = 0
    do i = 1,constructor%npar
       constructor%nprop_uni(:,i)=cpar%cs_spec_uni_nprop(:,i,id_abs)
       do j = 1,constructor%poltype(i)
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

    if (any(constructor%pol_pixreg_type(:,:) == 3)) then
       k = 0
       do i = 1,constructor%npar
          do j = 1,constructor%poltype(i)
             if (constructor%npixreg(j,i) > k) k = constructor%npixreg(j,i)
          end do
       end do
       allocate(constructor%nprop_pixreg(k,3,constructor%npar))
       allocate(constructor%npix_pixreg(k,3,constructor%npar))
       allocate(constructor%proplen_pixreg(k,3,constructor%npar))
       allocate(constructor%B_pp_fr(constructor%npar))
    end if

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
         & constructor%nmaps, constructor%pol)

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
          if (j > constructor%nmaps) cycle
          if (cpar%cs_spec_proplen_init(j,i,id_abs) > 0.d0) then
             constructor%pol_proplen(i)%p%map(:,j)=cpar%cs_spec_proplen_init(j,i,id_abs)
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
          if (j > constructor%nmaps) cycle
          if (cpar%cs_spec_nprop_init(j,i,id_abs) > 0) then
             constructor%pol_nprop(i)%p%map(:,j)=cpar%cs_spec_nprop_init(j,i,id_abs)*1.d0
          end if
       end do

       ! limit nprop
       constructor%pol_nprop(i)%p%map = min(max(constructor%pol_nprop(i)%p%map,constructor%nprop_uni(1,i)*1.d0), &
            & constructor%nprop_uni(2,i)*1.d0)

       ! initialize pixel regions if relevant 
       if (any(constructor%pol_pixreg_type(1:constructor%poltype(i),i) == 3)) then

          info2  => comm_mapinfo(constructor%theta(i)%p%info%comm, constructor%nside, &
               & 2*constructor%nside, 1, .false.) 

          smooth_scale = constructor%smooth_scale(i)
          if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
             if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                !create beam for smoothing
                constructor%B_pp_fr(i)%p => comm_B_bl(cpar, info2, 1, 1, fwhm=cpar%fwhm_smooth(smooth_scale), &
                     & init_realspace=.false.)
             end if
          end if

          ! Check if 'fullsky' or 'none' 
          if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'fullsky') then
             constructor%ind_pixreg_map(i)%p => comm_map(info)
             constructor%ind_pixreg_map(i)%p%map = 1.d0
          else if (trim(cpar%cs_spec_pixreg_map(i,id_abs)) == 'none') then
             constructor%ind_pixreg_map(i)%p => comm_map(info)
             constructor%ind_pixreg_map(i)%p%map = 0.d0
          else
             ! Read map from FITS file
             constructor%ind_pixreg_map(i)%p => comm_map(info, trim(cpar%datadir) // '/' // trim(cpar%cs_spec_pixreg_map(i,id_abs)))
             if (constructor%poltype(i) > constructor%ind_pixreg_map(i)%p%info%nmaps) then
                write(*,fmt='(a,i2,a,i2,a,i2)') trim(constructor%indlabel(i))//' pixreg map has fewer maps (', & 
                     & constructor%pol_ind_mask(i)%p%info%nmaps,') than poltype (',constructor%poltype(i), &
                     & ') for component nr. ',id_abs
                stop
             else if (constructor%theta(i)%p%info%nside /= constructor%ind_pixreg_map(i)%p%info%nside) then
                write(*,fmt='(a,i4,a,i4,a,i2)') trim(constructor%indlabel(i))//' pixreg map has different nside (', & 
                     & constructor%pol_ind_mask(i)%p%info%nside,') than the spectral index map (', &
                     & constructor%theta(i)%p%info%nside, ') for component nr. ',id_abs
                stop
             end if
          end if

          !compute the average theta in each pixel region for the poltype indices that sample theta using pixel regions
          do j = 1,constructor%poltype(i)
             if (j > constructor%nmaps) cycle
             constructor%theta_pixreg(:,j,i)=constructor%p_gauss(1,i) !prior
             if (.not. constructor%pol_pixreg_type(j,i) == 3) cycle

             n=constructor%npixreg(j,i)
             allocate(sum_pix(n),sum_theta(n),sum_nprop(n),sum_proplen(n))
             sum_theta=0.d0
             sum_pix=0
             sum_nprop=0
             sum_proplen=0

             do k = 0,constructor%theta(i)%p%info%np-1
                do m = 1,n
                   if ( constructor%ind_pixreg_map(i)%p%map(k,j) > (m-0.5d0) .and. &
                        & constructor%ind_pixreg_map(i)%p%map(k,j) < (m+0.5d0) ) then
                      sum_theta(m)=sum_theta(m)+constructor%theta(i)%p%map(k,j)
                      sum_proplen(m)=sum_proplen(m)+constructor%pol_proplen(i)%p%map(k,j)
                      sum_nprop(m)=sum_nprop(m)+IDINT(constructor%pol_nprop(i)%p%map(k,j))
                      sum_pix(m)=sum_pix(m)+1
                      constructor%ind_pixreg_arr(k,j,i)=m !assign pixel region index 
                      exit
                   end if
                end do
             end do

             !allreduce
             call mpi_allreduce(MPI_IN_PLACE, sum_theta, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_proplen, n, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_nprop, n, MPI_INTEGER, MPI_SUM, info%comm, ierr)
             call mpi_allreduce(MPI_IN_PLACE, sum_pix, n, MPI_INTEGER, MPI_SUM, info%comm, ierr)
            
             do k = 1,n
                !if (cpar%myid == cpar%root) write(*,*) 'pixreg',k,'  -- numbe of pixels',sum_pix(k)
                if (sum_pix(k) > 0) then
                   constructor%theta_pixreg(k,j,i)=sum_theta(k)/(1.d0*sum_pix(k))
                   constructor%nprop_pixreg(k,j,i)=IDINT(sum_nprop(k)/(1.d0*sum_pix(k)))
                   constructor%proplen_pixreg(k,j,i)=sum_proplen(k)/(1.d0*sum_pix(k))
                else
                   constructor%theta_pixreg(k,j,i)=constructor%p_gauss(1,i) ! the prior as theta
                   constructor%nprop_pixreg(k,j,i)=0
                   constructor%proplen_pixreg(k,j,i)=1.d0
                end if
                constructor%npix_pixreg(k,j,i)=sum_pix(k)
             end do
             constructor%theta_pixreg(0,j,i)=constructor%p_gauss(1,i) !all pixels in region 0 has the prior as theta

             deallocate(sum_pix,sum_theta,sum_proplen,sum_nprop)
             

             !Should also assign (and smooth) theta map. This might be dangerous so we don't assign for now
             ! might put in a flag for this if wanted
             if (.false.) then
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
                if (cpar%num_smooth_scales > 0 .and. smooth_scale > 0) then
                   if (cpar%fwhm_postproc_smooth(smooth_scale) > 0.d0) then
                      if (p_min<=p_max) then !just a guaranty that we dont smooth for nothing

                         !spec. ind. map with 1 map (will be smoothed like zero spin map using the existing code)
                         tp => comm_map(info2)

                         do k = 0,info2%np-1
                            tp%map(k,1) = constructor%theta_pixreg(constructor%ind_pixreg_arr(k,j,i),j,i)
                         end do
                         !smooth single map as intensity (i.e. zero spin)
                         call smooth_map(info2, .false., &
                              & constructor%B_pp_fr(i)%p%b_l*0.d0+1.d0, tp, &  
                              & constructor%B_pp_fr(i)%p%b_l, tp_smooth)

                         do k = p_min,p_max
                            constructor%theta(i)%p%map(:,k) = tp_smooth%map(:,1)
                         end do
                         call tp_smooth%dealloc()
                         call tp%dealloc()

                      end if
                   end if
                end if !num smooth scales > 0
             end if
          end do !poltype
          
       end if !any pixreg_type == 3

    end do !npar

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_physdust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, j, num_u_int
    real(dp) :: SED_physdust, SED_norm, wav_in, wav_ref, c
    real(dp) :: log_umin, umin, log_umax, umax, uval, du, fdu

    if (nu < 2d9) then
       evalSED = 0.d0
       return
    end if
     
    num_u_int = 100 ! How many dU steps to perform
    c         = 2.99792458d8
     
    log_umin     = theta(1)
    log_umax     = self%log_umax
    SED_physdust = 0.d0
    SED_norm     = 0.d0
    wav_in       = c/nu               * 1d6 ! In um
    wav_ref      = c/self%nu_ref(pol) * 1d6 ! In um
    umin         = 10**log_umin
    umax         = 10**log_umax
     do i = 1, self%num_comp
        ! Delta function component
        SED_physdust = SED_physdust + (1.d0-self%gamma)* &
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,&
             & self%coeff_arr(:,:,:,:,i), log(wav_in), log_umin)))
        SED_norm = SED_norm + (1.d0-self%gamma)*&
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
             log(wav_ref),log_umin)))

        ! Distribution of U values
        ! See, e.g., Aniano et al 2012, Eqs. 9 & 10
        if (self%gamma /= 0.d0) then
           du = umin*((umax/umin)**(1.d0/(num_u_int-1.d0)) - 1d0) ! Formally dU/U
           do j = 1, num_u_int ! integrate over U distribution
              uval = umin*(umax/umin)**((j-1.d0)/(num_u_int-1.d0))
              if (self%alpha /= 1.d0) then
                 fdu = uval**(1.d0-self%alpha)*du*self%gamma*&
                      & (self%alpha-1.d0)/(umin**(1.d0-self%alpha) - umax**(1.d0-self%alpha))
              else
                 ! Special case of alpha = 1 to ensure integrability
                 fdu = du*self%gamma/log(umax/umin)
              endif
              SED_physdust = SED_physdust + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   & log(wav_in),log10(uval)))*fdu
              SED_norm = SED_norm + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   log(wav_ref),log10(uval)))*fdu
           enddo
        endif
     enddo
     evalSED = (SED_physdust / SED_norm) * (self%nu_ref(pol)/nu)**3 ! Normalize to reference in T units

  end function evalSED
  
end module comm_physdust_comp_mod
