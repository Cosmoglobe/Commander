module comm_data_mod
  use comm_utils
  use comm_fg_component_mod
  use comm_genvec_mod
  use comm_N_mult_mod
  use comm_beam_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2013, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b)                                    :: nside, npix, lmax, nmaps, numcomp, nspec
  integer(i4b)                                    :: nu_fullsky, nu_highlat
  integer(i4b),                           private :: numproc_per_band
  integer(i4b),                           private :: myid_data, comm_data, root
  integer(i4b),                           private :: myid_alms, comm_alms
  integer(i4b),                           private :: myid_chain, comm_chain
  type(planck_rng),                       private :: rng_handle, noise_handle
  logical(lgt)                                    :: polarization
  logical(lgt)                                    :: sample_T_modes
  integer(i4b)                                    :: map_size, map_id, lmax_corr, nside_corr
  integer(i4b)                                    :: num_fg_temp
  integer(i4b)                                    :: num_fg_signal, lmax_lowres
  logical(lgt)                                    :: sample_fg_temp, sample_fg_pix
  logical(lgt)                                    :: sample_temp_coeffs=.true.
  logical(lgt)                                    :: exclude_pos_fg_amp
  logical(lgt),                           private :: output_cmb_freq_map, output_mixmat
  logical(lgt),                           private :: enforce_zero_cl
  real(dp)                                        :: fwhm_lowres, corr_chisq_threshold
  real(dp),                               private :: reg_noise, reg_scale
  character(len=128),                     private :: operation

  real(dp),     allocatable, dimension(:,:)       :: cmbmap, residual, foreground_map, tempamp, mask_lowres
  real(dp),     allocatable, dimension(:,:,:)     :: cmbmaps_lowres, residuals_lowres
  real(dp),     allocatable, dimension(:,:,:,:)   :: fg_temp_lowres
  real(dp),     allocatable, dimension(:,:)       :: signal_map, mask_corr, mask_calib
  real(dp),     allocatable, dimension(:,:,:)     :: fg_pix_spec_response
  real(dp),     allocatable, dimension(:,:,:)     :: fg_temp
  real(dp),     allocatable, dimension(:,:,:)     :: inv_N_lowres, inv_N_scaled
  logical(lgt), allocatable, dimension(:,:)       :: mask
  integer(i4b), pointer,     dimension(:)         :: pixels

contains

  subroutine initialize_data_mod(paramfile, comm_chain_in, comm_data_in, comm_alms_in, &
       & seed, map_id_in)
    implicit none

    character(len=128),                intent(in)           :: paramfile
    integer(i4b),                      intent(in)           :: comm_data_in, comm_alms_in
    integer(i4b),                      intent(in)           :: comm_chain_in, seed
    integer(i4b),                      intent(in), optional :: map_id_in

    integer(i4b)       :: ordering, temp_i, num_md_temp, i, j, k, l, ierr, ind
    integer(i4b)       :: myid_world, map_set, tot_numpix, unit
    logical(lgt)       :: bad, sample_foreground_templates, exist, sample_inside_mask
    real(dp)           :: S_eff, fwhm_degrade, my_reg_noise, my_reg_scale
    character(len=1)   :: i1
    character(len=2)   :: i2
    character(len=3)   :: i3
    character(len=2)   :: map_text, temp_text
    character(len=4)   :: nside_text, real_text
    character(len=128) :: maskfile, maskfile_calib, cmbfile, rmsfile, filename
    character(len=128) :: paramtext, maskfile_corr
    character(len=512) :: temp_amp_file
    character(len=32)  :: label
    real(dp)           :: nu_eff
    real(dp),           allocatable, dimension(:)     :: temp_map
    real(dp),           allocatable, dimension(:,:)   :: noisemap, noisemap_in, map_in, map
    real(dp),           allocatable, dimension(:,:)   :: matrix, matrix2, alms1
    real(dp),           allocatable, dimension(:,:)   :: foreground_map_in, map_lowres!, my_beam
    real(dp),           pointer,     dimension(:,:)   :: pixwin
    logical(lgt),       allocatable, dimension(:,:)   :: mask_in
    character(len=128), allocatable, dimension(:)     :: temp_names
    integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

    unit        = getlun()
    root        = 0
    comm_data = comm_data_in
    call mpi_comm_rank(comm_data, myid_data, ierr)
    call mpi_comm_size(comm_data, numband, ierr)
    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    call mpi_comm_size(comm_alms, numproc_per_band, ierr)
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)

    call mpi_comm_rank(MPI_COMM_WORLD, myid_world, ierr)

    call rand_init(rng_handle, seed)

    if (myid_alms == root) then

       ! Read general parameters
       call int2string(map_id_in, map_text)
       call get_parameter(paramfile, 'NSIDE',                     par_int=nside)
       call get_parameter(paramfile, 'LMAX',                      par_int=lmax)
       call get_parameter(paramfile, 'POLARIZATION',              par_lgt=polarization)
       call get_parameter(paramfile, 'MASKFILE',                  par_string=maskfile)
       call get_parameter(paramfile, 'MASKFILE_CALIB'//map_text,  par_string=maskfile_calib)
       call get_parameter(paramfile, 'OUTPUT_CMB_FREQUENCY_MAPS', par_lgt=output_cmb_freq_map)
       call get_parameter(paramfile, 'ENFORCE_ZERO_CL',           par_lgt=enforce_zero_cl)
       call get_parameter(paramfile, 'OPERATION',                 par_string=operation)
       call get_parameter(paramfile, 'CORR_CHISQ_THRESHOLD',      par_dp=corr_chisq_threshold)
       call get_parameter(paramfile, 'LMAX_CORR',                 par_int=lmax_corr)
       call get_parameter(paramfile, 'MASKFILE_CORR',             par_string=maskfile_corr)
       call get_parameter(paramfile, 'NSIDE_CORR',                par_int=nside_corr)
       call get_parameter(paramfile, 'TEMPLATE_AMP_INPUT',        par_string=temp_amp_file)
       call get_parameter(paramfile, 'OUTPUT_MIXING_MATRIX',      par_lgt=output_mixmat)
       call get_parameter(paramfile, 'REGULARIZATION_NOISE',      par_dp=reg_noise)
       call get_parameter(paramfile, 'REG_NOISE_SCALE'//map_text, par_dp=reg_scale)
       call get_parameter(paramfile, 'SAMPLE_INSIDE_MASK',        par_lgt=sample_inside_mask)
       reg_noise = reg_noise * reg_scale

       map_id             = map_id_in
       npix               = 12*nside**2
       numcomp            = lmax2numcomp(lmax)
       if (polarization) then
          nmaps = 3
          nspec = 6
       else
          nmaps = 1
          nspec = 1
       end if

       if (polarization) then
          call get_parameter(paramfile, 'SAMPLE_TT_MODES',  par_lgt=sample_T_modes)
       else
          sample_T_modes = .true.
       end if

       ! Read foreground information
       call get_parameter(paramfile, 'NUM_FG_TEMP', &
            & par_int=num_fg_temp)
       sample_fg_temp = num_fg_temp > 0
       if (sample_fg_temp) then
          allocate(temp_names(num_fg_temp))

          ! Read foreground template filenames
          do j = 1, num_fg_temp
             call int2string(j,temp_text)
             paramtext = 'FG_TEMPLATE' // map_text // '_' // temp_text
             call get_parameter(paramfile, trim(paramtext), &
                  & par_string=temp_names(j))
          end do

       else
          num_fg_temp = 0
       end if

       if (sample_T_modes) then
          num_md_temp        = 4
          num_fg_temp   = num_fg_temp + num_md_temp
       else
          num_md_temp        = 0
       end if
       
       call get_parameter(paramfile, 'NUM_FG_SIGNAL_COMPONENTS', par_int=num_fg_signal)
       sample_fg_pix = num_fg_signal > 0

       allocate(mask_corr(0:npix-1,1))
       call read_map(maskfile_corr, mask_corr)
       where (mask_corr > 0.5)
          mask_corr = 1.d0
       elsewhere
          mask_corr = 0.d0
       end where

    end if

    call mpi_bcast(map_id,              1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(polarization,        1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(nside,               1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(npix,                1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(lmax,                1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(nspec,               1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(numcomp,             1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(nmaps,               1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(sample_T_modes,      1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(enforce_zero_cl,     1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(maskfile,            len(maskfile),      MPI_CHARACTER, root, comm_alms, ierr)    
    call mpi_bcast(maskfile_calib,      len(maskfile_calib), MPI_CHARACTER, root, comm_alms, ierr)    
    call mpi_bcast(operation,           len(operation),     MPI_CHARACTER, root, comm_alms, ierr)    
    call mpi_bcast(temp_amp_file,       len(temp_amp_file), MPI_CHARACTER, root, comm_alms, ierr)    
    call mpi_bcast(num_fg_temp,    1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(num_md_temp,         1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(sample_fg_temp,      1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(sample_fg_pix,       1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(num_fg_signal,       1,   MPI_INTEGER,   root, comm_alms, ierr)    
    call mpi_bcast(output_cmb_freq_map, 1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(output_mixmat,       1,   MPI_LOGICAL,   root, comm_alms, ierr)    
    call mpi_bcast(reg_noise,           1,   MPI_DOUBLE_PRECISION,root, comm_alms, ierr)    
    call mpi_bcast(sample_inside_mask,  1,   MPI_LOGICAL,   root, comm_alms, ierr)    

    call get_pixels(pixels)
    map_size = size(pixels)

    allocate(cmbmap(0:map_size-1, nmaps))
    allocate(residual(0:map_size-1, nmaps))
    allocate(signal_map(0:map_size-1, nmaps))
    allocate(foreground_map(0:map_size-1, nmaps))
    allocate(mask(0:map_size-1, nmaps))
    allocate(mask_calib(0:map_size-1, nmaps))
    allocate(mask_lowres(0:npix-1, nmaps))
       
    foreground_map = 0.d0

    ! Read mask file
    allocate(mask_in(0:npix-1, nmaps), map_in(0:npix-1,nmaps))
    call read_map(maskfile, map_in)
    where (map_in > 0.5d0) 
       mask_lowres = 1.d0
       mask_in     = .true.
    elsewhere
       mask_lowres = 0.d0
       mask_in     = .false.
    end where
    nu_highlat  = numband * count(map_in > 0.5d0) 
    if (sample_inside_mask) then
       nu_fullsky  = numband * npix
    else
       nu_fullsky = nu_highlat
    end if
    mask        = mask_in(pixels,:)

    call read_map(maskfile_calib, map_in)
    mask_calib = map_in(pixels,:)
    where (mask_calib > 0.5d0)
       mask_calib = 1.d0
    elsewhere
       mask_calib = 0.d0
    end where
    deallocate(map_in, mask_in)


    ! Initialize templates
    if (num_fg_temp > 0) then

       allocate(fg_temp(0:map_size-1,nmaps,num_fg_temp))

       fg_temp = 0.d0
       if (sample_T_modes) then
          call initialize_mono_and_dipole_temp(nside, fg_temp(:,1,1:4), pixels)
       end if

       if (sample_fg_temp) then
          allocate(foreground_map_in(0:npix-1,nmaps))
          
          if (myid_alms == root) then
                
             do j = 1, num_fg_temp-num_md_temp
                call read_map(temp_names(j), foreground_map_in)
                if (nmaps==3 .and. (.not. sample_T_modes)) foreground_map_in(:,1) = 0.d0
                call mpi_bcast(foreground_map_in, size(foreground_map_in), &
                     & MPI_DOUBLE_PRECISION, root, comm_alms, ierr)
                fg_temp(:,:,num_md_temp+j) = foreground_map_in(pixels,:)
             end do
             
             deallocate(temp_names)
             
          else
             
             do j = 1, num_fg_temp-num_md_temp
                call mpi_bcast(foreground_map_in, size(foreground_map_in), &
                     & MPI_DOUBLE_PRECISION, root, comm_alms, ierr)
                fg_temp(:,:,num_md_temp+j) = foreground_map_in(pixels,:)
             end do
             
          end if
          deallocate(foreground_map_in)

          if (polarization .and. (.not. sample_T_modes)) fg_temp(:,1,:) = 0.d0

       end if

    end if

    if (sample_fg_pix) allocate(fg_pix_spec_response(0:map_size-1, nmaps, num_fg_comp))

    ! Get data and parameters for low-resolution spectral index sampling; root only
    allocate(fg_temp_lowres(0:npix-1, nmaps, numband, num_fg_temp))
    if (myid_chain == root) then

       call get_parameter(paramfile, 'FWHM_LOWRES', par_dp=fwhm_lowres)
       call get_parameter(paramfile, 'LMAX_LOWRES', par_int=lmax_lowres)

       allocate(cmbmaps_lowres(0:npix-1, nmaps, numband))
       allocate(inv_N_lowres(0:npix-1, nmaps, numband))
       allocate(inv_N_scaled(0:npix-1, nmaps, numband))
       allocate(residuals_lowres(0:npix-1, nmaps, numband))

       fg_temp_lowres = 0.d0
       inv_N_lowres  = 0.d0
       do i = 1, numband
          call int2string(i, map_text)

          paramtext = 'NOISE_RMS' // map_text
          call get_parameter(paramfile, trim(paramtext), par_string=rmsfile)
          paramtext = 'REG_NOISE_SCALE' // map_text
          call get_parameter(paramfile, trim(paramtext), par_dp=my_reg_scale)
          paramtext = 'REGULARIZATION_NOISE'
          call get_parameter(paramfile, trim(paramtext), par_dp=my_reg_noise)
          my_reg_noise = my_reg_noise * my_reg_scale

          allocate(map_in(0:npix-1,nmaps))
          call read_map(rmsfile, map_in)
          if (my_reg_noise > 0.d0) then
             map_in = my_reg_noise
          else if (my_reg_noise < 0.d0) then
             map_in = sqrt(map_in**2 + my_reg_noise**2)
          end if
          do j = 0, npix-1
             do k = 1, nmaps
                if ((mask_lowres(j,k) > 0.5d0 .or. sample_inside_mask) .and. map_in(j,k) > 0.d0) then
                   inv_N_lowres(j,k,i) = 1.d0 / map_in(j,k)**2
                end if
             end do
          end do
          deallocate(map_in)

          if (sample_T_modes) then
             call initialize_mono_and_dipole_temp(nside, fg_temp_lowres(:,1,i,1:4))
          end if
       
          if (num_fg_temp > num_md_temp) then
             allocate(map_in(0:npix-1,nmaps))

             ! Read foreground templates
             do j = 1, num_fg_temp-num_md_temp
                call int2string(j,temp_text)
                paramtext = 'FG_TEMPLATE' // map_text // '_' // temp_text
                call get_parameter(paramfile, trim(paramtext), par_string=filename)
                call read_map(filename, map_in)
                if (nmaps==3 .and. (.not. sample_T_modes)) map_in(:,1) = 0.d0
                fg_temp_lowres(:,:,i,num_md_temp+j) = map_in
             end do

             deallocate(map_in)
          end if

        end do

     endif

     do j = 1, numband
        do i = 1, num_fg_temp
           call mpi_bcast(fg_temp_lowres(:,:,j,i), size(fg_temp_lowres(:,:,j,i)),   MPI_DOUBLE_PRECISION,  &
                & root, comm_chain, ierr)
        end do
     end do

     if (sample_fg_pix) then
        if (.not. allocated(inv_N_lowres)) allocate(inv_N_lowres(0:npix-1, nmaps, numband))
        if (.not. allocated(inv_N_scaled)) allocate(inv_N_scaled(0:npix-1, nmaps, numband))
        call mpi_bcast(inv_N_lowres, size(inv_N_lowres), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
        inv_N_scaled = inv_N_lowres
     end if

     allocate(tempamp(num_fg_temp,numband))
     inquire(file=trim(temp_amp_file), exist=exist)
     if (exist) then
        open(unit, file=trim(temp_amp_file))
        do i = 1, numband
           read(unit,*) label, tempamp(:,i)
        end do
        close(unit)
     else
        tempamp = 0.d0
     end if

  end subroutine initialize_data_mod


  subroutine read_map_realization(paramfile, realization_id)
    implicit none

    character(len=*), intent(in) :: paramfile
    integer(i4b),     intent(in) :: realization_id

    integer(i4b)       :: ierr, i, j, k
    character(len=2)   :: map_text
    character(len=4)   :: realization_text
    character(len=128) :: cmbfile, paramtext
    real(dp)           :: my_reg_scale, my_reg_noise
    real(dp), allocatable, dimension(:,:) :: map_in, weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms

    call int2string(map_id, map_text)
    call int2string(realization_id, realization_text)
    
    ! Get filenames
    paramtext = 'MAP' // map_text // '_' // realization_text
    call get_parameter(paramfile, trim(paramtext), par_string=cmbfile)

    ! Read the CMB data
    allocate(map_in(0:npix-1,nmaps))
    if (myid_alms == root) then
       call read_map(cmbfile, map_in)
       call rand_init(noise_handle, map_id)
       do j = 1, nmaps
          do i = 0, npix-1
             map_in(i,j) = map_in(i,j) + reg_noise * rand_gauss(noise_handle)
          end do
       end do
    end if
    call mpi_bcast(map_in, size(map_in), MPI_DOUBLE_PRECISION, root, comm_alms, ierr)
    cmbmap = map_in(pixels,:)
    if (polarization .and. (.not. sample_T_modes)) cmbmap(:,1) = 0.d0
    deallocate(map_in)

    if (myid_chain == root) then
       if (num_fg_comp > 0) then
          allocate(map_in(0:npix-1,nmaps))
          do i = 1, numband
             call int2string(i, map_text)
             paramtext = 'REG_NOISE_SCALE' // map_text
             call get_parameter(paramfile, trim(paramtext), par_dp=my_reg_scale)
             paramtext = 'REGULARIZATION_NOISE'
             call get_parameter(paramfile, trim(paramtext), par_dp=my_reg_noise)
             my_reg_noise = my_reg_noise * my_reg_scale
             paramtext = 'MAP' // map_text // '_' // realization_text
             call get_parameter(paramfile, trim(paramtext), par_string=cmbfile)
             call read_map(cmbfile, map_in)
             call rand_init(noise_handle, i)
             if (my_reg_noise /= 0.d0) then
                do k = 1, nmaps
                   do j = 0, npix-1
                      map_in(j,k) = map_in(j,k) + my_reg_noise * rand_gauss(noise_handle)
                   end do
                end do
             end if
             cmbmaps_lowres(:,:,i) = map_in       
!             do j = 1, num_fg_temp
!                if (fix_temp(j,i)) then
!                   cmbmaps_lowres(:,:,i) = cmbmaps_lowres(:,:,i) - tempamp(j,i) * fg_temp_lowres(:,:,i,j)
!                end if
!             end do
          end do
          deallocate(map_in)
       end if
    end if

!!$    if (allocated(tempamp)) then
!!$       do i = 1, num_fg_temp
!!$          if (fix_temp(i,map_id)) then
!!$             cmbmap = cmbmap - tempamp(i,map_id) * fg_temp(:,:,i)
!!$          end if
!!$       end do
!!$    end if

  end subroutine read_map_realization


  subroutine cleanup_data_mod
    implicit none

    if (allocated(cmbmap))                deallocate(cmbmap)
    if (allocated(signal_map))            deallocate(signal_map)
    if (allocated(residual))              deallocate(residual)
    if (allocated(foreground_map))        deallocate(foreground_map)
    if (allocated(fg_temp))               deallocate(fg_temp)
    if (allocated(mask))                  deallocate(mask)
    if (allocated(fg_pix_spec_response))  deallocate(fg_pix_spec_response)
    if (allocated(cmbmaps_lowres))        deallocate(cmbmaps_lowres)
    if (allocated(inv_N_lowres))          deallocate(inv_N_lowres)
    if (allocated(inv_N_scaled))          deallocate(inv_N_scaled)
    if (allocated(residuals_lowres))      deallocate(residuals_lowres)
    if (allocated(fg_temp_lowres))        deallocate(fg_temp_lowres)
    if (allocated(mask_calib))            deallocate(mask_calib)
    deallocate(pixels)

  end subroutine cleanup_data_mod



  subroutine compute_lowres_residual(s)
    implicit none

    type(genvec), intent(in) :: s

    integer(i4b) :: i, j, l, m, ierr, numpix
    real(dp)     :: avg
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: gb, map
    real(dp),     pointer,     dimension(:,:), save   :: pixwin

    residuals_lowres = cmbmaps_lowres

!!$    if (.not. enforce_zero_cl) then
!!$       ! Subtract the CMB part
!!$       allocate(gb(0:lmax,nmaps), alms(nmaps,0:lmax_lowres,0:lmax_lowres), map(0:npix-1,nmaps))
!!$       call gaussbeam(fwhm_lowres, lmax_lowres, gb)
!!$       if (.not. (associated(pixwin))) call read_pixwin(nside, nmaps, pixwin)
!!$       call convert_real_to_complex_alms(s%cmb_amp, alms)
!!$       
!!$       do i = 1, nmaps
!!$          do l = 0, lmax_lowres
!!$             alms(i,l,0:l) = alms(i,l,0:l) * gb(l,i) * pixwin(l,i)
!!$          end do
!!$       end do
!!$       
!!$       if (nmaps == 1) then
!!$          call alm2map(nside, lmax_lowres, lmax_lowres, alms, map(:,1))
!!$       else
!!$          call alm2map(nside, lmax_lowres, lmax_lowres, alms, map)
!!$       end if
!!$
!!$       do i = 1, numband
!!$          residuals_lowres(:,:,i) = residuals_lowres(:,:,i) - spec2data(i, cmb=.true.) * map
!!$       end do
!!$    end if

    ! Subtract the foreground templates
    do i = 1, numband
       do j = 1, num_fg_temp
          residuals_lowres(:,:,i) = residuals_lowres(:,:,i) - s%temp_amp(j,i) * fg_temp_lowres(:,:,i,j)
       end do
    end do

    !deallocate(alms, gb, map)

  end subroutine compute_lowres_residual

  subroutine set_pix_cmb_equal_to_Cl_cmb(s)
    implicit none

    type(genvec), intent(inout) :: s

    integer(i4b) :: i, j, l, m, ierr, numpix, comp
    real(dp)     :: avg
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: gb, map
    real(dp),     pointer,     dimension(:,:), save   :: pixwin

    ! Check if we have a CMB component
    do comp = 1, num_fg_comp
       if (trim(fg_components(comp)%type) == 'cmb') exit
    end do
    if (comp > num_fg_comp) return

    ! If so, compute the smoothed CMB field
    allocate(gb(0:lmax,nmaps), alms(nmaps,0:lmax_lowres,0:lmax_lowres), map(0:npix-1,nmaps))
    call gaussbeam(fwhm_lowres, lmax_lowres, gb)
    if (.not. (associated(pixwin))) call read_pixwin(nside, nmaps, pixwin)
    call convert_real_to_complex_alms(s%cmb_amp, alms)
    
    do i = 1, nmaps
       do l = 0, lmax_lowres
          alms(i,l,0:l) = alms(i,l,0:l) * gb(l,i) * pixwin(l,i)
       end do
    end do

    if (nmaps == 1) then
       call alm2map(nside, lmax_lowres, lmax_lowres, alms, map(:,1))
    else
       call alm2map(nside, lmax_lowres, lmax_lowres, alms, map)
    end if

    where (mask_lowres > 0.5d0)
       s%fg_amp(:,:,comp) = map / compute_ant2thermo_single(fg_components(comp)%nu_ref)
    end where

    deallocate(alms, gb, map)

  end subroutine set_pix_cmb_equal_to_Cl_cmb


  subroutine compute_residual(s, subtract_signal)
    implicit none

    type(genvec), intent(in),    optional :: s
    logical(lgt), intent(in),    optional :: subtract_signal

    integer(i4b) :: i, j, ierr
    logical(lgt) :: sub_signal
    character(len=2) :: band
    real(dp), allocatable, dimension(:,:) :: a, my_coeff
 
    if (myid_chain == root) sub_signal = subtract_signal
    call mpi_bcast(sub_signal,    1, MPI_LOGICAL, root, comm_chain, ierr)
  
    residual = cmbmap
    if (sub_signal) residual = residual - signal_map

    ! Distribute foreground amplitudes
    if (num_fg_temp > 0) then
       allocate(my_coeff(num_fg_temp, numband))
       if (myid_chain == root) my_coeff = s%temp_amp
       call mpi_bcast(my_coeff, size(my_coeff), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       do j = 1, num_fg_temp
          if (fix_temp(j,map_id)) residual = residual - my_coeff(j,map_id) * fg_temp(:,:,j)
       end do
       deallocate(my_coeff)
    end if

    if (exclude_pos_fg_amp) then
       allocate(a(0:npix-1,nmaps))
       do i = 1, num_fg_comp
          if (fg_components(i)%enforce_positive_amplitude) then
             if (present(s)) a = s%fg_amp(:,:,i)
             call mpi_bcast(a, size(a), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
             residual = residual - fg_pix_spec_response(:,:,i) * a(pixels,:)
          end if
       end do
       deallocate(a)
    end if

!!$    call int2string(map_id, band)
!!$    call write_map('res'//band//'.fits', residual)
!!$    call mpi_finalize(ierr)
!!$    stop
    
  end subroutine compute_residual



  subroutine compute_signal_map(s)
    implicit none

    type(genvec), optional :: s

    integer(i4b) :: i, j, l, m, ierr
    real(dp),     allocatable, dimension(:,:) :: alms_work
    real(dp),     allocatable, dimension(:,:) :: my_coeff
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: status

    ! Compute the CMB part
    if (.not. enforce_zero_cl) then
       if (myid_alms == root) then
          allocate(alms_work(numcomp, nmaps))
          if (myid_data == root) alms_work = s%cmb_amp
          call mpi_bcast(alms_work, size(alms_work), MPI_DOUBLE_PRECISION, root, comm_data, ierr)
          call multiply_with_beam(alms_work)
          call convert_harmonic_to_real_space(signal_map, alms_work)
          deallocate(alms_work)
          signal_map = signal_map * spec2data(map_id, cmb=.true.)
       else
          call convert_harmonic_to_real_space(signal_map)
       end if
    else
       signal_map = 0.d0
    end if

    ! Add templates
    if (num_fg_temp > 0) then
       allocate(my_coeff(num_fg_temp, numband))
       if (myid_chain == root) my_coeff = s%temp_amp
       call mpi_bcast(my_coeff, size(my_coeff), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       do j = 1, num_fg_temp
          signal_map = signal_map + my_coeff(j,map_id) * fg_temp(:,:,j)
       end do
       deallocate(my_coeff)
    end if

  end subroutine compute_signal_map


  subroutine output_component_map(map_id, s_in, chain_in, iter_in, chain_dir)
    implicit none

    integer(i4b)               :: map_id
    integer(i4b),     optional :: chain_in, iter_in
    type(genvec),     optional :: s_in
    character(len=*), optional :: chain_dir

    integer(i4b) :: i, j, k, l, m, ierr, c, iter, chain
    character(len=512) :: filename, dir
    character(len=2)   :: band, comp
    character(len=5)   :: iter_text
    real(dp),     allocatable, dimension(:,:) :: alms_work
    real(dp),     allocatable, dimension(:,:) :: my_coeff
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: status
    real(dp),     allocatable, dimension(:,:) :: signal
    type(genvec) :: s
    real(dp), allocatable, dimension(:,:) :: outmap, outmap2

    if (.not. any(fg_components%output_freq_comp_maps) .and. .not. output_cmb_freq_map) return
    call allocate_genvec(s)
    if (myid_chain == root) then 
       call genvec_set_equal(s_in, s)
       chain = chain_in
       dir = chain_dir
       iter = iter_in
    end if
    call mpi_bcast(chain, 1, MPI_INTEGER, root, comm_chain, ierr)    
    if (chain /= 1) then
       call deallocate_genvec(s)
       return
    end if

    call bcast_genvec(comm_chain, s)
    call mpi_bcast(dir,  len(dir), MPI_CHARACTER, root, comm_chain, ierr)    
    call mpi_bcast(iter, 1, MPI_INTEGER, root, comm_chain, ierr)    

    if (output_cmb_freq_map) then
       allocate(signal(0:map_size-1,nmaps))
       signal = 0.d0

       if (num_fg_temp > 0) then
          do j = 1, num_fg_temp
             signal = signal + s%temp_amp(j,map_id) * fg_temp(:,:,j)
          end do
       end if

       do k = 1, num_fg_signal
          if (trim(fg_components(k)%type) == 'cmb') cycle
          signal = signal + s%fg_amp(pixels,:,k) * fg_pix_spec_response(:,:,k)
       end do

       ! Compute residual
       signal = cmbmap - signal

       allocate(outmap(0:npix-1,nmaps), outmap2(0:npix-1,nmaps))
       outmap = 0.d0
       outmap(pixels,:) = signal
       call mpi_reduce(outmap, outmap2, size(outmap), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
       if (myid_alms == root) then
          call int2string(map_id, band)
          call int2string(iter, iter_text)
          filename = trim(dir) // '/' // 'cmb_freq_band' // band // '_k' // iter_text // '.fits'
          outmap2 = outmap2 / bp(map_id)%gain / ant2data(map_id) * bp(map_id)%a2t
          call write_map(filename, outmap2, unit='uK_cmb', comptype='CMB')             
       end if
       deallocate(outmap, outmap2)
       deallocate(signal)

    end if

    do c = 1, num_fg_comp

       if (fg_components(c)%output_freq_comp_maps) then
          allocate(outmap(0:npix-1,nmaps), outmap2(0:npix-1,nmaps))
          outmap           = 0.d0
          outmap(pixels,:) = s%fg_amp(pixels,:,c) * fg_pix_spec_response(:,:,c) !* ant2data(map_id)
          call mpi_reduce(outmap, outmap2, size(outmap), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
          if (myid_alms == root) then
             call int2string(map_id, band)
             call int2string(c, comp)
             call int2string(iter, iter_text)
             filename = trim(dir) // '/' // trim(fg_components(c)%label) // '_band' // band // '_k' // &
                  & iter_text // '.fits'
             call write_map(filename, outmap2, unit=fg_components(c)%amp_unit, &
                  & comptype=fg_components(c)%label, nu_ref=fg_components(c)%nu_ref)
          end if
          deallocate(outmap, outmap2)
       end if
       
    end do

    call deallocate_genvec(s)

  end subroutine output_component_map

  subroutine output_mixing_matrix(map_id, chain_in, iter_in, chain_dir)
    implicit none

    integer(i4b)               :: map_id
    integer(i4b),     optional :: chain_in, iter_in
    character(len=*), optional :: chain_dir

    integer(i4b) :: i, j, k, l, m, ierr, c, iter, chain
    character(len=512) :: filename, dir
    character(len=2)   :: band, comp, chain_text
    character(len=5)   :: iter_text
    real(dp),     allocatable, dimension(:,:) :: alms_work
    real(dp),     allocatable, dimension(:,:) :: my_coeff
    integer(i4b), dimension(MPI_STATUS_SIZE)  :: status
    real(dp),     allocatable, dimension(:,:) :: signal
    real(dp), allocatable, dimension(:,:) :: outmap, outmap2

    if (.not. output_mixmat) return
    if (myid_chain == root) then 
       chain = chain_in
       dir = chain_dir
       iter = iter_in
    end if
    call mpi_bcast(chain, 1, MPI_INTEGER, root, comm_chain, ierr)    
    !if (chain /= 1) return
    call mpi_bcast(dir,  len(dir), MPI_CHARACTER, root, comm_chain, ierr)    
    call mpi_bcast(iter, 1, MPI_INTEGER, root, comm_chain, ierr)    

    allocate(outmap(0:npix-1,nmaps), outmap2(0:npix-1,nmaps))
    do i = 1, num_fg_comp
       outmap = 0.d0
       outmap(pixels,:) = fg_pix_spec_response(:,:,i)
       call mpi_reduce(outmap, outmap2, size(outmap), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
       if (myid_alms == root) then
          call int2string(chain, chain_text)
          call int2string(map_id, band)
          call int2string(iter, iter_text)
          call int2string(i, comp)
          filename = trim(dir) // '/' // 'mixmat_comp' // comp // '_band' // band // '_chain' // &
               & chain_text // '_k' // iter_text // '.fits'
          call write_map(filename, outmap2, comptype=fg_components(i)%label)
       end if
    end do
    deallocate(outmap, outmap2)

  end subroutine output_mixing_matrix


  subroutine degrade_map(lmax, map_in, map_out)
    implicit none

    integer(i4b),               intent(in)  :: lmax
    real(dp), dimension(0:,1:), intent(in)  :: map_in
    real(dp), dimension(0:,1:), intent(out) :: map_out

    integer(i4b) :: i, l, m, nside_in, nside_out, npix_in, npix_out, nmaps, numcomp
    real(dp),                        dimension(2)     :: zbounds = 0.d0
    real(dp),           pointer,     dimension(:,:)   :: weights
    complex(dpc),       allocatable, dimension(:,:,:) :: alms_cmplx
    real(dp),           allocatable, dimension(:,:)   :: alms

    npix_in   = size(map_in(:,1))
    npix_out  = size(map_out(:,1))
    nmaps     = size(map_in(0,:))
    nside_in  = nint(sqrt(real(npix_in,sp)/12.))
    nside_out = nint(sqrt(real(npix_out,sp)/12.))
    numcomp   = lmax2numcomp(lmax)

    allocate(alms_cmplx(nmaps, 0:lmax, 0:lmax))
    allocate(alms(numcomp, nmaps))
    allocate(weights(2*nside_in+1,nmaps))

    call read_ringweights(nside_in, weights)
    if (nmaps == 1) then
       call map2alm(nside_in, lmax, lmax, map_in(:,1), alms_cmplx, zbounds, weights)
    else
       call map2alm(nside_in, lmax, lmax, map_in, alms_cmplx, zbounds, weights)
    end if

    call convert_complex_to_real_alms(alms_cmplx, alms)

    call multiply_with_beam(alms, deconvolve=.true.)
    call multiply_with_beam_lowres(alms)
    call convert_real_to_complex_alms(alms, alms_cmplx)

    if (nmaps == 1) then
       call alm2map(nside_out, lmax, lmax, alms_cmplx, map_out(:,1))
    else
       call alm2map(nside_out, lmax, lmax, alms_cmplx, map_out)
    end if

    deallocate(alms, alms_cmplx, weights)

  end subroutine degrade_map


  subroutine output_signal_correlations(s, chisq_map, chain, iter, chain_dir)
    implicit none

    type(genvec),     intent(in) :: s
    integer(i4b),     intent(in) :: chain, iter
    character(len=*), intent(in) :: chain_dir
    real(dp), allocatable, dimension(:,:) :: chisq_map

    integer(i4b) :: i, j, k, l, m, q
    real(dp)     :: zbounds(2), n, corr, mu1, mu2, sigma1, sigma2
    complex(dpc) :: corr_cmplx, auto1, auto2
    real(dp), allocatable, dimension(:) :: mu, sigma
    real(dp), allocatable, dimension(:,:) :: weights, mask, corr_map
    real(dp), allocatable, dimension(:,:,:) :: maps
    complex(dpc), allocatable, dimension(:,:,:,:) :: alms
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmb
    character(len=512) :: filename
    character(len=2) :: chain_text
    character(len=5) :: iter_text
    character(len=2) :: c1_text, c2_text

    allocate(alms(nmaps,0:lmax_corr,0:lmax_corr, 0:num_fg_comp), weights(1:2*nside, nmaps), &
         & mask(0:npix-1,1), maps(0:npix-1,1,0:num_fg_comp))
    weights = 1.d0
    zbounds = 0.d0
    mask    = mask_corr

    maps(:,:,1:num_fg_comp) = s%fg_amp

    if (.not. enforce_zero_cl) then
       allocate(alms_cmb(nmaps,0:lmax,0:lmax))
       call convert_real_to_complex_alms(s%cmb_amp, alms_cmb)
       alms_cmb(:,0:1,:) = cmplx(0.d0,0.d0)
       do l = 2, lmax
          alms_cmb(1,l,:) = alms_cmb(1,l,:) * beam(l,1)
       end do
       if (nmaps == 1) then
          call alm2map(nside, lmax, lmax, alms_cmb, maps(:,1,0))
       else
          call alm2map(nside, lmax, lmax, alms_cmb, maps(:,:,0))
       end if
       deallocate(alms_cmb)
    end if

    ! Compute mean and standard deviation
    allocate(mu(0:num_fg_comp), sigma(0:num_fg_comp))
    mu = 0.d0
    n  = sum(mask(:,1))
    do j = 0, num_fg_comp
       mu(j) = sum(mask(:,1) * maps(:,1,j)) / n
    end do

    sigma = 0.d0
    do j = 0, num_fg_comp
       sigma(j) = sqrt(sum(mask(:,1)*(maps(:,1,j)-mu(j))**2) / (n-1))
    end do    

    ! Compute alms
    alms = 0.d0
    do i = 0, num_fg_comp
       if (nmaps == 1) then
          call map2alm(nside, lmax_corr, lmax_corr, (maps(:,1,i)-mu(i))*mask(:,1), &
               & alms(:,:,:,i), zbounds, weights)
       else
          call map2alm(nside, lmax_corr, lmax_corr, (maps(:,:,i)-mu(i))*mask, &
               & alms(:,:,:,i), zbounds, weights)
       end if
    end do

    call convert_ring2nest(nside, mask(:,1))
    do i = 0, num_fg_comp
       call convert_ring2nest(nside, maps(:,1,i))
    end do

    ! Output real-space correlation coefficients
    call int2string(chain, chain_text)
    call int2string(iter, iter_text)
    filename = trim(chain_dir) // '/fg_pix_corr_c' // chain_text // '_k' // iter_text // '.dat'
    open(58,file=trim(filename), recl=1024)
    do i = 0, num_fg_comp
       do j = 0, i
          corr = sum(mask(:,1) * (maps(:,1,i)-mu(i)) * (maps(:,1,j)-mu(j))) / (n-1)
          corr = corr / (sigma(i)*sigma(j))
          write(58,fmt='(f6.3,a)',advance="no") corr, ' '
       end do
       write(58,*) 
    end do
    close(58)

    ! Output cross spectra
    do i = 0, num_fg_comp
       call int2string(i,c1_text)
       do j = i+1, num_fg_comp
          call int2string(j,c2_text)
          filename = trim(chain_dir) // '/fg_cross_corr_' // c1_text // '_' // c2_text // '_c' // &
               & chain_text // '_k' // iter_text // '.dat'
          open(58,file=trim(filename))
          do l = 2, lmax_corr
             corr_cmplx = alms(1,l,0,i)*alms(1,l,0,j)
             auto1      = alms(1,l,0,i)*alms(1,l,0,i)
             auto2      = alms(1,l,0,j)*alms(1,l,0,j)
             do m = 1, l
                corr_cmplx = corr_cmplx + alms(1,l,m,i) * conjg(alms(1,l,m,j)) + &
                     & alms(1,l,m,j) * conjg(alms(1,l,m,i))
                auto1      = auto1      + 2.d0 * alms(1,l,m,i) * conjg(alms(1,l,m,i))
                auto2      = auto2      + 2.d0 * alms(1,l,m,j) * conjg(alms(1,l,m,j))
             end do
             write(58,*) l, real(corr_cmplx,sp) / sqrt(real(auto1,dp)*real(auto2,dp))
          end do
          close(58)
       end do
    end do

    ! Output cross-maps
    allocate(corr_map(0:12*nside_corr**2-1,1))
    q = (nside / nside_corr)**2
    do i = 0, num_fg_comp
       call int2string(i,c1_text)
       do j = i+1, num_fg_comp
          call int2string(j,c2_text)
          filename = trim(chain_dir) // '/fg_cross_map_' // c1_text // '_' // c2_text // '_c' // &
               & chain_text // '_k' // iter_text // '.fits'

          corr_map = -1.6375d30
          do k = 0, 12*nside_corr**2-1
             n = sum(mask(k*q:(k+1)*q-1,1))
             if (n < 3) cycle
             mu1    = sum(mask(k*q:(k+1)*q-1,1) * maps(k*q:(k+1)*q-1,1,i)) / n
             mu2    = sum(mask(k*q:(k+1)*q-1,1) * maps(k*q:(k+1)*q-1,1,j)) / n
             sigma1 = sqrt(sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,i)-mu1)**2) / (n-1))
             sigma2 = sqrt(sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,j)-mu2)**2) / (n-1))
             corr_map(k,1) = sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,i)-mu1) * &
                  & (maps(k*q:(k+1)*q-1,1,j)-mu2)) / (n-1) / (sigma1*sigma2)
          end do
          call convert_nest2ring(nside_corr, corr_map(:,1))
          call write_map(filename, corr_map, comptype='Correlation coefficient')             
       end do
    end do

    deallocate(alms, mu, sigma, weights, mask, corr_map, maps)

  end subroutine output_signal_correlations


  subroutine enforce_zero_cmb_md(s, mask)
    implicit none

    type(genvec),                intent(inout) :: s
    real(dp),     dimension(0:), intent(in)    :: mask

    integer(i4b) :: i, j, k, ref_comp
    real(dp)     :: A(4,4), b(4)
    real(dp),     allocatable, dimension(:,:)   :: map, alms_real, weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms

    ! Compute mask kernel
    A = 0.d0
    do i = 0, npix-1
       if (mask(i) == 1.d0) then
          do j = 1, 4
             do k = 1, 4
                A(j,k) = A(j,k) + fg_temp_lowres(i,1,1,j)*fg_temp_lowres(i,1,1,k)
             end do
          end do
       end if
    end do
    call invert_matrix(A)

    ! Compute total map
    allocate(map(0:npix-1,nmaps))
    map = 0.d0
    if (enforce_zero_cl) then
       do i = 1, num_fg_comp
          if (trim(fg_components(i)%type) == 'cmb') then
             ref_comp = i
             map      = s%fg_amp(:,:,i)
          end if
       end do
    else
       allocate(alms(nmaps,0:lmax,0:lmax))
       call convert_real_to_complex_alms(s%cmb_amp, alms)
       call alm2map(nside, lmax, lmax, alms(1:1,0:lmax,0:lmax), map(:,1))       
    end if

    ! Compute CMB monopole and dipole coefficients
    b = 0.d0
    do i = 0, npix-1
       if (mask(i) == 1.d0) b = b + fg_temp_lowres(i,1,1,1:4) * map(i,1)
    end do
    b = matmul(A, b)
       
    ! Subtract from CMB solution
    do i = 1, 4
       map(:,1) = map(:,1) - b(i) * fg_temp_lowres(:,1,1,i)
    end do

    if (enforce_zero_cl) then
       do i = 1, num_fg_comp
          if (trim(fg_components(i)%type) == 'cmb') s%fg_amp(:,:,i) = map
       end do
    else
       map = 0.d0
       do i = 1, 4
          map(:,1) = map(:,1) + b(i) * fg_temp_lowres(:,1,1,i)
       end do
       deallocate(alms)
       allocate(weights(2*nside+1,1), alms(1,0:lmax,0:lmax), alms_real(numcomp,1))
       alms = cmplx(0.d0, 0.d0)
       weights = 1.d0
       call map2alm(nside, lmax, lmax, map(:,1), alms, [0.d0,0.d0], weights)
       call convert_complex_to_real_alms(alms, alms_real)
       s%cmb_amp(1:4,1) = s%cmb_amp(1:4,1) - alms_real(1:4,1)
       deallocate(weights, alms, alms_real)

    end if

    ! Update template coefficients
    do i = 1, numband
       if (enforce_zero_cl) then
          ! CMB map is in brightness temperature at a (sharp) reference frequency
          s%temp_amp(1:4,i) = s%temp_amp(1:4,i) + b * compute_ant2thermo_single(fg_components(ref_comp)%nu_ref) &
               & * spec2data(i, cmb=.true.)/bp(i)%gain
       else
          ! CMB map is in thermodynamic temperature
          s%temp_amp(1:4,i) = s%temp_amp(1:4,i) + b * spec2data(i, cmb=.true.)/bp(i)%gain
       end if
    end do

    deallocate(map)

  end subroutine enforce_zero_cmb_md

  ! Solve (F^t invN F)^1 F^t invN d, and compute (F^t invN F)^1 
  subroutine output_ml_map_engine(myid_chain, paramfile)    
    implicit none

    integer(i4b),     intent(in)           :: myid_chain
    character(len=*), intent(in), optional :: paramfile

    integer(i4b) :: i, j, k, l, n, i2, j2, k2, p1, p2, ierr, unit, ntot
    integer(i4b), allocatable, dimension(:) :: numval
    character(len=512) :: chain_dir
    real(dp)     :: scale
    real(dp), allocatable, dimension(:)   :: Ft_invN_d, buffer1, x
    real(dp), allocatable, dimension(:,:) :: F, Ft_invN_F, buffer2, map, cov, rms

    if (myid_chain == root) call get_parameter(paramfile, 'CHAIN_DIRECTORY',    par_string=chain_dir)

    allocate(numval(nmaps))
    do i = 1, nmaps
       numval(i) = count(mask(:,i))
    end do
    ntot   = sum(numval)
    n      = ntot*num_fg_comp
    unit   = getlun()

    if (myid_chain == root) write(*,*) 'Size of ML system = ', n

    allocate(F(ntot,numcomp*ntot), Ft_invN_d(n), Ft_invN_F(n,n), x(npix*nmaps), map(0:npix-1,nmaps))
    allocate(rms(0:npix-1,nmaps))
    
    ! Set up linear system
    if (myid_chain == root) write(*,*) 'Setting up Ft*invN*F matrix'
    map = residual
    call multiply_by_inv_N(map)
    i = 1
    do j = 1, num_fg_comp
       do k = 1, nmaps
          do p1 = 0, npix-1
             if (mask(p1,k)) then
                Ft_invN_d(i) = fg_pix_spec_response(p1,k,j) * map(p1,k)
                i            = i+1
             end if
          end do
       end do
    end do

    i = 1
    do j = 1, num_fg_comp
       if (myid_chain == root) write(*,*) 'Component =', j, num_fg_comp
       do k = 1, nmaps
          do p1 = 0, npix-1
             if (mask(p1,k)) then

                map = 0.d0
                map(p1,k) = fg_pix_spec_response(p1,k,j)
                call multiply_by_inv_N(map)
                
                i2 = 1
                do j2 = 1, num_fg_comp
                   do k2 = 1, nmaps
                      do p2 = 0, npix-1
                         if (mask(p2,k2)) then
                            Ft_invN_F(i,i2) = map(p2,k2) * fg_pix_spec_response(p2,k2,j2)
                            i2              = i2+1
                         end if
                      end do
                   end do
                end do
                
                i = i+1                
             end if
          end do
       end do
    end do
    
    ! Reduce system
    if (myid_chain == root) write(*,*) 'Reducing linear system across frequencies'
    allocate(buffer1(n))
    call mpi_reduce(Ft_invN_d, buffer1, size(buffer1), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
    if (myid_chain == root) Ft_invN_d = buffer1
    deallocate(buffer1)

    allocate(buffer2(n,n))
    call mpi_reduce(Ft_invN_F, buffer2, size(buffer2), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
    if (myid_chain == root) Ft_invN_F = buffer2
    deallocate(buffer2)
    
    ! Output CMB map and covariance
    if (myid_chain == root) then

       if (myid_chain == root) write(*,*) 'Inverting Ft*invN*F'
       call invert_matrix(Ft_invN_F, cholesky=.true.)
       Ft_invN_d = matmul(Ft_invN_F, Ft_invN_d)
       do i = 1, num_fg_comp
          if (myid_chain == root) write(*,*) 'Outputting map and covariance matrix'

          scale = ant2unit(fg_components(i)%amp_unit, fg_components(i)%nu_ref, &
               & band_ref=fg_components(i)%ref_band)
          
          ! Output map
          call conv_mask2map(Ft_invN_d((i-1)*ntot+1:i*ntot)*scale, map)
          call write_map(trim(chain_dir)//'/ML_'//trim(fg_components(i)%label)//'_map.fits', map, &
               & unit=fg_components(i)%amp_unit, comptype=fg_components(i)%label, &
               & nu_ref=fg_components(i)%nu_ref)
          
          ! Output covariance matrix
          open(unit,file=trim(chain_dir)//'/ML_'//trim(fg_components(i)%label)//'_N.unf',form='unformatted')
          write(unit) npix*nmaps
          write(unit) 1 ! Ordering
          write(unit) nmaps ! Polarization status
          p1 = (i-1)*ntot+1
          do j = 1, nmaps
             do k = 0, npix-1
                if (mask(k,j)) then
                   call conv_mask2map(Ft_invN_F(p1,(i-1)*ntot+1:i*ntot), map)
                   do l = 1, nmaps
                      x((l-1)*npix+1:l*npix) = map(:,l) 
                   end do
                   where (x /= -1.6375d30)
                      x = x*scale**2
                   end where
                   p1 = p1+1
                else
                   x = 0.d0
                end if
                write(unit) x
             end do
          end do
          write(unit) .false.
          close(unit)

          ! Output rms
          rms = -1.6375d30
          l   = (i-1)*ntot+1
          do j = 1, nmaps
             do k = 0, npix-1
                if (mask(k,j)) then
                   rms(k,j) = sqrt(Ft_invN_F(l,l)) * scale
                   l        = l+1
                end if
             end do
          end do
          call write_map(trim(chain_dir)//'/ML_'//trim(fg_components(i)%label)//'_rms.fits', rms, &
               & unit=fg_components(i)%amp_unit, comptype=fg_components(i)%label, nu_ref=fg_components(i)%nu_ref)
       end do
    end if

    deallocate(Ft_invN_d, Ft_invN_F, x, map, rms)


  end subroutine output_ml_map_engine

  subroutine conv_mask2map(linmap, map)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: linmap
    real(dp), dimension(0:,1:), intent(out) :: map

    integer(i4b) :: i, j, k

    map = -1.6375d30
    i   = 1
    do j = 1, nmaps
       do k = 0, npix-1
          if (mask(k,j)) then
             map(k,j) = linmap(i)
             i        = i+1
          end if
       end do
    end do

  end subroutine conv_mask2map



end module comm_data_mod
