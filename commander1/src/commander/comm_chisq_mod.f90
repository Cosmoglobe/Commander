module comm_chisq_mod
  use comm_data_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b),                            private :: numproc_per_band
  integer(i4b),                            private :: myid_data, comm_data, root, ierr
  integer(i4b),                            private :: myid_alms, comm_alms
  integer(i4b),                            private :: myid_chain, comm_chain
  integer(i4b),                            private :: num_bp_step
  logical(lgt),                            private :: sample_inside_mask, output_band_chisq, enforce_zero_cl
  logical(lgt),                            private :: output_residuals
  character(len=512),                      private :: chain_dir, operation
  integer(i4b), allocatable, dimension(:,:), private :: g_l

  real(dp),     allocatable, dimension(:,:)       :: g_gauss
  real(dp),     allocatable, dimension(:,:)       :: procmask, procmask_full

contains

  subroutine initialize_chisq_mod(paramfile, comm_chain_in, comm_alms_in, comm_data_in)
    implicit none

    character(len=128),                intent(in) :: paramfile
    integer(i4b),                      intent(in) :: comm_chain_in, comm_alms_in, comm_data_in

    integer(i4b) :: i
    logical(lgt)       :: exist, apply_gain_corr
    character(len=2)   :: i_text
    character(len=256) :: paramtext, procmaskfile
    real(dp), allocatable, dimension(:,:) :: mask_in

    root        = 0
    comm_data = comm_data_in
    call mpi_comm_rank(comm_data, myid_data, ierr)
    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    call mpi_comm_size(comm_alms, numproc_per_band, ierr)
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)

    call get_parameter(paramfile, 'SAMPLE_INSIDE_MASK',     par_lgt=sample_inside_mask)
    call get_parameter(paramfile, 'OUTPUT_BAND_CHISQ',      par_lgt=output_band_chisq)
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL',        par_lgt=enforce_zero_cl)
    call get_parameter(paramfile, 'CHAIN_DIRECTORY',        par_string=chain_dir)
    call get_parameter(paramfile, 'OPERATION',              par_string=operation)
    call get_parameter(paramfile, 'APPLY_GAIN_CORRECTIONS', par_lgt=apply_gain_corr)
    call get_parameter(paramfile, 'NUM_BP_MCMC_SUBSTEPS',   par_int=num_bp_step)
    call get_parameter(paramfile, 'PROCESSING_MASK',        par_string=procmaskfile)

    allocate(procmask(0:map_size-1,nmaps), procmask_full(0:npix-1,nmaps))
    call read_map(procmaskfile, procmask_full)
    where (procmask_full > 0.5)
       procmask_full = 1.d0
    elsewhere
       procmask_full = 0.d0
    end where
    procmask = procmask_full(pixels,:)
    !deallocate(mask_in)

    allocate(g_gauss(numband,2), g_l(numband,2))
    if (apply_gain_corr) then
       call get_parameter(paramfile, 'OUTPUT_GAIN_RESIDUALS', par_lgt=output_residuals)
       do i = 1, numband
          call int2string(i, i_text)
          paramtext = 'GAIN_PRIOR_MEAN' // i_text
          call get_parameter(paramfile, trim(paramtext), par_dp=g_gauss(i,1))
          paramtext = 'GAIN_PRIOR_RMS' // i_text
          call get_parameter(paramfile, trim(paramtext), par_dp=g_gauss(i,2))
          paramtext = 'GAIN_LMIN' // i_text
          call get_parameter(paramfile, trim(paramtext), par_int=g_l(i,1))
          paramtext = 'GAIN_LMAX' // i_text
          call get_parameter(paramfile, trim(paramtext), par_int=g_l(i,2))
       end do
    else
       g_gauss(:,1) = 1.d0
       g_gauss(:,2) = 0.d0
       g_l          = -1
    end if

  end subroutine initialize_chisq_mod

  subroutine compute_total_chisq(map_id, s, output_stats, index_map, chisq_map, chisq_highlat, chisq_fullsky, &
       & chisq_rms, chisq_band, nu_band, chisq_band_fullsky, chain, iter)
    implicit none

    integer(i4b),                      intent(in)              :: map_id
    integer(i4b),                      intent(in),    optional :: chain
    logical(lgt),                      intent(in),    optional :: output_stats
    type(genvec),                      intent(in),    optional :: s
    real(dp),     dimension(0:,1:,1:), intent(in),    optional :: index_map
    real(dp),                          intent(out),   optional :: chisq_fullsky, chisq_highlat, chisq_rms
    real(dp),     dimension(0:,1:),    intent(out),   optional :: chisq_map
    real(dp),     dimension(1:),       intent(out),   optional :: chisq_band, nu_band, chisq_band_fullsky
    integer(i4b),                      intent(in),    optional :: iter

    integer(i4b) :: i, j, k, l, m, p, ierr, c, testpix(4)
    logical(lgt) :: output_chisq_rms, use_index_map, output_chisq_fullsky, output_chisq_band
    logical(lgt) :: output_chisq_highlat, output_stats_, output_chisq, output_chisq_band_fullsky
    real(dp)     :: my_chisq, chisq_tot, f, r, monopole
    character(len=512) :: dir, filename
    character(len=2)   :: band, chain_text
    character(len=5)   :: iter_text
    real(dp),     allocatable, dimension(:)     :: my_chisq_band, buffer
    real(dp),     allocatable, dimension(:,:,:) :: fg_amp, ind_map
    real(dp),     allocatable, dimension(:,:)   :: alms_work, my_signal, my_coeff, mask_output
    real(dp),     allocatable, dimension(:,:)   :: my_chisq_map, chisq_map_tot, outmap
    real(dp),     allocatable, dimension(:,:)   :: residual, invN_residual1, invN_residual2
    integer(i4b), dimension(MPI_STATUS_SIZE)    :: status
    type(fg_params) :: fg_par

    if (myid_chain == root) then
       output_stats_ = .false.; if (present(output_stats)) output_stats_ = output_stats
       output_chisq_rms           = present(chisq_rms)
       output_chisq_fullsky       = present(chisq_fullsky)
       output_chisq_highlat       = present(chisq_highlat)
       output_chisq_band          = present(chisq_band) 
       output_chisq_band_fullsky  = present(chisq_band_fullsky)
       use_index_map              = present(index_map)
    end if
    call mpi_bcast(output_chisq_rms,          1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(output_chisq_fullsky,      1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(output_chisq_highlat,      1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(output_stats_,             1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(output_chisq_band,         1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(output_chisq_band_fullsky, 1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(use_index_map,             1, MPI_LOGICAL, root, comm_chain, ierr)


    allocate(my_signal(0:map_size-1,nmaps))
    allocate(my_chisq_map(0:npix-1,nmaps))
    allocate(chisq_map_tot(0:npix-1,nmaps))
    allocate(residual(0:map_size-1,nmaps))
    allocate(invN_residual1(0:map_size-1,nmaps), invN_residual2(0:map_size-1,nmaps))
    if (use_index_map) then
       allocate(ind_map(0:npix-1,nmaps,num_fg_par))
       if (myid_chain == root) ind_map = index_map
       call mpi_bcast(ind_map, size(ind_map), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
    end if

    ! Compute the CMB part
    my_signal = 0.d0
    if (.not. enforce_zero_cl) then
       if (myid_alms == root) then
          allocate(alms_work(numcomp, nmaps))
          if (myid_data == root) alms_work = s%cmb_amp
          call mpi_bcast(alms_work, size(alms_work), MPI_DOUBLE_PRECISION, root, comm_data, ierr)
          call multiply_with_beam(alms_work)
          call convert_harmonic_to_real_space(my_signal, alms_work)
          deallocate(alms_work)
       else
          call convert_harmonic_to_real_space(my_signal)
       end if
       my_signal = my_signal * spec2data(map_id, cmb=.true.)
       where (.not. mask)
          my_signal = 0.d0
       end where
    end if

    ! Distribute foreground amplitudes
    if (num_fg_temp > 0) then
       allocate(my_coeff(num_fg_temp, numband))
       if (myid_chain == root) my_coeff = s%temp_amp
       call mpi_bcast(my_coeff, size(my_coeff), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       !write(*,fmt='(i5,5f16.3)') map_id, my_coeff(:,map_id)
       do j = 1, num_fg_temp
          my_signal = my_signal + my_coeff(j,map_id) * fg_temp(:,:,j)
       end do
       monopole = 0.d0; if (sample_T_modes) monopole = my_coeff(1,map_id)
       deallocate(my_coeff)
    end if

    if (sample_fg_pix) then
       allocate(fg_amp(0:npix-1,nmaps, num_fg_signal))
       if (myid_chain == root) fg_amp = s%fg_amp
       call mpi_bcast(fg_amp, size(fg_amp), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       !write(*,*) 'LEAVING FREE-FREE!!'
       do k = 1, num_fg_signal
          do j = 1, nmaps
             do i = 0, map_size-1
                if (.not. enforce_zero_cl .and. trim(fg_components(k)%type) == 'cmb' .and. mask(i,j)) cycle
                !if (trim(fg_components(k)%type) == 'freefree') cycle
                if (use_index_map) then
                   call reorder_fg_params(ind_map(pixels(i),:,:), fg_par)
                   f = get_effective_fg_spectrum(fg_components(k), map_id, fg_par%comp(k)%p(j,:), &
                        & pixel=pixels(i), pol=j)
                else
                   f = fg_pix_spec_response(i,j,k)
                end if
                my_signal(i,j) = my_signal(i,j) + f * fg_amp(pixels(i),j,k)
             end do
          end do
       end do
       if (allocated(fg_par%comp)) call deallocate_fg_params(fg_par)
       deallocate(fg_amp)
    end if

    ! Compute the chi-square
    my_chisq_map = 0.d0
    do i = 0, map_size-1
       residual(i,:) = cmbmap(i,:) - my_signal(i,:)
    end do
    residual = residual

!!$    call int2string(map_id, band)
!!$    call write_map('chires'//band//'.fits', residual)
!!$    call mpi_finalize(ierr)
!!$    stop

    call multiply_by_sqrt_inv_N(residual, invN_residual1)
    if (sample_inside_mask) then
       mask_state = INSIDE_MASK
       call multiply_by_sqrt_inv_N(residual, invN_residual2)
    else
       invN_residual2 = 0.d0
    end if
    mask_state = OUTSIDE_MASK

    do j = 1, nmaps
       do i = 0, map_size-1
!       my_chisq_map(pixels(i),1) = sum(residual(i,:)*invN_residual1(i,:)) + &
!            & sum(residual(i,:)*invN_residual2(i,:))
          my_chisq_map(pixels(i),j) = invN_residual1(i,j)**2 + invN_residual2(i,j)**2
       end do
    end do

    ! First output per-band chi-squares
    if (myid_chain == root) output_chisq = present(iter)
    call mpi_bcast(output_chisq, 1, MPI_LOGICAL, root, comm_chain, ierr)    

    if (output_chisq_band) then
       allocate(my_chisq_band(numband), buffer(numband))
       ! Return chisq per band
       my_chisq_band = 0.d0
       my_chisq_band(map_id) = sum(invN_residual1**2 * procmask)
       call mpi_reduce(my_chisq_band, buffer, numband, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
       if (myid_chain == root) chisq_band = buffer
       ! Return number of degrees of freedom per band
       my_chisq_band = 0.d0
       my_chisq_band(map_id) = count(invN_residual1 /= 0.d0 * procmask)
       call mpi_reduce(my_chisq_band, buffer, numband, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
       if (myid_chain == root) nu_band = buffer
       deallocate(my_chisq_band, buffer)
    end if

    if (output_chisq_band_fullsky) then
       allocate(my_chisq_band(numband), buffer(numband))
       ! Return chisq per band
       my_chisq_band = 0.d0
       my_chisq_band(map_id) = sum(invN_residual1**2 * procmask) + sum(invN_residual2**2 * procmask)
       call mpi_reduce(my_chisq_band, buffer, numband, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
       if (myid_chain == root) chisq_band_fullsky = buffer
       deallocate(my_chisq_band, buffer)
    end if

    if (output_stats_ .and. output_band_chisq .and. output_chisq) then

       if (myid_chain == 0) then
          c   = chain
          i   = iter
       end if
       call mpi_bcast(i,   1,        MPI_INTEGER,   root, comm_chain, ierr)    
       call mpi_bcast(c,   1,        MPI_INTEGER,   root, comm_chain, ierr)    
       
       call int2string(map_id, band)
       call int2string(i,      iter_text)
       call int2string(c,      chain_text)
       filename = trim(chain_dir) // '/chisq_' // trim(bp(map_id)%label) // '_c' // chain_text // '_k' // iter_text // '.fits'
       call mpi_reduce(my_chisq_map, chisq_map_tot, size(my_chisq_map), MPI_DOUBLE_PRECISION, &
            & MPI_SUM, root, comm_alms, ierr)
       
       allocate(mask_output(0:npix-1,nmaps))
       mask_output = 1.d0
       if (.not. sample_inside_mask) then
          where (my_chisq_map == 0.d0) 
             mask_output = 0.d0
          end where
       end if
       where (mask_output == 0.d0) 
          chisq_map_tot = -1.6375d30
       end where
       if (myid_alms == root) call write_map(filename, chisq_map_tot, comptype='Chisq', ttype='Chisq')

       allocate(outmap(0:npix-1,nmaps))
       outmap = 0.d0
       outmap(pixels,:) = residual
       filename = trim(chain_dir) // '/residual_' // trim(bp(map_id)%label) // '_c' // chain_text // '_k' // iter_text // '.fits'
       call mpi_reduce(outmap, chisq_map_tot, size(my_chisq_map), MPI_DOUBLE_PRECISION, &
            & MPI_SUM, root, comm_alms, ierr)
       where (mask_output == 0.d0) 
          chisq_map_tot = -1.6375d30
       end where
       if (myid_alms == root) call write_map(filename, chisq_map_tot, comptype='Residual', &
            & unit=bp(map_id)%unit, nu_ref=bp(map_id)%nu_c)
       deallocate(outmap)

       allocate(outmap(0:npix-1,nmaps))
       outmap = 0.d0
       do i = 1, nmaps
          if (i == 1) then
             !outmap(pixels,i) = residual(:,i) / (cmbmap(:,i)-monopole)
          else
             !outmap(pixels,i) = residual(:,i) / cmbmap(:,i)
          end if
       end do
       filename = trim(chain_dir) // '/frac_residual_' // trim(bp(map_id)%label) // '_c' // chain_text &
            & // '_k' // iter_text // '.fits'
       call mpi_reduce(outmap, chisq_map_tot, size(my_chisq_map), MPI_DOUBLE_PRECISION, &
            & MPI_SUM, root, comm_alms, ierr)
       where (mask_output == 0.d0) 
          chisq_map_tot = -1.6375d30
       end where
       if (myid_alms == root) call write_map(filename, chisq_map_tot, comptype='Chisq', ttype='Chisq')
       deallocate(outmap, mask_output)
       
    end if

    ! Then output the co-added chisquare map
    call mpi_reduce(my_chisq_map, chisq_map_tot, size(my_chisq_map), MPI_DOUBLE_PRECISION, &
         & MPI_SUM, root, comm_chain, ierr)
    if (myid_chain == root .and. present(chisq_map))   chisq_map = chisq_map_tot
    if (myid_chain == root .and. output_chisq_fullsky) chisq_fullsky = sum(chisq_map_tot*procmask_full)

    ! Output high-latitude value, if requested
    if (output_chisq_highlat) then
!       call mpi_reduce(sum(residual*invN_residual1), chisq_tot, 1, MPI_DOUBLE_PRECISION, &
!            & MPI_SUM, root, comm_chain, ierr)
       call mpi_reduce(sum(invN_residual1**2 * procmask), chisq_tot, 1, MPI_DOUBLE_PRECISION, &
            & MPI_SUM, root, comm_chain, ierr)
       if (myid_chain == root) chisq_highlat = chisq_tot
    end if

    ! Compute chisquare defined by RMS if requested
    if (output_chisq_rms) then
       !call multiply_by_inv_N(residual, invN_residual1, N_format='rms')
       call multiply_by_sqrt_inv_N(residual, invN_residual1, N_format='rms')
       !my_chisq = sum(residual*invN_residual1)
       my_chisq = sum(invN_residual1**2 * procmask)
       if (sample_inside_mask) then
          mask_state = INSIDE_MASK
          !call multiply_by_inv_N(residual, invN_residual2, N_format='rms')
          call multiply_by_sqrt_inv_N(residual, invN_residual2, N_format='rms')
          mask_state = OUTSIDE_MASK
          !my_chisq = my_chisq + sum(residual*invN_residual2)
          my_chisq = my_chisq + sum(invN_residual2**2 * procmask)
       end if
       call mpi_reduce(my_chisq, chisq_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)
       if (myid_chain == root) chisq_rms = chisq_tot
    end if

    if (allocated(ind_map)) deallocate(ind_map)
    deallocate(residual)
    deallocate(invN_residual1)
    deallocate(invN_residual2)
    deallocate(my_signal)
    deallocate(my_chisq_map)
    deallocate(chisq_map_tot)

  end subroutine compute_total_chisq


  subroutine sample_gain(handle, map_id, s, chain, iter)
    implicit none

    type(planck_rng),                  intent(inout)           :: handle
    integer(i4b),                      intent(in)              :: map_id
    type(genvec),                      intent(in),    optional :: s
    integer(i4b),                      intent(in),    optional :: iter, chain

    integer(i4b) :: i, j, k, l, m, p, ierr, c, testpix(4), lmin, lmax, mychain, myiter, dl_low, dl_high, ncomp
    integer(i4b) :: numpair, numval
    real(dp)     :: my_chisq, chisq_tot, f, r, mu, sigma, my_mu, my_sigma, gain_new, old_gain
    character(len=512) :: filename
    character(len=2)   :: band, chain_text
    character(len=5)   :: iter_text
    character(len=5)   :: id_text
    real(dp)           :: MAX_DELTA_G = 0.05d0, md(4)
    integer(i4b), allocatable, dimension(:)     :: mask2map
    real(dp),     allocatable, dimension(:)     :: my_gain, buffer, cl1, cl2, slope
    real(dp),     allocatable, dimension(:,:,:) :: fg_amp, ind_map
    real(dp),     allocatable, dimension(:,:)   :: alms_work, my_signal, my_coeff, map, weights
    real(dp),     allocatable, dimension(:,:)   :: my_chisq_map, chisq_map_tot
    real(dp),     allocatable, dimension(:,:)   :: residual, invN_signal
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: alms1, alms2
    integer(i4b), dimension(MPI_STATUS_SIZE)    :: status
    type(fg_params) :: fg_par

    dl_low  = 10
    dl_high = 25

    allocate(my_signal(0:map_size-1,nmaps), invN_signal(0:map_size-1,nmaps))
    allocate(residual(0:map_size-1,nmaps))

    ! Compute residual
    residual = cmbmap
    if (num_fg_temp > 0) then
       allocate(my_coeff(num_fg_temp, numband))
       if (myid_chain == root) my_coeff = s%temp_amp
       call mpi_bcast(my_coeff, size(my_coeff), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       do j = 1, num_fg_temp
          residual = residual - my_coeff(j,map_id) * fg_temp(:,:,j)
       end do
       deallocate(my_coeff)
    end if

    my_signal = 0.d0
    if (sample_fg_pix) then
       allocate(fg_amp(0:npix-1,nmaps, num_fg_signal))
       if (myid_chain == root) fg_amp = s%fg_amp
       call mpi_bcast(fg_amp, size(fg_amp), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       do k = 1, num_fg_signal
          if (k == bp(map_id)%comp_calib) cycle
          do j = 1, nmaps
             do i = 0, map_size-1
                my_signal(i,j) = my_signal(i,j) + fg_pix_spec_response(i,j,k) * fg_amp(pixels(i),j,k)
             end do
          end do
       end do
       if (bp(map_id)%comp_calib > 0) residual = residual - my_signal
       if (allocated(fg_par%comp)) call deallocate_fg_params(fg_par)
       deallocate(fg_amp)
    end if

    ! Compute sky signal
    my_signal = 0.d0
    if (.not. enforce_zero_cl) then
       if (myid_alms == root) then
          allocate(alms_work(numcomp, nmaps))
          if (myid_data == root) alms_work = s%cmb_amp
          call multiply_with_beam(alms_work)
          call convert_harmonic_to_real_space(my_signal, alms_work)
          deallocate(alms_work)
       else
          call convert_harmonic_to_real_space(my_signal)
       end if
       my_signal = my_signal * spec2data(map_id, cmb=.true.)
       where (.not. mask)
          my_signal = 0.d0
       end where
    end if

    if (sample_fg_pix) then
       allocate(fg_amp(0:npix-1,nmaps, num_fg_signal))
       if (myid_chain == root) fg_amp = s%fg_amp
       call mpi_bcast(fg_amp, size(fg_amp), MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
       do k = 1, num_fg_signal
          do j = 1, nmaps
             do i = 0, map_size-1
                if (bp(map_id)%comp_calib > 0 .and. k /= bp(map_id)%comp_calib) cycle
                if (.not. enforce_zero_cl .and. trim(fg_components(k)%type) == 'cmb' .and. mask(i,j)) cycle
                my_signal(i,j) = my_signal(i,j) + fg_pix_spec_response(i,j,k) * fg_amp(pixels(i),j,k)
             end do
          end do
       end do
       if (allocated(fg_par%comp)) call deallocate_fg_params(fg_par)
       deallocate(fg_amp)
    end if

    ! Divide by old gain
    my_signal = my_signal / bp(map_id)%gain

    ! Apply mask
    my_signal = my_signal * mask_calib * procmask
    residual  = residual  * mask_calib * procmask

    ! Remove monopole and dipole if calib_comp is CMB
!!$    if (bp(map_id)%comp_calib > 0) then
!!$       if (trim(fg_components(bp(map_id)%comp_calib)%type) == 'cmb') then
!!$          call remove_dipole(nside, my_signal(:,1), 1, 2, md, [-1.d0,1.d0], mask=mask_calib(:,1))
!!$          call remove_dipole(nside, residual(:,1),  1, 2, md, [-1.d0,1.d0], mask=mask_calib(:,1))
!!$       end if
!!$    end if

!!$    if (myid_chain == 10) then
!!$       call write_map('ka2_res.fits', residual, comptype='Residual', unit=bp(map_id)%unit, nu_ref=bp(map_id)%nu_c)
!!$       call write_map('ka2_sig.fits', my_signal, comptype='Residual', unit=bp(map_id)%unit, nu_ref=bp(map_id)%nu_c)
!!$    end if
!!$    call mpi_finalize(ierr)
!!$    stop


    allocate(my_gain(numband))
    lmin = g_l(map_id,1)
    lmax = g_l(map_id,2)
    my_gain  = 0.d0
    if (lmin > 0 .and. lmax > 0) then
       ! Compute cross-spectrum
       ncomp = (lmax+dl_high+1)**2
       allocate(alms(1,0:lmax+dl_high,0:lmax+dl_high), weights(2*nside,1), alms1(ncomp,1), alms2(ncomp,1))
       allocate(cl1(lmin-dl_low:lmax+dl_high), cl2(lmin-dl_low:lmax+dl_high))
       weights = 1.d0
       alms    = 0.d0
       call map2alm(nside, lmax+dl_high, lmax+dl_high, my_signal(:,1), alms, [0.d0,0.d0], weights)
       !alms(:,0:lmin-1,:) = cmplx(0.d0,0.d0)
       call apodize_alms(alms, lmin, lmax, dl_low, dl_high)
       call alm2map(nside, lmax+dl_high, lmax+dl_high, alms, my_signal(:,1))
       call convert_complex_to_real_alms(alms, alms1)
       call map2alm(nside, lmax+dl_high, lmax+dl_high, residual(:,1), alms, [0.d0,0.d0], weights)
       !alms(:,0:lmin-1,:) = cmplx(0.d0,0.d0)
       call apodize_alms(alms, lmin, lmax, dl_low, dl_high)
       call alm2map(nside, lmax+dl_high, lmax+dl_high, alms, residual(:,1))
       call convert_complex_to_real_alms(alms, alms2)

       cl1 = 0.d0; cl2 = 0.d0
       i = lmin**2+1
       do l = lmin, lmax
          do m = -l, l
             cl1(l) = cl1(l) + alms1(i,1)*alms2(i,1)
             cl2(l) = cl2(l) + alms1(i,1)*alms1(i,1)
             i      = i+1
          end do
!!$          cl1(l) = cl1(l) / real(2*l+1,dp)
!!$          cl2(l) = cl2(l) / real(2*l+1,dp)
       end do

       if (g_gauss(map_id,2) == 0.d0) then
          my_gain(map_id) = bp(map_id)%gain
       else
!          my_gain(map_id) = min(max(mean(cl1/cl2), bp(map_id)%gain-MAX_DELTA_G), bp(map_id)%gain+MAX_DELTA_G)
!          my_gain(map_id) = max(min(my_gain(map_id),5.d0),0.2d0)
          my_gain(map_id) = mean(cl1(lmin:lmax)/cl2(lmin:lmax))
       end if


!!$       call int2string(map_id, band)
!!$       open(58,file='cl'//band//'.dat', recl=1024)
!!$       do l = lmin, lmax
!!$          write(58,*) l, cl1(l), cl2(l)
!!$       end do
!!$       close(58)
!!$
!!$       call int2string(map_id, band)
!!$       open(58,file='alms'//band//'.dat', recl=1024)
!!$       i = lmin**2+1
!!$       do l = lmin, lmax
!!$          do m = -l, l
!!$             write(58,*) i, alms1(i,1), alms2(i,1)
!!$             i = i+1
!!$          end do
!!$       end do
!!$       close(58)

       deallocate(cl1, cl2, alms, weights, alms1, alms2)

    else if (lmin == 0 .and. lmax == 0) then
       ! Correlate with TT plot
       numval  = count(mask_calib(:,1) > 0.5d0)
       numpair = numval*(numval-1)/2
       if (numval > 10000) then
          write(*,*) 'Error -- do not use TT plot approach with more than 10,000 pixels'
          stop
       end if
       allocate(mask2map(numval), slope(numpair))
       
       j = 1
       do i = 0, npix-1
          if (mask_calib(i,1) > 0.5d0) then
             mask2map(j) = i
             j           = j+1
          end if
       end do
       
       if (myid_alms == root) then
          if (g_gauss(map_id,2) == 0.d0) then
             my_gain(map_id) = bp(map_id)%gain
          else
             k = 1
             do i = 1, numval
                do j = i+1, numval
                   slope(k) = (residual(mask2map(i),1)-residual(mask2map(j),1)) / &
                        & (my_signal(mask2map(i),1)-my_signal(mask2map(j),1))
                   k        = k+1
                end do
             end do
             ! Only allow relatively small changes between steps, and not outside the range from 0.2 to 5
             my_gain(map_id) = median(slope)
             my_gain(map_id) = min(max(my_gain(map_id), bp(map_id)%gain-MAX_DELTA_G), bp(map_id)%gain+MAX_DELTA_G)
             my_gain(map_id) = max(min(my_gain(map_id), 5.d0),0.05d0)
             
          end if
       end if

       deallocate(mask2map, slope)

    else
       ! Correlate in pixel space with standard likelihood fit

       ! Sample gain
       invN_signal = 0.d0
       do j = 1, nmaps
          do i = 0, map_size-1
             if (mask_calib(i,j) > 0.5d0) then
                invN_signal(i,j) = my_signal(i,j) * sum(invN_rms(i,j,:))
             end if
          end do
       end do

       invN_signal = invN_signal * mask_calib
       my_sigma = sum(my_signal * invN_signal)
       my_mu    = sum(residual  * invN_signal)
       call mpi_reduce(my_mu,    mu,    1, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
       call mpi_reduce(my_sigma, sigma, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
       if (myid_alms == root) then
          if (g_gauss(map_id,2) == 0.d0) then
             my_gain(map_id) = bp(map_id)%gain
          else
             ! Compute mu and sigma from likelihood term
             mu    = mu / sigma
             sigma = sqrt(1.d0 / sigma)
             if (trim(operation) == 'optimize') sigma = 0.d0
             if (g_gauss(map_id,2) > 0.d0) then
                ! Add prior
                mu = (mu*g_gauss(map_id,2)**2 + g_gauss(map_id,1) * sigma**2) / (g_gauss(map_id,2)**2 + sigma**2)
                sigma = sqrt(sigma**2 * g_gauss(map_id,2)**2 / (sigma**2 + g_gauss(map_id,2)**2))
             end if
             gain_new = mu + sigma * rand_gauss(handle)
             ! Only allow relatively small changes between steps, and not outside the range from 0.2 to 5
             my_gain(map_id) = min(max(gain_new, bp(map_id)%gain-MAX_DELTA_G), bp(map_id)%gain+MAX_DELTA_G)
             my_gain(map_id) = max(min(my_gain(map_id),5.d0),0.05d0)

          end if
       end if
    end if

!    write(*,*) map_id, my_gain(map_id)
!!$

!!$    invN_signal = residual-bp(map_id)%gain*my_signal
!!$    call multiply_by_inv_N(invN_signal)
!!$    write(*,*) myid_chain, 'full chisq foer  = ', sum((residual-bp(map_id)%gain*my_signal)*invN_signal)
!!$
!!$    invN_signal = residual-my_gain(map_id)*my_signal
!!$    call multiply_by_inv_N(invN_signal)
!!$    write(*,*) myid_chain, 'full chisq etter = ', sum((residual-my_gain(map_id)*my_signal)*invN_signal)

!!$    if (.true. .or. map_id == 8) then
!!$       invN_signal = residual-bp(map_id)%gain*my_signal
!!$       call multiply_by_inv_N(invN_signal)
!!$       my_chisq = sum((residual-bp(map_id)%gain*my_signal)*invN_signal*mask_calib)
!!$       write(*,*) myid_chain, 'cut chisq foer  = ', chisq
!!$       call mpi_reduce(my_chisq, chisq_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)       
!!$       
!!$       invN_signal = residual-my_gain(map_id)*my_signal
!!$       call multiply_by_inv_N(invN_signal)
!!$       my_chisq = sum((residual-my_gain(map_id)*my_signal)*invN_signal*mask_calib)
!!$       write(*,*) myid_chain, 'cut chisq etter = ', sum((residual-my_gain(map_id)*my_signal)*invN_signal*mask_calib)
!!$       write(*,*) myid_chain, 'cut residual etter  = ', sum(abs(residual-bp(map_id)%gain*my_signal))
!!$    end if

    ! Collect results
    call mpi_allreduce(my_gain, bp%gain, numband, MPI_DOUBLE_PRECISION, MPI_SUM, comm_chain, ierr)

    if (output_residuals) then
       if (myid_chain == root) then
          myiter  = iter
          mychain = chain
       end if
       call mpi_bcast(myiter,  1, MPI_INTEGER, root, comm_chain, ierr)
       call mpi_bcast(mychain, 1, MPI_INTEGER, root, comm_chain, ierr)
       call int2string(map_id, band)
       call int2string(mychain, chain_text)
       call int2string(myiter,    iter_text)
       call write_map(trim(chain_dir)//'/gainres_c'//chain_text//'_b'//band//'_k'//&
            & iter_text // '.fits', residual, comptype='Residual', unit=bp(map_id)%unit, nu_ref=bp(map_id)%nu_c)
       call write_map(trim(chain_dir)//'/gainsig_c'//chain_text//'_b'//band//'_k'//&
            & iter_text // '.fits', my_signal, comptype='Signal model', &
            & unit=bp(map_id)%unit, nu_ref=bp(map_id)%nu_c)
    end if

    deallocate(residual, invN_signal, my_signal, my_gain)

!!$    call mpi_finalize(ierr)
!!$    stop

  end subroutine sample_gain



  subroutine sample_bandpass(handle, map_id, s, index_map, chain, iter)
    implicit none

    type(planck_rng),                  intent(inout)           :: handle
    integer(i4b),                      intent(in)              :: map_id
    type(genvec),                      intent(in),    optional :: s
    real(dp), dimension(0:,1:,1:),     intent(in),    optional :: index_map
    integer(i4b),                      intent(in),    optional :: iter, chain

    integer(i4b) :: i, j, k, l, m, p, ierr, ngain
    real(dp)     :: gold, gmin, gmax
    logical(lgt) :: accept
    character(len=512) :: filename
    character(len=2)   :: band, chain_text
    character(len=5)   :: iter_text
    character(len=5)   :: id_text
    logical(lgt), allocatable, dimension(:)     :: done
    real(dp),     allocatable, dimension(:)     :: delta_prop, delta0
    real(dp),     allocatable, dimension(:)     :: dgain, gain_prop, gain_old
    real(dp),     allocatable, dimension(:,:)   :: chisq_gain
    real(dp),     allocatable, dimension(:)     :: chisq0, chisq_prop, chisq_curr
    integer(i4b), dimension(MPI_STATUS_SIZE)    :: status

    gold  = 2.d0 !1.618034d0

    allocate(delta_prop(numband), chisq0(numband), chisq_prop(numband), delta0(numband), chisq_curr(numband))
    allocate(done(numband), dgain(ngain), chisq_gain(ngain,numband), gain_prop(numband), gain_old(numband))

    if (myid_chain == root) then
       call compute_bandpass_chisq(bp%delta, chisq0, s, index_map)
    else
       call compute_bandpass_chisq(bp%delta, chisq0)
    end if

    if (num_bp_step > 0) then
       accept = .false.
       do k = 1, num_bp_step
          
          delta0 = bp%delta
          if (myid_chain == root) then
             do i = 1, numband
                delta_prop(i) = delta0(i) + bp(i)%delta_rms * rand_gauss(handle)
             end do
          end if
          call mpi_bcast(delta_prop, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
          
          call update_tau(delta_prop)
          call update_eff_fg_spectrum
          if (myid_chain == root) then
             call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, index_map)
             call compute_total_chisq(map_id, s, chisq_band_fullsky=chisq_prop)
          else
             call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)
             call compute_total_chisq(map_id)
          end if
          
          ! Apply Metropolis rule for each band separately
          if (myid_chain == root) then
             do i = 1, numband
                if (bp(i)%delta_rms <= 0.d0) cycle
                if (exp(-0.5d0*(chisq_prop(i)-chisq0(i))) > rand_uni(handle)) then
                   write(*,fmt='(a,a,2f8.3,f16.3)') bp(i)%label, ' accept ', &
                        & delta0(i), delta_prop(i), chisq0(i)-chisq_prop(i)                
                   bp(i)%delta = delta_prop(i)
                   chisq0(i)   = chisq_prop(i)
                   accept      = .true.
                else
                   write(*,fmt='(a,i4,a,2f8.3,f16.3)') 'bp band = ', i, ' reject ', &
                        & delta0(i), delta_prop(i), chisq0(i)-chisq_prop(i)
                   bp(i)%delta = delta0(i)
                end if
             end do
          end if
          call mpi_bcast(bp%delta, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
          
       end do
       call mpi_bcast(accept, 1, MPI_LOGICAL, root, comm_chain, ierr)

    else if (num_bp_step < 0) then

       ! Do a rough golden ratio search
       accept     = .true.
       done       = bp%delta_rms == 0.d0
       where (.not. done)
          delta_prop = 1.d-4
       elsewhere 
          delta_prop = 0.d0
       end where

       if (myid_chain == root) then
          call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop, s, index_map)
       else
          call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop)
       end if
       where (chisq_prop > chisq0)
          delta_prop = -delta_prop
       end where

       if (myid_chain == root) then
          call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop, s, index_map)
       else
          call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop)
       end if
       where (chisq_prop > chisq0)
          delta_prop = 0.d0
          chisq_prop = chisq0
          done       = .true.
       end where

       do i = 1, 13
          if (all(done)) exit

          gain_old   = bp%gain
          chisq_curr = chisq_prop
          do j = 1, numband
             if (.not. done(j)) then
                delta_prop(j) = delta_prop(j) * gold
             end if
          end do
          if (myid_chain == root) then
             call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop, s, index_map)
             call sample_gain(handle, map_id, s, 0, 0)
          else
             call compute_bandpass_chisq(bp%delta + delta_prop, chisq_prop)
             call sample_gain(handle, map_id)
          end if
          where (done)
             chisq_prop = chisq_curr
             bp%gain    = gain_old
          end where

          if (myid_chain == root) then
             call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, index_map)
             call compute_total_chisq(map_id, s, chisq_band_fullsky=chisq_prop)
          else
             call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)
             call compute_total_chisq(map_id)
          end if
          call mpi_bcast(chisq_prop, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)

          do j = 1, numband
             if (done(j)) then
                ! Already done; do nothing
             else if (chisq_curr(j)-chisq_prop(j) < 0.d0) then
                ! Revert to previous best fit point
                delta_prop(j) = delta_prop(j) / gold
                chisq_prop(j) = chisq_curr(j)
                bp(j)%gain    = gain_old(j)
                done(j)       = .true.               ! Mark as done
             else
                chisq_curr(j) = chisq_prop(j) ! Mark this point as current best fit
             end if
          end do
          if (myid_chain == root) write(*,*) '--------------------------------------------------------------'
          if (myid_chain == root) write(*,*) 'iteration  = ', i, count(.not. done), chisq0(22)-chisq_curr(22)
          if (myid_chain == root) write(*,*) 'delta      = ', delta_prop(22)
          if (myid_chain == root) write(*,*) 'gain       = ', bp(22)%gain
       end do

       if (myid_chain == root) then
          write(*,*) 'Final corrections = ', real(delta_prop,sp)
          write(*,*) 'Final gains       = ', real(bp%gain,sp)
       end if

       ! Apply final corrections
       bp%delta = bp%delta + delta_prop


    else

       write(*,*) 'Unsupported operation mode = ', trim(operation)
       call mpi_finalize(ierr)
       stop

    end if

    if (accept) then
       ! Prepare for further evaluations
       call update_tau(bp%delta)
       call update_eff_fg_spectrum
       if (myid_chain == root) then
          call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, index_map)
       else
          call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)
       end if
    end if

    deallocate(delta_prop, chisq0, chisq_prop, delta0, done, dgain, chisq_gain, gain_prop, gain_old)

  end subroutine sample_bandpass
  

  subroutine compute_bandpass_chisq(delta, chisq, s, index_map)
    implicit none

    real(dp), dimension(:), intent(in)  :: delta
    real(dp), dimension(:), intent(out) :: chisq
    type(genvec),                      intent(in), optional :: s
    real(dp),     dimension(0:,1:,1:), intent(in), optional :: index_map
    
    ! Find search direction
    call mpi_bcast(delta, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)
    call update_tau(delta, overwrite=.false.)
    call update_eff_fg_spectrum
    if (myid_chain == root) then
       call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response, index_map)
       call compute_total_chisq(map_id, s, chisq_band_fullsky=chisq)
    else
       call update_fg_pix_response_map(map_id, pixels, fg_pix_spec_response)
       call compute_total_chisq(map_id)
    end if
    call mpi_bcast(chisq, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)

  end subroutine compute_bandpass_chisq

  subroutine apodize_alms(alms, lmin, lmax, dl_low, dl_high)
    implicit none

    complex(dpc), dimension(1:,0:,0:), intent(inout) :: alms
    integer(i4b),                      intent(in)    :: lmin, lmax, dl_low, dl_high

    integer(i4b) :: l, m

    alms(1,0:lmin-dl_low-1,:) = cmplx(0.d0, 0.d0)
!!$    do l = lmin-dl_low, lmin+dl_low
!!$       alms(1,l,:) = alms(1,l,:) * 0.5d0*(1-cos(pi*real((l-(lmin-dl_low)),dp)/real(2*dl_low,dp)))
!!$    end do
!!$
!!$    do l = lmax-dl_high, lmax+dl_high
!!$       alms(1,l,:) = alms(1,l,:) * 0.5d0*(1-cos(pi*real((l-(lmin-dl_low)),dp)/real(2*dl_low,dp)))
!!$    end do
    alms(1,lmax+1:lmax+dl_high,:) = cmplx(0.d0, 0.d0)

  end subroutine apodize_alms


end module comm_chisq_mod
