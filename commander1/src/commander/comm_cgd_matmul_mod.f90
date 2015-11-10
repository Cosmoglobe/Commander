module comm_cgd_matmul_mod
  use comm_data_mod
  use comm_S_mult_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b),       private :: myid_alms, comm_alms, root, myid_chain, comm_chain
  logical(lgt),       private :: enforce_zero_cl
  type(planck_rng),   private :: rng_handle
  character(len=128), private :: operation
  
contains

  ! Driver routine
  subroutine compute_Au(vec_in, vec_out)
    implicit none

    type(genvec), intent(inout) :: vec_in
    type(genvec), intent(inout) :: vec_out

    integer(i4b) :: ierr

    ! Multiply with beam and inverse noise matrix
    call mpi_bcast(1, 1, MPI_INTEGER, root, comm_chain, ierr)
    call compute_At_invN_A_x(vec_in, vec_out)

  end subroutine compute_Au


  ! **********************************************************************************

  subroutine initialize_cgd_matmul_mod(paramfile, comm_alms_in, comm_chain_in, seed, &
       & npix_in, nmaps_in, lmax_in)
    implicit none

    integer(i4b),       intent(in) :: comm_alms_in, comm_chain_in
    character(len=128), intent(in) :: paramfile
    integer(i4b),       intent(in) :: seed, npix_in, nmaps_in, lmax_in

    integer(i4b) :: ierr

    lmax     = lmax_in
    npix     = npix_in
    nmaps    = nmaps_in
    numcomp  = lmax2numcomp(lmax)
       
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL', par_lgt=enforce_zero_cl)
    call get_parameter(paramfile, 'OPERATION',       par_string=operation)

    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    root = 0

    call rand_init(rng_handle, seed)

  end subroutine initialize_cgd_matmul_mod


  subroutine clean_up_cgd_matmul_mod
    implicit none

  end subroutine clean_up_cgd_matmul_mod



  subroutine compute_signal_rhs(add_wf_in, add_S_omega_in, add_N_omega_in, rhs, fluct)
    implicit none

    logical(lgt), intent(in),    optional :: add_wf_in, add_S_omega_in, add_N_omega_in
    type(genvec), intent(inout), optional :: rhs
    type(genvec), intent(inout), optional :: fluct

    integer(i4b) :: i, j, k, l, m, ierr
    real(dp)     :: omega
    real(dp), allocatable, dimension(:,:)   :: map, map_sum, alms, eta
    logical(lgt) :: add_wf, add_S_omega, add_N_omega
    type(genvec) :: v1, v2

    call allocate_genvec(v1)
    call allocate_genvec(v2)
    call nullify_genvec(v1)
    if (myid_chain == root) then
       call nullify_genvec(rhs)
       add_wf      = add_wf_in
       add_S_omega = add_S_omega_in
       add_N_omega = add_N_omega_in
    end if
    call mpi_bcast(add_wf,      1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(add_S_omega, 1, MPI_LOGICAL, root, comm_chain, ierr)
    call mpi_bcast(add_N_omega, 1, MPI_LOGICAL, root, comm_chain, ierr)

    allocate(map(0:map_size-1, nmaps), eta(0:map_size-1,nmaps))
    allocate(map_sum(0:map_size-1, nmaps))
    if (add_wf) then
       call multiply_by_inv_N(residual, map_sum)
    else
       map_sum = 0.d0
    end if

    if (.not. trim(operation) == 'optimize' .and. add_N_omega) then
       ! Draw the random field in pixel space
       do i = 1, nmaps
          do j = 0, map_size-1
             eta(j,i)  = rand_gauss(rng_handle)
          end do
       end do
       call multiply_by_sqrt_inv_N(eta)
       map_sum = map_sum + eta
    end if

    ! Add CMB terms
    if (.not. enforce_zero_cl .and. mask_state == OUTSIDE_MASK) then
       if (myid_alms == root) then
          allocate(alms(numcomp,nmaps))
          call convert_real_to_harmonic_space(map_sum, alms)
          call multiply_with_beam(alms, transpose=.true.)
          v1%cmb_amp = spec2data(map_id, cmb=.true.) * alms * real(npix,dp)/(4.d0*pi)
          deallocate(alms)
       else
          call convert_real_to_harmonic_space(map_sum)
       end if
    end if

    ! Add template terms; only outside mask
    if (mask_state == OUTSIDE_MASK .and. sample_temp_coeffs) then
       do i = 1, num_fg_temp
          if (.not. fix_temp(i,map_id)) v1%temp_amp(i,map_id) = sum(fg_temp(:,:,i) * map_sum)
       end do
    end if

    ! The pixel foreground part
    do i = 1, num_fg_signal
       if (fg_components(i)%enforce_positive_amplitude) cycle
       if (.not. enforce_zero_cl .and. trim(fg_components(i)%type) == 'cmb') cycle
       v1%fg_amp(pixels,:,i) = fg_pix_spec_response(:,:,i) * map_sum
    end do

    ! Collect results from all processors
    call reduce_genvec(comm_chain, v1, v2)

    if (myid_chain == root) then

       call multiply_by_sqrt_S(.true., v2%cmb_amp)
       call genvec_set_equal(v2, rhs)

       ! Add CMB prior term
       if ((.not. enforce_zero_cl) .and. mask_state == OUTSIDE_MASK .and. &
            & (trim(operation) /= 'optimize') .and. add_S_omega) then
          do j = 1, nmaps
             do i = 1, numcomp
                omega = rand_gauss(rng_handle)
                rhs%cmb_amp(i,j) = rhs%cmb_amp(i,j) + omega
                if (present(fluct)) fluct%cmb_amp(i,j) = omega
             end do
          end do
       end if

       ! Apply constraints on free template amplitudes
       !call apply_P_constraint(rhs)
    
    end if

!!$       do l = 2, 3
!!$          do m = -l, l
!!$             i = lm2ind(l,m)
!!$             write(*,*) rhs%cmb_amp(i,1)
!!$          end do
!!$       end do
!!$
!!$    
!!$    open(58,file='rhs_old.dat')
!!$    do i = 1, size(rhs%cmb_amp,1)
!!$       write(58,*) rhs%cmb_amp(i,1)
!!$    end do
!!$    close(58)
!!$    call mpi_finalize(ierr)
!!$    stop
    
    deallocate(map, map_sum, eta)
    call deallocate_genvec(v1)
    call deallocate_genvec(v2)

  end subroutine compute_signal_rhs



  subroutine compute_At_invN_A_x(coeff_in, coeff_out)
    implicit none

    type(genvec), intent(inout), optional :: coeff_in
    type(genvec), intent(inout), optional :: coeff_out

    real(dp)     :: t1, t2
    integer(i4b) :: i, j, k, l
    type(genvec) :: v1, v2
    real(dp), allocatable, dimension(:,:) :: map, map_sum, alms1

    call wall_time(t1)

    call allocate_genvec(v1)
    call allocate_genvec(v2)
    if (myid_chain == root) then
       call genvec_set_equal(coeff_in, v1)
       !call apply_P_trans_constraint(v1)
       call multiply_by_sqrt_S(.false., v1%cmb_amp)
    end if
    call bcast_genvec(comm_chain, v1)
    call nullify_genvec(v2)

    allocate(map(0:map_size-1,nmaps), map_sum(0:map_size-1,nmaps))
    map_sum = 0.d0

    ! ========================================
    ! Add signal terms
    !========================================
    ! Add CMB term
    if (.not. enforce_zero_cl .and. mask_state == OUTSIDE_MASK) then
       if (myid_alms == root) then
          allocate(alms1(numcomp,nmaps))
          alms1 = spec2data(map_id, cmb=.true.) * v1%cmb_amp
          call multiply_with_beam(alms1)
          call wall_time(t1)
          call convert_harmonic_to_real_space(map, alms1)
          call wall_time(t2)
!          write(*,*) 'Wall SHT1 = ', t2-t1
          deallocate(alms1)
       else
          call convert_harmonic_to_real_space(map)
       end if
       map_sum = map_sum + map 
    end if

    if (mask_state == OUTSIDE_MASK .and. sample_temp_coeffs) then
       call wall_time(t1)
       ! Add free template term
       do l = 1, num_fg_temp
          if (.not. fix_temp(l,map_id)) map_sum = map_sum + v1%temp_amp(l,map_id) * fg_temp(:,:,l)
       end do
       call wall_time(t2)
!       write(*,*) 'Wall temp1 = ', t2-t1
    end if

    ! Add pixel foreground term
    do l = 1, num_fg_signal
       if (fg_components(l)%enforce_positive_amplitude) cycle
       if (.not. enforce_zero_cl .and. trim(fg_components(l)%type) == 'cmb') cycle
       map_sum = map_sum + fg_pix_spec_response(:,:,l) * v1%fg_amp(pixels,:,l)
    end do

    call wall_time(t1)
    call multiply_by_inv_N(map_sum)
    call wall_time(t2)
!    write(*,*) 'Wall invN = ', t2-t1


    ! ========================================
    ! Collect results
    !========================================

    ! Collect CMB term
    if (.not. enforce_zero_cl .and. mask_state == OUTSIDE_MASK) then
       if (myid_alms == root) then
          allocate(alms1(numcomp,nmaps))
          call wall_time(t1)
          call convert_real_to_harmonic_space(map_sum, alms1)
          call wall_time(t2)
!          write(*,*) 'Wall SHT2 = ', t2-t1
          call multiply_with_beam(alms1, transpose=.true.)
          v2%cmb_amp = v2%cmb_amp + spec2data(map_id, cmb=.true.) * alms1 * real(npix,dp) / (4.d0*pi)
          deallocate(alms1)
       else
          call convert_real_to_harmonic_space(map_sum)
       end if
    end if

    if (mask_state == OUTSIDE_MASK .and. sample_temp_coeffs) then
       call wall_time(t1)
       ! Collect free template term
       do i = 1, num_fg_temp
          if (.not. fix_temp(i,map_id)) then
             v2%temp_amp(i,map_id) = v2%temp_amp(i,map_id) + sum(fg_temp(:,:,i)*map_sum)
          end if
       end do
       call wall_time(t2)
!       write(*,*) 'Wall temp2 = ', t2-t1
    end if

    ! The pixel foreground part
    do i = 1, num_fg_signal
       if (fg_components(i)%enforce_positive_amplitude) cycle
       if (.not. enforce_zero_cl .and. trim(fg_components(i)%type) == 'cmb') cycle
       v2%fg_amp(pixels,:,i) = v2%fg_amp(pixels,:,i) + fg_pix_spec_response(:,:,i) * map_sum
    end do

    ! Collect results from all processors
    call wall_time(t1)
    call reduce_genvec(comm_chain, v2, v1)
    call wall_time(t2)
!    write(*,*) 'Wall reduce = ', t2-t1

    ! Apply constraints on free template amplitudes
    if (myid_chain == root) then
       call genvec_set_equal(v1, coeff_out)
       !call apply_P_constraint(coeff_out)

       ! Add unit vector from prior
       call multiply_by_sqrt_S(.true., coeff_out%cmb_amp)
       coeff_out%cmb_amp = coeff_out%cmb_amp + coeff_in%cmb_amp
    end if

    deallocate(map, map_sum)
    call deallocate_genvec(v1)
    call deallocate_genvec(v2)

    call wall_time(t2)
    !if (myid_chain == root) write(*,*) 'CPU time = ', t2-t1

  end subroutine compute_At_invN_A_x

end module comm_cgd_matmul_mod

! LocalWords:  rhs
