module comm_cgd_precond_mod
  use comm_cgd_matmul_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  real(dp),     allocatable, dimension(:,:,:)            :: inv_N_diag
  real(dp),     allocatable, dimension(:),       private :: inv_N_lowl
  real(dp),     allocatable, dimension(:,:,:,:), private :: M_fg_pix
  real(dp),     allocatable, dimension(:,:,:),   private :: A_diagonal
  real(dp),     allocatable, dimension(:),       private :: A_lowl
  real(dp),     allocatable, dimension(:,:),     private :: A_lowl_template, A_lowl2

  integer(i4b),       private :: verbosity, numcomp_noise, precond_type
  integer(i4b),       private :: lmax_precond, numcomp_pre, lmax_noise, numcomp_block
  integer(i4b),       private :: numprocs, ierr, root = 0
  integer(i4b),       private :: myid_data, myid_alms, comm_data, comm_alms
  integer(i4b),       private :: myid_chain, comm_chain
  logical(lgt),       private :: enforce_zero_cl

contains


  ! Solve preconditioner equation, M*v = u 
  subroutine compute_invM_u(u, invM_u)
    implicit none

    type(genvec), intent(inout) :: u
    type(genvec), intent(inout) :: invM_u

    integer(i4b)     :: i, j, k, l, m, n, ind, info, nrhs, lmin_diag
    character(len=1) :: check
    real(dp), allocatable, dimension(:,:)   :: int_u, A_tmp

    call genvec_set_equal(u, invM_u)

    ! Apply constraints on free template amplitudes
    !call apply_P_trans_constraint(invM_u)

    if (mask_state == OUTSIDE_MASK) then

       if (.not. enforce_zero_cl) then
          if (precond_type == 1) then
             lmin_diag = 0
          else
             lmin_diag = lmax_precond+1
          end if
          
          ! Precondition the high-l CMB components with diagonal preconditioner
          do l = lmin_diag, lmax
             do m = -l, l
                ind = lm2ind(l,m)
                invM_u%cmb_amp(ind,:) = matmul(A_diagonal(:,:,ind), u%cmb_amp(ind,:))
             end do
          end do
       end if

       ! Precondition the low-l CMB and template coefficients
       n    = numcomp_block
       if (n > 0) then
          nrhs = 1
          allocate(int_u(n,1))
          call convert_genvec2lin_precond(u, int_u(:,1))
          
          call DPOTRS('u', n, nrhs, A_lowl2, n, int_u, n, info)
          call convert_lin_precond2genvec(int_u(:,1), invM_u)             
          deallocate(int_u)
       end if
             
    end if

    if (sample_fg_pix) then
       do i = 0, npix-1
          do j = 1, nmaps
             invM_u%fg_amp(i,j,:) = matmul(M_fg_pix(:,:,i,j), u%fg_amp(i,j,:))
          end do
       end do
    end if

    ! Apply constraints on free template amplitudes
    !call apply_P_constraint(invM_u)

  end subroutine compute_invM_u



  ! Routine for initializing the pre-conditioner before each main call to the CGD solver
  subroutine initialize_preconditioner(precond_type_in)
    implicit none

    integer(i4b), intent(in), optional :: precond_type_in

    integer(i4b)     :: i, j, k, l, m, l1, l2, m1, m2, n, ii, jj, kk, ll, q
    integer(i4b)     :: info, pos, row, col, ind1, ind2, ind, status
    real(dp)         :: t1, t2

    real(dp),     allocatable, dimension(:,:)     :: my_At_invN_A, At_invN_A_sum, A_mask, V, inv_A, alms, map
    real(dp),     allocatable, dimension(:)       :: my_At_invN_A_lowl, At_invN_A_lowl_sum, W, inv_N
    real(dp),     allocatable, dimension(:,:,:,:) :: my_M_fg_pix, M_fg_pix_sum
    integer(i4b), allocatable, dimension(:)       :: mask

    type(genvec) :: z

    if (myid_alms == root) then

       if (myid_data == root) then
          precond_type = precond_type_in
       end if
       call mpi_bcast(precond_type, 1,  MPI_LOGICAL, root, comm_data, ierr)

       if (.not. enforce_zero_cl) then
          allocate(my_At_invN_A(numcomp,nmaps))
          allocate(At_invN_A_sum(numcomp,nmaps))
          
          ! Compute the diagonal part
          my_At_invN_A = inv_N_diag(:,:,mask_state)
          call multiply_with_beam(my_At_invN_A, use_transfer=.true.)
          call multiply_with_beam(my_At_invN_A, use_transfer=.true.)
          
          call mpi_reduce(my_At_invN_A, At_invN_A_sum, size(my_At_invN_A), &
               & MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_data, ierr)    
          
          if (myid_data == root) then
             A_diagonal = 0.d0
             do l = 0, lmax
                do m = -l, l
                   ind = lm2ind(l, m)
                   do i = 1, nmaps
                      A_diagonal(i,i,ind) = At_invN_A_sum(ind,i)
                   end do
                   call compute_Lt_A_L(l, A_diagonal(:,:,ind)) 
                   do i = 1, nmaps
                      A_diagonal(i,i,ind) = 1.d0 + A_diagonal(i,i,ind)
                   end do
                   call invert_matrix(A_diagonal(:,:,ind))
                end do
             end do
          end if
          
          deallocate(my_At_invN_A, At_invN_A_sum)
       end if

       if (precond_type == 2 .and. myid_data == root) then          

          call allocate_genvec(z)

          ! Multiply with sqrt(C_l) matrix from the right
          A_lowl2 = 0.d0
          do i = 1, numcomp_block 
             call convert_lin_precond2genvec(A_lowl_template(:,i), z)
             call multiply_by_sqrt_S(.false., z%cmb_amp, lmax_int=lmax_precond)
             call convert_genvec2lin_precond(z, A_lowl2(:,i))
          end do

          ! Multiply with sqrt(C_l) matrix from the left
          do i = 1, numcomp_block 
             call convert_lin_precond2genvec(A_lowl2(i,:), z)
             call multiply_by_sqrt_S(.true., z%cmb_amp, lmax_int=lmax_precond)
             call convert_genvec2lin_precond(z, A_lowl2(i,:))
          end do

          ! Set up mask
          q = 0
          do i = 1, numcomp_block
             if (sum(abs(A_lowl2(i,:))) /= 0.d0) q = q+1
          end do
          allocate(mask(q))
          q = 0
          do i = 1, numcomp_block
             if (sum(abs(A_lowl2(i,:))) /= 0.d0) then
                q = q+1
                mask(q) = i
             end if
          end do

          ! Add unity on the diagonal of the CMB entries
          do i = 1, nmaps*numcomp_pre
             A_lowl2(i,i) = A_lowl2(i,i) + 1.d0
          end do


 !         allocate(W(numcomp_block), V(numcomp_block,numcomp_block), inv_A(numcomp_block,numcomp_block))
 !         call eigen_decomp(A_lowl2, W, V, status)
 !         if (status == 0) then
 !            open(90,file='eigenvals_precond.dat')
 !            do i = 1, numcomp_block
 !               write(90,*) i, W(i)
 !            end do
 !            close(90)
!             allocate(alms(numcomp,1), map(0:npix-1,1))
!             alms(:,1) = V(:,1)
!             call convert_harmonic_to_real_space(map, alms)
!             call write_map('eigenmode.fits', map)
             
!             do i = 1, numcomp_block
!                if (W(i)/W(numcomp_block) > 1.d-12) then
!                   W(i) = 1.d0 / W(i)
!                else
!                   W(i) = 0.d0
!                end if
!             end do
             ! Re-compose inverse matrix
!             do i = 1, numcomp_block
!                inv_A(i,:) = W(i) * V(:,i)
!             end do
!             call matmul_gen(inv_A, V, A_lowl2)
!             deallocate(V, W, inv_A)
!             call invert_singular_matrix(A_lowl2, 1.d-12)
!          else
!             write(*,*) 'Error in eigen-decomposing preconditioner'
!             stop
!          end if
!          stop

          if (q == 0) then

             ! The power spectrum is all zero 
             A_lowl2 = 0.d0
             do i = 1, numcomp_block
                A_lowl2(i,i) = 1.d0
             end do

          else

             ! Compute the Cholesky factorization
             allocate(A_mask(q,q))
             A_mask = A_lowl2(mask,mask)
             call DPOTRF('u', q, A_mask, q, info)

             if (info /= 0) then
                !             do i = 1, size(A_lowl2(:,1))
                !                write(*,*) i, real(A_lowl2(i,:),sp)
                !             end do
                write(*,*)
                allocate(W(q))
                A_mask = A_lowl2(mask,mask)
                open(58,file='eigenvals_precond.dat')
                call get_eigenvalues(A_mask, W)
                do i = 1, q
                   write(58,*) i, W(i)
                end do
                close(58)

                open(58,file='A_mat.dat', recl=10000)
                do i = 1, q
                   write(58,*) i, real(A_mask(i,:),sp)
                end do
                close(58)

                write(*,*) 
                write(*,*) 'Preconditioner matrix not positive definite. Info = ', info
                stop
             end if

             A_lowl2 = 0.d0
             do i = 1, size(A_lowl2,1)
                A_lowl2(i,i) = 1.d0
             end do

             A_lowl2(mask,mask) = A_mask
             deallocate(A_mask)

          end if

          if (allocated(mask)) deallocate(mask)
          call deallocate_genvec(z)

       end if

    end if

    ! Set up the pixel foreground preconditioner
    if (sample_fg_pix) then

       allocate(my_M_fg_pix(num_fg_signal, num_fg_signal, 0:npix-1, nmaps)) 
       allocate(M_fg_pix_sum(num_fg_signal, num_fg_signal, 0:npix-1, nmaps)) 
       allocate(inv_N(nmaps))
       
       my_M_fg_pix = 0.d0
       do k = 1, num_fg_signal
          if (fg_components(k)%enforce_positive_amplitude) cycle
          do l = 1, num_fg_signal
             if (fg_components(l)%enforce_positive_amplitude) cycle
             do i = 0, map_size-1
                call get_inv_N_sub_matrix(map_id, mask_state, i, inv_N)
                do j = 1, nmaps
                   my_M_fg_pix(k,l,pixels(i),j) = fg_pix_spec_response(i,j,k) * &
                        & inv_N(j) * fg_pix_spec_response(i,j,l)
                end do
             end do
                
          end do
       end do
       deallocate(inv_N)

       call mpi_reduce(my_M_fg_pix, M_fg_pix_sum, size(my_M_fg_pix), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_chain, ierr)           


       if (myid_chain == root) then

          M_fg_pix = M_fg_pix_sum
          do i = 0, npix-1
             do j = 1, nmaps
                do l = 1, num_fg_signal
                   if (M_fg_pix(l,l,i,j) <= 0.d0) then 
                      M_fg_pix(l,:,i,j) = 0.d0
                      M_fg_pix(:,l,i,j) = 0.d0
                      M_fg_pix(l,l,i,j) = 1.d0
                   end if
                end do
                call invert_matrix(M_fg_pix(:,:,i,j))
             end do
          end do
       end if

       deallocate(my_M_fg_pix)
       deallocate(M_fg_pix_sum)

    end if

  end subroutine initialize_preconditioner



  ! ******************************************
  ! *          Initialization routines       *
  ! ******************************************


  ! Routine for pre-computing preconditioner quantities
  subroutine pre_initialize_preconditioner(chain_dir, comm_chain_in, comm_data_in, comm_alms_in, comm_master, &
       & paramfile)
    implicit none

    integer(i4b),        intent(in)  :: comm_chain_in, comm_data_in, comm_alms_in, comm_master
    character(len=128),  intent(in)  :: paramfile, chain_dir

    integer(i4b) :: i, j, k, l, m, ii, jj, row, col, pos, incx, incy, ind, unit
    real(dp)     :: alpha, beta, t1, t2
    logical(lgt) :: polarization, exist, ok
    character(len=128) :: method, filename

    real(dp),     allocatable, dimension(:,:,:)   :: alms_noise
    real(dp),     allocatable, dimension(:,:)     :: fg_alms
    real(dp),     allocatable, dimension(:,:)     :: map
    real(dp),                  dimension(1:2)     :: zbounds
    real(dp),     allocatable, dimension(:,:)     :: weights, inv_N
    complex(dpc), allocatable, dimension(:,:,:)   :: tmp_alms
    integer(i4b), dimension(MPI_STATUS_SIZE)      :: status

    unit          = getlun()
    comm_data     = comm_data_in
    comm_alms     = comm_alms_in
    comm_chain    = comm_chain_in
    
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    call mpi_comm_rank(comm_data,  myid_data,  ierr)
    call mpi_comm_rank(comm_alms,  myid_alms,  ierr)
    call mpi_comm_size(comm_alms,  numprocs,   ierr)

    if (myid_alms == root) then
       
       call get_parameter(paramfile, 'VERBOSITY',             par_int=verbosity)
       call get_parameter(paramfile, 'LMAX_PRECOND',          par_int=lmax_precond)
       call get_parameter(paramfile, 'POLARIZATION',          par_lgt=polarization)
       call get_parameter(paramfile, 'PRECOND_MATRIX_METHOD', par_string=method)
       call get_parameter(paramfile, 'ENFORCE_ZERO_CL',       par_lgt=enforce_zero_cl)
       if (enforce_zero_cl) then
          lmax_precond = -1
          numcomp      = 0
          numcomp_pre  = 0
       else
          lmax_precond = max(lmax_precond,2)
          if (lmax_precond > lmax) lmax_precond = lmax
          numcomp         = lmax2numcomp(lmax)
          numcomp_pre     = lmax2numcomp(lmax_precond)
       end if
       numcomp_block = nmaps*numcomp_pre + num_fg_temp*numband

       allocate(inv_N_lowl(numcomp_pre*(numcomp_pre+1)/2))
       allocate(inv_N_diag(1:numcomp,nmaps,2))

       allocate(A_diagonal(nmaps, nmaps, 1:numcomp))
       allocate(A_lowl2(numcomp_block, numcomp_block))
       allocate(A_lowl_template(numcomp_block, numcomp_block))
       if (sample_fg_pix) allocate(M_fg_pix(num_fg_signal, num_fg_signal, 0:npix-1, nmaps))

    end if

    call mpi_bcast(lmax,             1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(npix,             1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(nmaps,            1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(nside,            1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(numcomp,          1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(lmax_precond,     1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(numcomp_pre,      1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(numcomp_block,    1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(verbosity,        1,  MPI_INTEGER,   root, comm_alms, ierr)
    call mpi_bcast(polarization,     1,  MPI_LOGICAL,   root, comm_alms, ierr)
    call mpi_bcast(enforce_zero_cl,  1,  MPI_LOGICAL,   root, comm_alms, ierr)
    call mpi_bcast(method,         128,  MPI_CHARACTER, root, comm_alms, ierr)


    ! Pre-compute noise alms
    lmax_noise    = 2*lmax+5
    numcomp_noise = lmax2numcomp(lmax_noise)
    zbounds       = 0.d0

    allocate(alms_noise(numcomp_noise, nmaps, 2))

    ! Compute a_lms of inverse noise covariance maps
    if (myid_alms == root) then

       allocate(tmp_alms(1,0:lmax_noise,0:lmax_noise))
       allocate(weights(1:2*nside, 1))
!       call read_ringweights(nside, weights)
       weights = 1.d0

       allocate(inv_N(0:npix-1, nmaps))

       do i = 1, 2
          call get_noise_map(i, inv_N)
          
          ! Convert to thermodynamic units
          inv_N = inv_N * (spec2data(map_id, cmb=.true.)/bp(map_id)%gain)**2

          call map2alm(nside, lmax_noise, lmax_noise, inv_N(:,1), &
               & tmp_alms, zbounds, weights)
          
          call convert_complex_to_real_alms(tmp_alms, alms_noise(:,1:1,i))
          
          if (nmaps == 3) then
             inv_N(:,1) = 0.5d0*(inv_N(:,2)+inv_N(:,3))
             call map2alm(nside, lmax_noise, lmax_noise, inv_N(:,1), &
                  & tmp_alms, zbounds, weights)
             call convert_complex_to_real_alms(tmp_alms, alms_noise(:,2:2,i))
             call convert_complex_to_real_alms(tmp_alms, alms_noise(:,3:3,i))
          end if
       end do

       deallocate(inv_N)
       deallocate(tmp_alms)
       deallocate(weights)

    else

       do i = 1, 2
          call get_noise_map(i)
       end do

    end if

    call mpi_bcast(alms_noise,  size(alms_noise),  MPI_DOUBLE_PRECISION, &
         & root, comm_alms, ierr)


    ! Pre-compute high-l noise matrices in spherical harmonic space
    if (.not. enforce_zero_cl) then
       do i = 1, 2
          call initialize_inv_N_diag(i, alms_noise(:,:,i))
       end do
    end if

    ! Precompute low-l + template block
    if (myid_alms == root .and. myid_data == root) then
       filename = trim(chain_dir) // '/precond_mat.unf'
       ok       = .false.
       inquire(file=trim(filename), exist=exist)
       if (exist) then
          open(unit,file=trim(filename),form='unformatted')
          read(unit) i
          ok = (i == size(A_lowl_template))
          if (ok) read(unit) A_lowl_template
          close(unit)
       end if
       if (.not. ok) then
          call initialize_inv_N_lowl(comm_master)
          open(unit,file=trim(filename),form='unformatted')
          write(unit) size(A_lowl_template)
          write(unit) A_lowl_template
          close(unit)
       end if
    end if

    if (myid_chain /= root .and. allocated(A_lowl2))         deallocate(A_lowl2)
    if (myid_chain /= root .and. allocated(A_lowl_template)) deallocate(A_lowl_template)
    deallocate(alms_noise)

  end subroutine pre_initialize_preconditioner


  subroutine cleanup_precond_mod
    implicit none

    if (allocated(A_diagonal))          deallocate(A_diagonal)
    !if (allocated(A_lowl))              deallocate(A_lowl)
    if (allocated(A_lowl2))             deallocate(A_lowl2)
    if (allocated(A_lowl_template))     deallocate(A_lowl_template)
    if (allocated(inv_N_diag))          deallocate(inv_N_diag)
    if (allocated(inv_N_lowl))          deallocate(inv_N_lowl)
    if (allocated(M_fg_pix))            deallocate(M_fg_pix)

  end subroutine cleanup_precond_mod



  ! Routine for computing the diagonal elements of the inverse noise matrix, ie., setting
  ! up the diagonal part of the preconditioner
  subroutine initialize_inv_N_diag(mask_state, alms_noise)
    implicit none

    integer(i4b),                   intent(in) :: mask_state
    real(dp),     dimension(1:,1:), intent(in) :: alms_noise

    real(dp)     :: one_over_sqrt_4pi, l0min, l0max, l1min, l1max
    integer(i4b) :: ind1, ind2, ier, twolmaxp2
    integer(i4b) :: pos, row, col
    integer(i4b) :: i, j, k, l, m, lp, l1, l2, l3, m1, m2, m3, ind
    complex(dpc) :: val
    real(dp),     allocatable, dimension(:)     :: threej_symbols, threej_symbols_m0, sqrt_2lp1

    integer(i4b), dimension(MPI_STATUS_SIZE) :: status

    real(dp),     allocatable, dimension(:,:)     :: my_inv_N_diag
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_noise_cmplx

    twolmaxp2     = 2*lmax_noise+2

    allocate(threej_symbols(twolmaxp2))
    allocate(threej_symbols_m0(twolmaxp2))
    allocate(sqrt_2lp1(0:lmax_noise))
    allocate(alms_noise_cmplx(nmaps,0:lmax_noise,0:lmax_noise))

    call convert_real_to_complex_alms(alms_noise, alms_noise_cmplx)

    ! Precompute square root of 2l+1
    do l1 = 0, lmax_noise
       sqrt_2lp1(l1) = sqrt(2.d0*real(l1,dp)+1.d0)
    end do
    one_over_sqrt_4pi = 1.d0 / sqrt(4.d0*pi)

    ! Compute the high-l diagonal elements
    allocate(my_inv_N_diag(1:numcomp,nmaps))
    my_inv_N_diag = 0.d0

    do l = myid_alms, lmax, numprocs
       
       call DRC3JJ(real(l,dp), real(l,dp), 0.d0, 0.d0, l0min, l0max, &
            & threej_symbols_m0, twolmaxp2, ier)
             
       do m = 0, l
             
          call DRC3JJ(real(l,dp), real(l,dp), real(-m,dp), real(m,dp), l1min, l1max, &
               & threej_symbols, twolmaxp2, ier)
             
          do i = 1, nmaps

             val = cmplx(0.d0,0.d0)
             do l2 = 1, nint(l1max-l1min)+1
                lp = nint(l1min) + l2 - 1
                if (lp > lmax_noise) exit
                ind = lm2ind(lp,0)
                
                val = val + alms_noise_cmplx(i,lp,0) * sqrt_2lp1(lp) * &
                     & threej_symbols(l2) * threej_symbols_m0(lp-nint(l0min)+1)
             end do
             
             val = val * sqrt_2lp1(l)**2 * one_over_sqrt_4pi
             val = val * real(12*nside**2,dp) / (4.d0*pi)
             
             if (mod(m,2)/=0) val = -val
             
             ind = lm2ind(l,m)
             my_inv_N_diag(ind,i) = real(val,dp)

             ind = lm2ind(l,-m)
             my_inv_N_diag(ind,i) = real(val,dp)

          end do

       end do
    end do

!    call mpi_barrier(comm_alms, ierr)

    if (myid_alms == root) then
    
       ! Insert my own data
       inv_N_diag(:,:,mask_state) = my_inv_N_diag
          
       do j = 1, numprocs-1

          call mpi_recv(my_inv_N_diag, size(my_inv_N_diag), MPI_DOUBLE_PRECISION, j, &
               & 61, comm_alms, status, ierr)    
          
          do l = j, lmax, numprocs
             do m = -l, l
                ind = lm2ind(l,m)
                inv_N_diag(ind,:,mask_state) = my_inv_N_diag(ind,:)
             end do
          end do

       end do

    else
       
       call mpi_send(my_inv_N_diag, size(my_inv_N_diag), MPI_DOUBLE_PRECISION, &
            & root, 61, comm_alms, ierr)
       
    end if

    deallocate(my_inv_N_diag)
    deallocate(sqrt_2lp1)
    deallocate(threej_symbols)
    deallocate(threej_symbols_m0)
    deallocate(alms_noise_cmplx)

  end subroutine initialize_inv_N_diag



  ! Routine for computing a low-l block. Used for setting up the low-l part of
  ! the preconditioner
  subroutine initialize_inv_N_lowl(comm_master)
    implicit none

    integer(i4b), intent(in) :: comm_master

    integer(i4b) :: ind1, ind2, myid_master, numprocs_master
    integer(i4b) :: ierr, pos, row, col
    integer(i4b) :: i, j, k, l, m, l1, l2, l3, m1, m2, m3
    real(dp)     :: t1, t2
    integer(i4b), dimension(MPI_STATUS_SIZE) :: status

    real(dp), allocatable, dimension(:)   :: z
    real(dp), allocatable, dimension(:,:) :: map
    real(dp), allocatable, dimension(:,:) :: alms, buffer
    type(genvec) :: x, y

    call mpi_comm_rank(comm_master, myid_master, ierr)
    call mpi_comm_size(comm_master, numprocs_master, ierr)

    call allocate_genvec(x)
    call allocate_genvec(y)
    allocate(z(numcomp_block))

    ! Initialize to zero
    A_lowl_template = 0.d0
    
    ! Compute CMB columns
    col = 0
    do i = 1, nmaps
       do l = 0, lmax_precond
          do m = -l, l
             col = col+1
             if (mod(col-1+myid_master,numprocs_master) /= 0) cycle
             if (.not. enforce_zero_cl) then

                if (m == -l) then
                   write(*,fmt='(a,i4,a,i4)') 'cgd_precond_mod -- computing A_temp for l = ', l, ' and i = ', i
                end if

                call nullify_genvec(x)
                ind1              = lm2ind(l,m)
                x%cmb_amp(ind1,i) = 1.d0
                call compute_Au(x, y)

                y%cmb_amp(ind1,i) = y%cmb_amp(ind1,i) - 1.d0 ! Subtract prior term

                ! Store result
                call convert_genvec2lin_precond(y, z)
                A_lowl_template(:,col) = z
             end if

          end do
       end do
    end do

    ! Compute foreground template columns
    do i = 1, numband
       do j = 1, num_fg_temp
          col = col+1
          if (mod(col-1+myid_master,numprocs_master) /= 0) cycle
          
          call nullify_genvec(x)
          x%temp_amp(j,i) = 1.d0
          call compute_Au(x, y)
          
          ! Store result
          call convert_genvec2lin_precond(y, z)
          A_lowl_template(:,col) = z
       end do
    end do


    ! Send all elements to all chains
    allocate(buffer(numcomp_block,numcomp_block))
    call mpi_allreduce(A_lowl_template, buffer, size(buffer), MPI_DOUBLE_PRECISION, MPI_SUM, comm_master, ierr)
    A_lowl_template = buffer
    deallocate(buffer)
    
    ! Divide out power spectrum
    do i = 1, numcomp_block
       call nullify_genvec(x)
       call convert_lin_precond2genvec(A_lowl_template(i,:), x)
       call multiply_by_inv_sqrt_S(.false., x%cmb_amp)
       call convert_genvec2lin_precond(x, A_lowl_template(i,:))
    end do
    
    do i = 1, numcomp_block
       call nullify_genvec(x)
       call convert_lin_precond2genvec(A_lowl_template(:,i), x)
       call multiply_by_inv_sqrt_S(.true., x%cmb_amp)
       call convert_genvec2lin_precond(x, A_lowl_template(:,i))
    end do

!!$    ! Check for numerically singular modes
!!$    do i = 1, numcomp_block
!!$       if (abs(A_lowl_template(i,i)) < 1.d-8) then
!!$          A_lowl_template(i,:) = 0.d0
!!$          A_lowl_template(:,i) = 0.d0
!!$       end if
!!$    end do
!!$
!!$    ! Check for symmetric matrix
!!$    do i = 1, numcomp_block
!!$       do j = i+1, numcomp_block
!!$          if (A_lowl_template(i,j) == 0.d0 .and. A_lowl_template(j,i) == 0.d0) cycle
!!$          if (abs(A_lowl_template(i,j)-A_lowl_template(j,i)) / (A_lowl_template(i,j)+A_lowl_template(j,i)) > 1.d-0 .and. &
!!$               & abs(A_lowl_template(i,j)) / sqrt(A_lowl_template(i,i)*A_lowl_template(j,j)) > 1.d-0) then
!!$             write(*,fmt='(a,i5,i5)') 'comm_cgd_precond_mod -- Asymmetric preconditioner for (i,j):', i, j
!!$             write(*,*) '     A(i,j) = ', A_lowl_template(i,j)
!!$             write(*,*) '     A(j,i) = ', A_lowl_template(j,i)
!!$             write(*,*) '     A(i,i) = ', A_lowl_template(i,i)
!!$             write(*,*) '     A(j,j) = ', A_lowl_template(j,j)
!!$             write(*,*) '     Fractional difference = ', abs(A_lowl_template(i,j)-A_lowl_template(j,i)) / &
!!$                  & (A_lowl_template(i,j)+A_lowl_template(j,i))
!!$             write(*,*) '     Fraction of diagonal  = ', abs(A_lowl_template(i,j)) / sqrt(A_lowl_template(i,i)*A_lowl_template(j,j))
!!$             stop
!!$          end if
!!$       end do
!!$    end do

  end subroutine initialize_inv_N_lowl



  subroutine convert_genvec2lin_precond(x, y)
    implicit none

    type(genvec),                intent(in)  :: x
    real(dp),     dimension(1:), intent(out) :: y

    integer(i4b) :: i, j, k

    ! Set up CMB part
    do i = 1, nmaps
       y((i-1)*numcomp_pre+1:i*numcomp_pre) = x%cmb_amp(1:numcomp_pre,i)
    end do
    
    ! Set up free template part
    k = nmaps*numcomp_pre + 1
    do j = 1, numband
       do i = 1, num_fg_temp
          y(k) = x%temp_amp(i,j)
          k    = k+1
       end do
    end do

  end subroutine convert_genvec2lin_precond

  subroutine convert_lin_precond2genvec(y, x)
    implicit none

    real(dp),     dimension(1:), intent(in)    :: y
    type(genvec),                intent(inout) :: x

    integer(i4b) :: i, j, k

    ! Set up CMB part
    do i = 1, nmaps
       x%cmb_amp(1:numcomp_pre,i) = y((i-1)*numcomp_pre+1:i*numcomp_pre)
    end do

    ! Set up free template part
    k = nmaps*numcomp_pre + 1
    do j = 1, numband
       do i = 1, num_fg_temp
          x%temp_amp(i,j) = y(k)
          k = k+1
       end do
    end do

  end subroutine convert_lin_precond2genvec


end module comm_cgd_precond_mod
 
