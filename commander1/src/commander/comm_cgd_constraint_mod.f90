module comm_cgd_constraint_mod
  use comm_data_mod
  use comm_N_mult_mod
  implicit none

  integer(i4b),     private :: num_V_orth_vec
  integer(i4b),     private :: comm_data, comm_chain, comm_alms, myid_alms, myid_chain, myid_data
  integer(i4b),     private :: root, ierr, seed
  type(planck_rng), private :: rng_handle
  integer(i4b)              :: num_unconst_fg_temp, num_fg_temp, num_constraints, num_fg_temp

  logical(lgt), allocatable, dimension(:),   private :: impose_fg_orthogonality

  real(dp), allocatable, dimension(:,:), private :: P, P_t

contains

  subroutine initialize_cgd_constraint_mod(paramfile, handle, comm_chain_in, comm_data_in, comm_alms_in)
    implicit none

    integer(i4b),       intent(in)    :: comm_chain_in, comm_data_in, comm_alms_in
    character(len=128), intent(in)    :: paramfile
    type(planck_rng),   intent(inout) :: handle 

    integer(i4b)       :: i, j, k, l, ind
    real(dp)           :: response, nu_eff, S_eff
    character(len=2)   :: i_text, j_text
    character(len=128) :: param_text, comp_type
    real(dp),     allocatable, dimension(:,:) :: map

    comm_data = comm_data_in
    call mpi_comm_rank(comm_data, myid_data, ierr)
    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    root = 0

    allocate(fix_free_temp(num_free_fg_temp,numband))
    allocate(impose_fg_orthogonality(num_fg_signal))

    num_constraints = 0
    do i = 1, numband
       call int2string(i, i_text)

       if (sample_T_modes) then
          param_text = 'FIX_MONOPOLE' // i_text
          call get_parameter(paramfile, trim(param_text), par_lgt=fix_free_temp(1,i))
          param_text = 'FIX_DIPOLE' // i_text
          call get_parameter(paramfile, trim(param_text), par_lgt=fix_free_temp(2,i))
          fix_free_temp(3:4,i) = fix_free_temp(2,i) 

          if (sample_free_fg_temp) then
             do j = 5, num_free_fg_temp
                call int2string(j-4, j_text)
                param_text = 'FIX_FREE_TEMP' // i_text // '_' // j_text
                call get_parameter(paramfile, trim(param_text), par_lgt=fix_free_temp(j,i))
             end do
          end if

       else

          do j = 1, num_free_fg_temp
             call int2string(j, j_text)
             param_text = 'FIX_FREE_TEMP' // i_text // '_' // j_text
             call get_parameter(paramfile, trim(param_text), par_lgt=fix_free_temp(j,i))
          end do

       end if
       

       do j = 1, num_free_fg_temp
          if (fix_free_temp(j,i)) num_constraints = num_constraints + 1
       end do
    end do

    do i = 1, num_fg_signal
       call int2string(i, i_text)

       param_text = 'COMP_TYPE' // i_text
       call get_parameter(paramfile, trim(param_text), par_string=comp_type) 
       
       param_text = 'IMPOSE_FG_ORTHOGONALITY' // i_text
       call get_parameter(paramfile, trim(param_text), par_lgt=impose_fg_orthogonality(i))

       if (trim(comp_type) == 'cmb' .and. impose_fg_orthogonality(i)) then
          if (myid_chain == root) then
             write(*,*) 'Warning: Turning off orthogonality constraints for CMB foreground component;'
             write(*,*) '         the CMB is already taken into account in the definition of the '
             write(*,*) '         orthogonality constraint (Eriksen et al. (2007) ApJ, 676, 10)'
          end if
          impose_fg_orthogonality(i) = .false.
       end if

       if (impose_fg_orthogonality(i)) num_constraints = num_constraints + num_free_fg_temp
    end do

    ! Set up global module variables
    num_constraints     = min(num_constraints, num_free_fg_temp*numband)
    num_fg_temp   = num_free_fg_temp * numband
    num_unconst_fg_temp = num_free_fg_temp * numband - num_constraints

    if (num_unconst_fg_temp > 0) then
       allocate(P(num_unconst_fg_temp, num_fg_temp))
       allocate(P_t(num_fg_temp, num_unconst_fg_temp))
    end if

    seed = 12345

    
    ! Update with default values
    call update_cgd_constraints

  end subroutine initialize_cgd_constraint_mod


  subroutine cleanup_cgd_constraint_mod
    implicit none

    if (allocated(P))                       deallocate(P)
    if (allocated(P_t))                     deallocate(P_t)
    if (allocated(fix_free_temp))           deallocate(fix_free_temp)
    if (allocated(impose_fg_orthogonality)) deallocate(impose_fg_orthogonality)

  end subroutine cleanup_cgd_constraint_mod


  subroutine update_cgd_constraints
    implicit none

    integer(i4b) :: i, j, k, l, ind
    real(dp), allocatable, dimension(:,:)   :: X, X_sum
    real(dp), allocatable, dimension(:,:)   :: map
    real(dp), allocatable, dimension(:,:)   :: A, B
    real(dp), allocatable, dimension(:,:,:) :: C, D, buffer

    if (num_unconst_fg_temp <= 0) return

    allocate(X(num_fg_temp, num_fg_temp))
    allocate(X_sum(num_fg_temp, num_fg_temp))
    allocate(map(0:map_size-1,nmaps))

    k = 1
    X = 0.d0
    do i = 1, num_fg_signal
       if (impose_fg_orthogonality(i)) then

          allocate(A(num_free_fg_temp, num_free_fg_temp))
          allocate(B(num_free_fg_temp, num_free_fg_temp))
          allocate(C(num_free_fg_temp, num_free_fg_temp, numband))
          allocate(D(num_free_fg_temp, num_free_fg_temp, numband))
          allocate(buffer(num_free_fg_temp, num_free_fg_temp, numband))
          
          A = 0.d0
          B = 0.d0
          C = 0.d0
          D = 0.d0
          do j = 1, num_free_fg_temp
             
             map = fg_temp_free(:,:,j)

             call multiply_by_inv_N(map)
             !map = map * (4.d0*pi / real(npix,dp))

             do l = 1, num_free_fg_temp
                A(j,l)        = sum(fg_temp_free(:,:,l) * map)
                C(j,l,map_id) = sum(fg_temp_free(:,:,l) * map)
                B(j,l)        = sum(fg_pix_spec_response(:,:,i) * fg_temp_free(:,:,l) * map)
                D(j,l,map_id) = sum(fg_pix_spec_response(:,:,i) * fg_temp_free(:,:,l) * map)
             end do

          end do

          call mpi_reduce(A, buffer(:,:,1), size(A), MPI_DOUBLE_PRECISION, &
               & MPI_SUM, root, comm_chain, ierr)
          if (myid_chain == root) A = buffer(:,:,1)

          call mpi_reduce(B, buffer(:,:,1), size(B), MPI_DOUBLE_PRECISION, &
               & MPI_SUM, root, comm_chain, ierr)
          if (myid_chain == root) B = buffer(:,:,1)

          call mpi_reduce(C, buffer, size(C), MPI_DOUBLE_PRECISION, &
               & MPI_SUM, root, comm_chain, ierr)
          if (myid_chain == root) C = buffer

          call mpi_reduce(D, buffer, size(D), MPI_DOUBLE_PRECISION, &
               & MPI_SUM, root, comm_chain, ierr)
          if (myid_chain == root) D = buffer


          if (myid_chain == root) then

             call invert_matrix(A)

             A = matmul(B, A)

             do j = 1, numband
                C(:,:,j) = D(:,:,j) - matmul(A, C(:,:,j))
             end do

             do j = 1, num_free_fg_temp
                do l = 1, numband
                   ind = num_free_fg_temp * (l-1)
                   X(ind+1:ind+num_free_fg_temp,k) = C(j,:,l)
                end do
                k = k+1             
             end do
          end if

          deallocate(A)
          deallocate(B)
          deallocate(C)
          deallocate(D)
          deallocate(buffer)

       end if
    end do


    if (myid_chain == root) then

       ! Normalize foreground basis vectors
       do i = 1, k-1
          X(:,i) = X(:,i) / sqrt(sum(X(:,i)*X(:,i)))
       end do

       ! Add free template constraints
       do i = 1, numband
          do j = 1, num_free_fg_temp
             if (fix_free_temp(j,i)) then
                X(num_free_fg_temp*(i-1)+j,k) = 1.d0
                k                             = k+1
             end if
          end do
       end do

       ! Add random vectors for remaining sub-space
       call rand_init(rng_handle, seed)
       do i = k, num_fg_temp
          do j = 1, num_fg_temp
             X(j,i) = rand_gauss(rng_handle)
          end do
       end do

       ! Do a Gram-Schmidt orthogonalization of the k:N subspace
       do j = 1, num_fg_temp
          do i = 1, j-1
             X(:,j) = X(:,j) - sum(X(:,j)*X(:,i)) * X(:,i)
          end do
          X(:,j) = X(:,j) / sqrt(sum(X(:,j)*X(:,j)))
       end do

       P_t = X(:,k:num_fg_temp)
       P   = transpose(P_t)

    end if

    deallocate(map)
    deallocate(X)
    deallocate(X_sum)

  end subroutine update_cgd_constraints


  subroutine apply_P_constraint(u)
    implicit none

    type(genvec), intent(inout) :: u

    integer(i4b) :: i, j, k
    real(dp), allocatable, dimension(:) :: amp

    if (num_unconst_fg_temp <= 0) then
       return
    end if

    allocate(amp(num_fg_temp))

    k = 1
    do i = 1, numband
       do j = 1, num_free_fg_temp
          amp(k) = u%temp_amp_free(j,i)
          k      = k+1
       end do
    end do

    u%temp_amp_free_const = matmul(P, amp)    

    deallocate(amp)

  end subroutine apply_P_constraint


  subroutine apply_P_trans_constraint(u)
    implicit none

    type(genvec), intent(inout) :: u

    integer(i4b) :: i, j, k
    real(dp), allocatable, dimension(:) :: amp

    if (num_unconst_fg_temp <= 0) then
       if (allocated(u%temp_amp_free)) u%temp_amp_free = 0.d0
       return
    end if

    allocate(amp(num_fg_temp))

    amp = matmul(P_t, u%temp_amp_free_const)

    k = 1
    do i = 1, numband
       do j = 1, num_free_fg_temp
          u%temp_amp_free(j,i) = amp(k)
          k                    = k+1
       end do
    end do

    deallocate(amp)

  end subroutine apply_P_trans_constraint


end module comm_cgd_constraint_mod
