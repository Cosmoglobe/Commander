module comm_N_mult_mod
  use healpix_types
  use comm_utils
  use rngmod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  integer(i4b), parameter :: OUTSIDE_MASK = 1
  integer(i4b), parameter :: INSIDE_MASK  = 2

  integer(i4b) :: mask_state

  integer(i4b),                              private :: comm_alms, myid_alms, root, numprocs, numband
  integer(i4b),                              private :: comm_chain, myid_chain, map_id
  integer(i4b),                              private :: nside, npix, nmaps, map_size, lmax, numcomp
  integer(i4b),                              private :: inv_N_order
  integer(i8b),                              private :: col_from, col_to, numpix
  type(planck_rng),                          private :: rng_handle

  character(len=128)                                     :: noise_format
  real(dp),     allocatable, dimension(:,:,:)            :: invN_rms, sqrt_invN_rms
  real(dp),     allocatable, dimension(:,:),     private :: invN_dense, sqrt_invN_dense
  real(dp),     allocatable, dimension(:,:,:),   private :: invN_pixblock, sqrt_invN_pixblock
  integer(i4b), pointer,     dimension(:),       private :: pixels
  integer(i4b), allocatable, dimension(:),       private :: mask2map
  real(dp),     allocatable, dimension(:)                :: noiseamp

  interface multiply_by_inv_N
     module procedure multiply_by_inv_N_one_map, multiply_by_inv_N_two_maps
  end interface 


  interface multiply_by_sqrt_inv_N
     module procedure multiply_by_sqrt_inv_N_one_map, multiply_by_sqrt_inv_N_two_maps
  end interface 

contains

  subroutine initialize_N_mult_mod(paramfile, comm_alms_in, comm_chain_in, map_id_in)
    implicit none

    integer(i4b),       intent(in) :: map_id_in, comm_alms_in, comm_chain_in
    character(len=128), intent(in) :: paramfile

    integer(i4b)       :: i, j, k, l, m, n, ind, ierr, map, pix, unit
    logical(lgt)       :: polarization, sample_inside_mask
    real(dp)           :: reg_noise, reg_scale
    character(len=2)   :: map_text, i_text
    character(len=128) :: noisefile, maskfile, paramtext
    real(dp),     allocatable, dimension(:,:) :: noisemap, map_in
    real(dp),     allocatable, dimension(:,:) :: mask

    unit = getlun()
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    call mpi_comm_size(comm_alms, numprocs, ierr)
    root = 0

    mask_state = OUTSIDE_MASK

    ! Read general parameters
    map_id = map_id_in
    call int2string(map_id, map_text)
    call get_parameter(paramfile, 'MASKFILE',                  par_string=maskfile)
    call get_parameter(paramfile, 'NOISE_FORMAT',              par_string=noise_format)       
    call get_parameter(paramfile, 'NSIDE',                     par_int=nside)       
    call get_parameter(paramfile, 'LMAX',                      par_int=lmax)       
    call get_parameter(paramfile, 'POLARIZATION',              par_lgt=polarization)       
    call get_parameter(paramfile, 'REGULARIZATION_NOISE',      par_dp=reg_noise)       
    call get_parameter(paramfile, 'REG_NOISE_SCALE'//map_text, par_dp=reg_scale)       
    call get_parameter(paramfile, 'SAMPLE_INSIDE_MASK',        par_lgt=sample_inside_mask)       
    call get_parameter(paramfile, 'NUMBAND',                   par_int=numband)       
    reg_noise = reg_noise * reg_scale
    npix    = 12*nside**2
    numcomp = (lmax+1)**2

    if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if

    ! Initialize noise scalings
    allocate(noiseamp(numband))
    noiseamp = 1.d0

    ! Find which pixels to include for this processor
    call get_pixels(pixels)
    map_size = size(pixels)
 
    ! Read mask file
    allocate(mask(0:npix-1,nmaps))
    call read_map(maskfile, mask)

    ! Read noise RMS maps -- this is always needed
    paramtext = 'NOISE_RMS' // map_text
    call get_parameter(paramfile, trim(paramtext), par_string=noisefile)
    
    allocate(invN_rms(0:map_size-1,nmaps,2))
    allocate(sqrt_invN_rms(0:map_size-1,nmaps,2))
    
    allocate(noisemap(0:npix-1,nmaps))
    call read_map(noisefile, noisemap)
    if (reg_noise > 0.d0) then
       noisemap = reg_noise
    else if (reg_noise < 0.d0) then
       noisemap = sqrt(noisemap**2 + reg_noise**2)
    end if

    invN_rms = 0.d0
    do j = 1, nmaps
       do i = 0, map_size-1
          if (mask(pixels(i),j) == 1.d0) then
             invN_rms(i,j,OUTSIDE_MASK) = 1.d0 / noisemap(pixels(i),j)**2
          else if (sample_inside_mask) then
             invN_rms(i,j,INSIDE_MASK)  = 1.d0 / noisemap(pixels(i),j)**2
          end if
       end do
       if (all(mask(:,j) == 0.d0)) invN_rms(:,j,:) = 0.d0
    end do
    sqrt_invN_rms = sqrt(invN_rms)
    deallocate(noisemap)

    if (trim(noise_format) == 'rms') then
    
       ! OK by default

    else if (trim(noise_format) == 'dense_matrix') then

       call int2string(map_id, map_text)
       call get_parameter(paramfile, 'INV_N_DENSE_ORDERING', par_int=inv_N_order)       

       ! Set up (linear) mask2map array
       n = sum(mask)
       allocate(mask2map(n))
       k = 1
       do j = 1, nmaps
          if (inv_N_order == 2) call convert_ring2nest(nside, mask(:,j))
          do i = 0, npix-1
             if (mask(i,j) == 1.d0) then
                mask2map(k) = (j-1)*npix + i
                k           = k+1
             end if
          end do
          if (inv_N_order == 2) call convert_nest2ring(nside, mask(:,j))
       end do

       ! Find which segment of the matrix this processor should do
       numpix   = sum(mask)
       col_from = numpix / numprocs * myid_alms+1
       col_to   = numpix / numprocs * (myid_alms+1)
       if (myid_alms == numprocs-1) col_to = numpix

       paramtext = 'INV_N_MAT' // map_text
       call get_parameter(paramfile, trim(paramtext), par_string=noisefile)
       allocate(invN_dense(numpix, col_from:col_to))
       open(unit,file=trim(noisefile),form='unformatted')
       read(unit) n
       if (n /= numpix) then
          write(*,*) 'ERROR: Number of pixels in covariance matrix file does not match mask', n, numpix
          call mpi_finalize(ierr)
          stop
       end if
       
       read(unit) i ! Ordering
       read(unit) i ! Polarization status
       write(*,fmt='(a,i8,i8)') 'Size of inverse noise covariance matrix         = ', n, col_to-col_from+1
       do i = 1, col_from-1
          read(unit) 
       end do
       do i = col_from, col_to
          read(unit) invN_dense(:,i)
          call mpi_barrier(mpi_comm_world, ierr)
       end do
       close(unit)

       paramtext = 'SQRT_INV_N_MAT' // map_text
       call get_parameter(paramfile, trim(paramtext), par_string=noisefile)
       allocate(sqrt_invN_dense(numpix, col_from:col_to))
       open(unit,file=trim(noisefile),form='unformatted')
       read(unit) n
       read(unit) i ! Ordering
       read(unit) i ! Polarization status
       write(*,fmt='(a,i8,i8)') 'Size of sqrt of inverse noise covariance matrix = ', n, col_to-col_from+1
       do i = 1, col_from-1
          read(unit) 
       end do
       do i = col_from, col_to
          read(unit) sqrt_invN_dense(:,i)
       end do
       close(unit)


    else 
       write(*,*) 'Error -- noise format ', trim(noise_format), ' not supported'
       call mpi_finalize(ierr)
       stop
    end if

  end subroutine initialize_N_mult_mod


  subroutine clean_up_N_mult_mod
    implicit none

  end subroutine clean_up_N_mult_mod


  subroutine set_mask_state(mask_state_in)
    implicit none

    integer(i4b), intent(in), optional :: mask_state_in

    integer(i4b) :: state, ierr

    if (myid_chain == root) state = mask_state_in
    call mpi_bcast(state, 1, MPI_INTEGER, root, comm_chain, ierr)    

    mask_state = STATE

  end subroutine set_mask_state



  ! Noise multiplication interface routines
  subroutine multiply_by_inv_N_one_map(map, N_format)
    implicit none

    real(dp),     dimension(0:,1:), intent(inout)  :: map
    character(len=*),               intent(in), optional :: N_format

    integer(i4b) :: npix, nmaps
    real(dp), allocatable, dimension(:,:) :: temp_map

    npix  = size(map(:,1))
    nmaps = size(map(0,:))

    allocate(temp_map(0:npix-1,nmaps))
    call multiply_by_noise_matrix(.false., map, temp_map, N_format)
    map = temp_map
    deallocate(temp_map)

  end subroutine multiply_by_inv_N_one_map

  subroutine multiply_by_inv_N_two_maps(map_in, map_out, N_format)
    implicit none

    real(dp),     dimension(0:,1:), intent(in)  :: map_in
    real(dp),     dimension(0:,1:), intent(out) :: map_out
    character(len=*),               intent(in), optional :: N_format

    call multiply_by_noise_matrix(.false., map_in, map_out, N_format)

  end subroutine multiply_by_inv_N_two_maps


  subroutine multiply_by_sqrt_inv_N_one_map(map, N_format)
    implicit none

    real(dp),     dimension(0:,1:), intent(inout)  :: map
    character(len=*),               intent(in), optional :: N_format

    integer(i4b) :: npix, nmaps
    real(dp), allocatable, dimension(:,:) :: temp_map

    npix  = size(map(:,1))
    nmaps = size(map(0,:))

    allocate(temp_map(0:npix-1,nmaps))
    call multiply_by_noise_matrix(.true., map, temp_map)
    map = temp_map
    deallocate(temp_map)

  end subroutine multiply_by_sqrt_inv_N_one_map
  
  subroutine multiply_by_sqrt_inv_N_two_maps(map_in, map_out, N_format)
    implicit none

    real(dp),     dimension(0:,1:), intent(in)  :: map_in
    real(dp),     dimension(0:,1:), intent(out) :: map_out
    character(len=*),               intent(in), optional :: N_format

    call multiply_by_noise_matrix(.true., map_in, map_out)

  end subroutine multiply_by_sqrt_inv_N_two_maps



  ! Noise multiplication main routine
  subroutine multiply_by_noise_matrix(do_sqrt, map_in, map_out, N_format)
    implicit none

    logical(lgt),                   intent(in)  :: do_sqrt 
    real(dp),     dimension(0:,1:), intent(in)  :: map_in
    real(dp),     dimension(0:,1:), intent(out) :: map_out
    character(len=*),               intent(in), optional :: N_format

    integer(i4b)     :: i, j, k, m, n, pix, map_size_in, ierr
    character(len=1) :: uplo
    character(len=128) :: N_format_
    real(dp)         :: t1, t2
    real(dp),     allocatable, dimension(:,:) :: map_glob
    real(dp),     allocatable, dimension(:)   :: map_lin, map_lin2, map_lin_glob

    if (present(N_format)) then
       N_format_ = N_format
    else
       N_format_ = noise_format
    end if

    !write(*,*) myid_chain, ', noiseamp = ', noiseamp(map_id)

    if (trim(N_format_) == 'rms') then
       
       ! Multiply with the inverse noise covariance matrix
       if (do_sqrt) then
          map_out = sqrt_invN_rms(:,:,mask_state)*map_in
       else
          map_out = invN_rms(:,:,mask_state)*map_in
       end if

    else if (trim(N_format_) == 'dense_matrix') then

       allocate(map_glob(0:npix-1,nmaps), map_lin_glob(0:npix*nmaps-1), map_lin2(numpix), map_lin(numpix))
       map_glob = 0.d0; map_glob(pixels,:) = map_in
       if (inv_N_order == 2) then
          do i = 1, nmaps
             call convert_ring2nest(nside, map_glob(:,i))
          end do
       end if
       call mpi_allreduce(reshape(map_glob,shape(map_lin_glob)), map_lin_glob, size(map_lin_glob), &
            & MPI_DOUBLE_PRECISION, MPI_SUM, comm_alms, ierr)
       map_lin2 = map_lin_glob(mask2map)

       ! Multiply with N^-1
       m     = size(invN_dense(:,1)) ! = numpix
       n     = size(invN_dense(1,:))
       if (n > 0 .and. m > 0) then
          if (numprocs == 1) then
             if (do_sqrt) then
                call dsymv('L', n, 1.d0, sqrt_invN_dense, n, map_lin2, 1, 0.d0, map_lin, 1)
             else
                call dsymv('L', n, 1.d0, invN_dense,      n, map_lin2, 1, 0.d0, map_lin, 1)
             end if
             map_lin2 = map_lin
          else
             if (do_sqrt) then
                call dgemv('N',m,n,1.d0,sqrt_invN_dense,m,map_lin2(col_from:col_to),1,0.d0,map_lin,1)
             else
                call dgemv('N',m,n,1.d0,invN_dense,     m,map_lin2(col_from:col_to),1,0.d0,map_lin,1)
             end if
             call mpi_allreduce(map_lin, map_lin2, m, MPI_DOUBLE_PRECISION, MPI_SUM, comm_alms, ierr)
          end if
       end if

       map_lin_glob           = 0.d0
       map_lin_glob(mask2map) = map_lin2
       map_glob = reshape(map_lin_glob,shape(map_glob))
       if (inv_N_order == 2) then
          do i = 1, nmaps
             call convert_nest2ring(nside, map_glob(:,i))
          end do
       end if
       map_out = map_glob(pixels,:)
       deallocate(map_lin, map_lin2, map_glob, map_lin_glob)

    end if

    if (do_sqrt) then
       map_out = map_out / noiseamp(map_id)
    else
       map_out = map_out / noiseamp(map_id)**2
    end if

  end subroutine multiply_by_noise_matrix


  subroutine get_noise_map(mask_state, inv_N)
    implicit none
    
    integer(i4b),                   intent(in)            :: mask_state
    real(dp),     dimension(0:,1:), intent(out), optional :: inv_N
    
    integer(i4b) :: ierr
    real(dp),     allocatable, dimension(:,:) :: map1, map2

    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps))
    map1 = 0.d0; map1(pixels,:) = invN_rms(:,:,mask_state)
    call mpi_reduce(map1, map2, size(map1), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm_alms, ierr)
    if (present(inv_N)) inv_N = map2 / noiseamp(map_id)**2
    deallocate(map1,map2)

  end subroutine get_noise_map

  subroutine get_inv_N_sub_matrix(band, mask_state, pixel, inv_N)
    implicit none
    
    integer(i4b),                intent(in)  :: band, mask_state, pixel
    real(dp),     dimension(1:), intent(out) :: inv_N

    inv_N = invN_rms(pixel,:,mask_state) / noiseamp(band)**2

  end subroutine get_inv_N_sub_matrix

  subroutine set_noiseamp(amps)
    implicit none

    real(dp), dimension(1:), intent(in), optional :: amps

    integer(i4b) :: ierr

    if (myid_chain == root) noiseamp = amps
    call mpi_bcast(noiseamp, numband, MPI_DOUBLE_PRECISION, root, comm_chain, ierr)

  end subroutine set_noiseamp

end module comm_N_mult_mod
