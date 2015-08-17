module comm_Cl_sampling_mod
  use comm_mp_mod
  use comm_Cl_util_mod
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

contains

  subroutine initialize_Cl_sampling_mod(paramfile)
    implicit none

    character(len=128),                   intent(in) :: paramfile

  end subroutine initialize_Cl_sampling_mod

  subroutine cleanup_Cl_sampling_mod
    implicit none

  end subroutine cleanup_Cl_sampling_mod
  
  

  ! Routine for sampling a power spectrum given a signal map
  subroutine sample_inv_wishart_spectrum(alms, rng_handle, cls)
    implicit none

    type(planck_rng),           intent(inout) :: rng_handle
    real(dp), dimension(1:,1:), intent(in)    :: alms
    real(dp), dimension(0:,1:), intent(inout) :: cls

    integer(i4b) :: bin, b, i, j, k, l, m, n, p, ind, b1, b2, col
    logical(lgt), allocatable, dimension(:,:) :: pattern
    integer(i4b), allocatable, dimension(:) :: i2p
    real(dp), allocatable, dimension(:,:)   :: y, y_t
    real(dp), allocatable, dimension(:,:)   :: C_b
    real(dp), allocatable, dimension(:,:)   :: sigma, s

    allocate(sigma(nmaps,nmaps), pattern(nmaps,nmaps))

    do bin = 1, num_Cl_bin
       if (.not. any(Cl_bin_stat(bin,:) == 'S') .and. .not. any(Cl_bin_stat(bin,:) == 'M')) cycle
       
       sigma = 0.d0
       n     = 0.d0
       do l = Cl_bins(bin,1), Cl_bins(bin,2)
          do m = -l, l
             ind = lm2ind(l,m)
             n   = n + 1
             do i = 1, nmaps
                do j = 1, nmaps
                   sigma(i,j) = sigma(i,j) + alms(ind,i) * alms(ind,j) * real(l*(l+1),dp) / (2.d0*pi)
                end do
             end do
          end do
       end do

       ! Set up sampling pattern
       k = 1
       do i = 1, nmaps
          do j = i, nmaps
             pattern(i,j) = Cl_bin_stat(bin,k) == 'S' .or. Cl_bin_stat(bin,k) == 'M'
             pattern(j,i) = pattern(i,j)
             k            = k+1
          end do
       end do

       ! Keep sampling until all elements have been treated
       do while (any(pattern))
          ! Find which elements to include this time
          do col = 1, nmaps
             if (any(pattern(:,col))) exit
          end do

          ! Extract the appropriate segment
          p = count(pattern(:,col))
          allocate(s(p,p), y(p,1), y_t(1,p), i2p(p), C_b(p,p))
          j = 1
          do i = 1, nmaps
             if (pattern(i,col)) then
                i2p(j) = i
                j      = j+1
             end if
          end do
          s = sigma(i2p,i2p)
          call invert_matrix(s)
          call cholesky_decompose_single(s)

          ! Draw sample
          C_b = 0.d0
          do i = 1, n - p - 1
             do j = 1, p
                y(j,1) = rand_gauss(rng_handle)
             end do
             y_t(1,:) = matmul(s, y(:,1))
             y(:,1)   = y_t(1,:)
             C_b      = C_b + matmul(y(:,1:1), y_t(1:1,:))
          end do
          call invert_matrix(C_b)

          ! Copy information over to output cls array
          do i = 1, p
             do j = i, p
                ind = i2p(i)*(1-i2p(i))/2 + (i2p(i)-1)*nmaps + i2p(j)
                cls(Cl_bins(bin,1):Cl_bins(bin,2),ind) = C_b(i,j) !/ real(l*(l+1),dp) * 2.d0*pi
             end do
          end do
          
          ! Remove current elements from pattern, and prepare for next round
          pattern(i2p,i2p) = .false.
          deallocate(s, y, y_t, i2p, C_b)

       end do
       
    end do

    deallocate(sigma, pattern)

  end subroutine sample_inv_wishart_spectrum


end module comm_Cl_sampling_mod
