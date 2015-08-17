program amr
  use healpix_types
  use rngmod
  use AMR_mod
  implicit none

  integer(i4b)        :: iargc, ierr, numprocs, myid, unit, root, num_chain, num_chain_per_realization
  integer(i4b)        :: i, j, k, l, m, s, lmax, nside, npix, nmaps, nspec
  integer(i4b)        :: base_seed
  type(planck_rng)    :: rng_handle

  real(dp), allocatable, dimension(:) :: samples

  j = 10000
  allocate(samples(j))

  call rand_init(rng_handle, 387491)
  open(58,file='samp.dat')
  do i = 1, j
     samples(i) = sample_AMR(rng_handle, lnL_gaussian)
     write(58,*) samples(i)
  end do
  write(*,*) mean(samples), sqrt(variance(samples))
  close(58)
  stop

contains


  function lnL_gaussian(x)
    use healpix_types
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: lnL_gaussian

    lnL_gaussian = (x-3.d0)**2 / 0.4d0**2
  end function lnL_gaussian

end program amr
