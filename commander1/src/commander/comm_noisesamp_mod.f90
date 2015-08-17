module comm_noisesamp_mod
  use comm_N_mult_mod
  use ARS_mod
  use comm_mp_mod
  implicit none

  integer(i4b),                          private :: band
  real(dp),                              private :: my_chisq, nu
  real(dp), allocatable, dimension(:,:)          :: N_prior

contains

  subroutine initialize_noisesamp_mod(paramfile)
    implicit none

    character(len=*), intent(in) :: paramfile

    integer(i4b)       :: i, j
    character(len=2)   :: i_text
    character(len=256) :: paramtext

    ! Read noise amplitude priors
    allocate(N_prior(numband,2))
    do i = 1, numband
       call int2string(i, i_text)
       paramtext = 'NOISEAMP_PRIOR_MEAN' // i_text
       call get_parameter(paramfile, trim(paramtext), par_dp=N_prior(i,1))
       paramtext = 'NOISEAMP_PRIOR_RMS' // i_text
       call get_parameter(paramfile, trim(paramtext), par_dp=N_prior(i,2))
    end do

  end subroutine initialize_noisesamp_mod


  subroutine sample_noiseamps(handle, s, noiseamp)
    implicit none

    type(planck_rng),                intent(inout) :: handle
    type(genvec),                    intent(in)    :: s
    real(dp),         dimension(1:), intent(inout) :: noiseamp

    integer(i4b) :: i, j, status
    real(dp)     :: alpha, beta, k, x, x_init(3), mu, sigma
    real(dp), allocatable, dimension(:) :: chisq_band, amps, nu_band

    allocate(chisq_band(numband), amps(numband), nu_band(numband))
    call compute_chisq(map_id, s, chisq_band=chisq_band, nu_band=nu_band)

    do i = 1, numband
       if (N_prior(i,2) == 0.d0) cycle
       band      = i
       my_chisq  = chisq_band(i) * noiseamp(i)**2
       nu        = nu_band(i)
       status    = 0
       mu        = sqrt(my_chisq / nu)
       sigma     = (2.d0 * my_chisq**2 / nu**3)**0.25d0
       x_init(2) = mu
       x_init(1) = max(mu - sigma, 1.d-6)
       x_init(3) =     mu + sigma

       status = 0
       noiseamp(i) = sample_InvSamp(handle, x_init, lnL_noiseamp, status=status)
       if (status /= 0) then
          write(*,*) 'comm_noiseamp_mod -- Error: noise amplitude sampler failed, status = ', status
          write(*,*) 'comm_noiseamp_mod --        Prior too tight?'
          write(*,*) 'chisq = ', my_chisq, nu
          open(58,file='P_noise.dat')
          do j = 2, 10000
             x = 0 + (1000-0)/(10000.d0-1.d0) * (j-1)
             write(58,*) x, lnL_noiseamp(x)
          end do
          close(58)
          stop
       end if
    end do
    deallocate(chisq_band, amps, nu_band)

    !stop

  end subroutine sample_noiseamps
  
  function lnL_noiseamp(x)
    use healpix_types
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: lnL_noiseamp

    if (x < 0.d0) then
       lnL_noiseamp = -1.d30
    else
       ! Posterior is a product of an inverse gamma likelihood and a Gaussian prior
       lnL_noiseamp = -0.5d0*((my_chisq/x**2) + 2.d0*nu*log(x))
       if (N_prior(band,2) > 0.d0) then
          lnL_noiseamp = lnL_noiseamp -0.5d0 * ((x-N_prior(band,1))**2 / N_prior(band,2)**2)
       end if
    end if

  end function lnL_noiseamp
  
end module comm_noisesamp_mod
