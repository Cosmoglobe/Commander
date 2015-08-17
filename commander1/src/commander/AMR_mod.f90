module AMR_mod
  use healpix_types
  use rngmod
  implicit none


contains

  function sample_AMR(handle, lnL)
    implicit none

    type(planck_rng) :: handle
    real(dp)         :: sample_AMR
    interface
       function lnL(x)
         use healpix_types
         implicit none
         real(dp), intent(in) :: x
         real(dp)             :: lnL
       end function lnL
    end interface

    sample_AMR = lnL(0.d0)

  end function sample_AMR


end module AMR_mod
