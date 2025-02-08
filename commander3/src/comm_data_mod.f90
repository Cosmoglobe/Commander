!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================
module comm_data_mod
  use comm_N_mod
  use comm_N_rms_mod
  implicit none

  type comm_data_set
     character(len=512)                  :: label, instlabel, unit, comp_sens, noise_format
     integer(i4b)                        :: id_abs
     real(dp)                            :: gain, gain_tmp, gain_prior(2)
     real(dp), allocatable, dimension(:) :: gain_sigmas
     character(len=128)                  :: gain_comp
     integer(i4b)                        :: gain_lmin, gain_lmax
     integer(i4b)                        :: ndet
     character(len=128)                  :: tod_type
     integer(i4b)                        :: tod_freq
     logical(lgt)                        :: pol_only, subtract_zodi
     logical(lgt)                        :: cr_active

     class(comm_mapinfo), pointer :: info      => null()
     class(comm_mapinfo), pointer :: rmsinfo   => null()
     class(comm_map),     pointer :: map       => null()
     !class(comm_map),     pointer :: map0      => null() !for TOD data if outputing to HDF
     class(comm_map),     pointer :: res       => null()
     class(comm_map),     pointer :: mask      => null()
     class(comm_map),     pointer :: procmask  => null()
     class(comm_map),     pointer :: gainmask  => null()
     class(comm_N),       pointer :: N         => null()
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data

  
contains

  subroutine initialize_data_mod(cpar)
    !
    ! Routine to initialise Commander3 data
    !
    implicit none
    type(comm_params), intent(in)    :: cpar

    integer(i4b)       :: i, n, nmaps, numband_tot

    ! Read all data sets
    numband = count(cpar%ds_active)
    numband_tot = cpar%numband
    allocate(data(numband))
    n = 0
    do i = 1, numband_tot
       if (.not. cpar%ds_active(i)) cycle
       n                      = n+1

       write(*,fmt='(a,i5,a,a)') ' |  Reading data set ', i, ' : ', trim(data(n)%label)
       nmaps = 1
       write(*,*) cpar%ds_nside(i), cpar%ds_lmax(i), cpar%ds_polarization(i), nmaps
       data(n)%info => comm_mapinfo(cpar%comm_chain, cpar%ds_nside(i), cpar%ds_lmax(i), &
            & nmaps, cpar%ds_polarization(i))
       data(n)%rmsinfo => data(n)%info

       ! Initialize mask structures
       data(n)%mask     => comm_map(data(n)%info)

       data(n)%N       => comm_N_rms(cpar, data(n)%rmsinfo, i)!, regnoise)

    end do


  end subroutine initialize_data_mod





end module comm_data_mod
