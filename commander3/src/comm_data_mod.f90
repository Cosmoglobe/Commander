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
  use comm_bp_mod
  use comm_noise_mod
  use comm_beam_mod
  use comm_tod_inst_mod
  use comm_N_mod
  use comm_N_rms_mod
  implicit none

  type comm_data_set
     character(len=512)                  :: label, instlabel, unit, comp_sens, noise_format
     integer(i4b)                        :: period, id_abs
     logical(lgt)                        :: sample_gain
     integer(i4b),       allocatable, dimension(:) :: gain_stat
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
     class(comm_tod),     pointer :: tod       => null()
     class(comm_N),       pointer :: N         => null()
     class(B_ptr),         allocatable, dimension(:) :: B
     class(comm_bp_ptr),   allocatable, dimension(:) :: bp
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_smooth
     type(comm_B_bl_ptr),  allocatable, dimension(:) :: B_postproc
     class(comm_N_ptr),     allocatable, dimension(:) :: N_smooth
   contains
     procedure :: RJ2data
     procedure :: chisq => get_chisq
  end type comm_data_set

  integer(i4b) :: numband
  type(comm_data_set), allocatable, dimension(:) :: data
  integer(i4b),        allocatable, dimension(:) :: ind_ds

  
contains

  subroutine initialize_data_mod(cpar, handle)
    !
    ! Routine to initialise Commander3 data
    !
    implicit none
    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b)       :: i, j, k, n, nmaps, numband_tot, ierr

    integer(i4b)                        :: n_dummy

    character(len=1) :: j_str

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

       data(n)%N       => comm_N_rms(cpar, data(n)%rmsinfo, n, i, 0, data(n)%mask, handle)!, regnoise)

    end do


  end subroutine initialize_data_mod

  function get_chisq(self)
    implicit none
    class(comm_data_set), intent(in)           :: self
    real(dp)                                   :: get_chisq

    integer(i4b) :: ierr
    real(dp)     :: chisq
    class(comm_map), pointer :: invN_res => null()

    get_chisq = 0d0    

  end function get_chisq

  function RJ2data(self, det)
    implicit none

    class(comm_data_set), intent(in)           :: self
    integer(i4b),         intent(in), optional :: det
    real(dp)                                   :: RJ2data

    integer(i4b) :: d

    d = 0; if (present(det)) d = det

    RJ2data = 1.d0
    
  end function RJ2data

  subroutine dump_unit_conversion(dir)
    implicit none
    character(len=*), intent(in) :: dir
    integer(i4b) :: i, q, unit

  end subroutine dump_unit_conversion


  subroutine smooth_map(info, alms_in, bl_in, map_in, bl_out, map_out, spinzero)
    implicit none
    class(comm_mapinfo),                      intent(in),   target :: info
    logical(lgt),                             intent(in)           :: alms_in
    real(dp),            dimension(0:,1:),    intent(in)           :: bl_in, bl_out
    class(comm_map),                          intent(inout)        :: map_in
    class(comm_map),                          intent(out), pointer :: map_out
    logical(lgt),                             intent(in), optional :: spinzero

    integer(i4b) :: i, j, l, b, lmax
    logical(lgt) :: spinzero_


  end subroutine smooth_map

end module comm_data_mod
