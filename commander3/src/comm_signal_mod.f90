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
module comm_signal_mod
  use comm_chisq_mod
  use comm_md_comp_mod
  use comm_line_comp_mod
  implicit none

  integer(i4b), parameter, private :: MAXSAMPGROUP     = 100


contains

  subroutine initialize_signal_mod(cpar)
    implicit none
    type(comm_params), intent(in) :: cpar

    integer(i4b) :: i, n
    class(comm_comp), pointer :: c => null()
    class(comm_comp), pointer :: c_two => null()
    logical(lgt) :: prior_exists


  end subroutine initialize_signal_mod

  subroutine dump_components(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(i4b) :: i, unit
    class(comm_comp), pointer :: c => null()

    
  end subroutine dump_components

  subroutine sample_amps_by_CG(cpar, samp_group, handle, handle_noise)
    implicit none

    type(comm_params), intent(in)    :: cpar
    integer(i4b),      intent(in)    :: samp_group
    type(planck_rng),  intent(inout) :: handle, handle_noise

    integer(i4b) :: stat, i, j, l, m, nactive
    real(dp)     :: Nscale = 1.d-4
    class(comm_comp), pointer :: c => null()
    character(len=32) :: cr_active_bands(100)
    real(dp),           allocatable, dimension(:) :: rhs, x, mask
    class(comm_map),     pointer :: res  => null()



  end subroutine sample_amps_by_CG


  subroutine sample_all_amps_by_CG(cpar, handle, handle_noise, cg_groups)
    !
    !
    !  Convenience function for performing amplitude sampling over
    !  all sampling groups
    !
    !
    implicit none

    type(comm_params), intent(in)            :: cpar
    type(planck_rng),  intent(inout)         :: handle, handle_noise
    character(len=512), intent(in), optional :: cg_groups


    integer(i4b)                          :: samp_group, i, n
    integer(i4b), dimension(MAXSAMPGROUP) :: group_inds
    character(len=3) :: toks(MAXSAMPGROUP)


    call timer%stop(TOT_AMPSAMP)

  end subroutine sample_all_amps_by_CG

  subroutine initPrecond(comm, samp_group)
    implicit none
    integer(i4b), intent(in) :: comm, samp_group
  end subroutine initPrecond

  subroutine add_to_complist(c)
    implicit none
    class(comm_comp), pointer :: c 
    
  end subroutine add_to_complist

  subroutine initialize_from_chain(cpar, handle, init_samp, init_from_output, first_call)
    implicit none
    type(comm_params), intent(in)           :: cpar
    type(planck_rng),  intent(inout)        :: handle
    integer(i4b),      intent(in), optional :: init_samp
    logical(lgt),      intent(in), optional :: init_from_output, first_call

    integer(i4b)              :: i, j, ext(2), initsamp, initsamp2
    character(len=4)          :: ctext
    character(len=6)          :: itext, itext2
    character(len=512)        :: chainfile, hdfpath
    class(comm_comp), pointer :: c => null()
    type(hdf_file) :: file, file2
    class(comm_N),      pointer :: N => null() 
    class(comm_map),    pointer :: rms => null()
    real(dp), allocatable, dimension(:,:) :: bp_delta, regnoise

    
    
  end subroutine initialize_from_chain


  subroutine sample_powspec(handle, ok)
    implicit none
    type(planck_rng),  intent(inout) :: handle
    logical(lgt),      intent(out)   :: ok

    class(comm_comp), pointer :: c => null()

      ok = .false.

  end subroutine sample_powspec


  subroutine sample_partialsky_tempamps(cpar, handle)
    implicit none

    type(comm_params), intent(in)    :: cpar
    type(planck_rng),  intent(inout) :: handle

    integer(i4b) :: i, ierr
    logical(lgt) :: skip
    real(dp)     :: vals(2), vals2(2), mu, sigma, amp, mu_p, sigma_p
    class(comm_map),           pointer :: res => null(), invN_T => null()
    class(comm_comp),          pointer :: c => null()
    class(comm_template_comp), pointer :: pt => null()
    

  end subroutine sample_partialsky_tempamps


  subroutine synchronize_bp_delta(initHDF)
    implicit none
    logical(lgt), intent(in) :: initHDF

    integer(i4b) :: i, j, ndet
    real(dp)     :: mu


  end subroutine synchronize_bp_delta
  

  subroutine sample_joint_alm_Cl(handle)
    implicit none
    type(planck_rng), intent(inout) :: handle

  end subroutine sample_joint_alm_Cl


end module comm_signal_mod
