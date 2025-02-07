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
module comm_cr_mod
  use comm_output_mod
  implicit none

!  private
!  public solve_cr_eqn_by_CG, cr_amp2x, cr_x2amp, cr_computeRHS, cr_matmulA, cr_invM

  interface cr_amp2x
     module procedure cr_amp2x_full
  end interface cr_amp2x

  interface cr_x2amp
     module procedure cr_x2amp_full
  end interface cr_x2amp

contains

  subroutine solve_cr_eqn_by_CG(cpar, samp_group, x, b, stat)
    implicit none
    type(comm_params),                intent(in)    :: cpar
    integer(i4b),                     intent(in)    :: samp_group
    real(dp),          dimension(1:), intent(out)   :: x
    real(dp),          dimension(1:), intent(in)    :: b
    integer(i4b),                     intent(out)   :: stat

    x = 0d0
    stat = 0

    
  end subroutine solve_cr_eqn_by_CG

  subroutine cr_compute_chisq(comm, samp_group, x, chisq)
    implicit none
    integer(i4b),            intent(in)  :: comm, samp_group
    real(dp), dimension(1:), intent(in)  :: x
    real(dp),                intent(out) :: chisq

    integer(i4b) :: j, k, n
    real(dp), allocatable, dimension(:)   :: x_out
    real(dp), allocatable, dimension(:,:) :: alm, pamp
    class(comm_comp),   pointer :: c => null()

      chisq = 0d0


  end subroutine cr_compute_chisq

  subroutine cr_amp2x_full(samp_group, x) 
    implicit none
    integer(i4b),           intent(in)  :: samp_group
    real(dp), dimension(:), intent(out) :: x

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    x = 0d0


  end subroutine cr_amp2x_full

  subroutine cr_x2amp_full(samp_group, x)
    implicit none
    integer(i4b),           intent(in) :: samp_group
    real(dp), dimension(:), intent(in) :: x

    integer(i4b) :: i, ind
    class(comm_comp), pointer :: c => null()

    
  end subroutine cr_x2amp_full

  ! ---------------------------
  ! Definition of linear system
  ! ---------------------------

  subroutine cr_computeRHS(operation, resamp_cmb, only_pol, handle, handle_noise, mask, samp_group, rhs)
    implicit none
    character(len=*),                            intent(in)             :: operation
    logical(lgt),                                intent(in)             :: resamp_cmb, only_pol
    type(planck_rng),                            intent(inout)          :: handle, handle_noise
    integer(i4b),                                intent(in)             :: samp_group
    real(dp),         allocatable, dimension(:), intent(in)             :: mask
    real(dp),         allocatable, dimension(:), intent(out)            :: rhs

    integer(i4b) :: i, j, l, m, k, n, ierr
    real(dp)     :: tmp
    class(comm_map),     pointer                 :: map  => null()
    class(comm_map),     pointer                 :: Tm   => null()
    class(comm_map),     pointer                 :: mu   => null()
    class(comm_comp),    pointer                 :: c    => null()
    class(comm_mapinfo), pointer                 :: info => null()
    real(dp),        allocatable, dimension(:,:) :: eta, Tp

    ! Initialize output vector
    allocate(rhs(ncr))
    rhs = 0.d0


  end subroutine cr_computeRHS

  function cr_matmulA(x, samp_group)
    implicit none

    real(dp),     dimension(1:),     intent(in)  :: x
    integer(i4b),                    intent(in)  :: samp_group
    real(dp),     dimension(size(x))             :: cr_matmulA

    real(dp)                  :: t1, t2
    integer(i4b)              :: i, j, l, myid, lmax
    class(comm_map),  pointer :: map => null(), pmap => null(), map_buff => null()
    class(comm_comp), pointer :: c => null()
    real(dp),        allocatable, dimension(:)   :: y, sqrtS_x
    real(dp),        allocatable, dimension(:,:) :: alm, m, pamp
    class(comm_mapinfo), pointer :: info => null()


    cr_matmulA = 0d0

  end function cr_matmulA

  function cr_invM(comm, x, samp_group)
    implicit none
    integer(i4b),                        intent(in) :: comm, samp_group
    real(dp),              dimension(:), intent(in) :: x
    real(dp), allocatable, dimension(:)             :: cr_invM

    integer(i4b) :: ierr
    logical(lgt) :: Q_is_active
    real(dp), allocatable, dimension(:,:) :: alm, alm0
    real(dp), allocatable, dimension(:)   :: Qx
    class(comm_comp), pointer :: c => null()

    if (.not. allocated(cr_invM)) allocate(cr_invM(size(x)))
  end function cr_invM


  subroutine applyDeflatePrecond(x, Qx)
    real(dp), dimension(:), intent(in)  :: x
    real(dp), dimension(:), intent(out) :: Qx

    integer(i4b) :: ierr
    real(dp), allocatable, dimension(:,:) :: alm, Qalm
    class(comm_comp), pointer :: c => null()

    Qx = 0.d0

    return

    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_def > -1) then
             allocate(alm(0:c%x%info%nalm-1,c%x%info%nmaps))
             allocate(Qalm(0:c%x%info%nalm-1,c%x%info%nmaps))
             call cr_extract_comp(c%id, x, alm)
             call c%applyDeflatePrecond(alm, Qalm)
             call cr_insert_comp(c%id, .false., Qalm, Qx)
             deallocate(alm, Qalm)
          end if
       end select
       c => c%nextComp()
    end do

  end subroutine applyDeflatePrecond


  subroutine update_precond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update
    class(comm_comp), pointer :: c => null()
    logical(lgt), save :: first_call = .true.

    ! Set up deflation preconditioner for CMB+diagonal only
    if (.false.) then
       c   => compList
       do while (associated(c))
          select type (c)
          class is (comm_diffuse_comp)
             call c%updateDeflatePrecond()
          end select
          c => c%nextComp()
       end do
    end if

    call updateDiffPrecond(samp_group, force_update)
    call updatePtsrcPrecond(samp_group)
    call updateTemplatePrecond(samp_group)

    !if (.not. first_call) return

    ! Set up low-l preconditioner for CMB+diagonal only
    c   => compList
    do while (associated(c))
       select type (c)
       class is (comm_diffuse_comp)
          if (c%lmax_pre_lowl > -1) then
             call c%updateLowlPrecond()
          end if
       end select
       c => c%nextComp()
    end do

    first_call = .false.

  end subroutine update_precond
  
end module comm_cr_mod
