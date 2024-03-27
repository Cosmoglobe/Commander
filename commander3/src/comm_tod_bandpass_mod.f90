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
module comm_tod_bandpass_mod
  use comm_tod_mod
  implicit none

contains

  subroutine sample_bp(tod, iter, delta, map_sky, handle, chisq_S)
    implicit none
    class(comm_tod),                          intent(inout)  :: tod
    integer(i4b),                             intent(in)     :: iter
    real(dp),            dimension(0:,1:,1:), intent(inout)  :: delta
    !type(shared_2d_sp),  dimension(0:,1:),    intent(inout)  :: smap_sky
    real(sp),            dimension(1:,1:,0:,1:), intent(inout)  :: map_sky
    type(planck_rng),                         intent(inout)  :: handle
    real(dp),            dimension(1:,1:),    intent(in)     :: chisq_S

    integer(i4b) :: i, k, ierr, ndelta, current
    logical(lgt) :: accept
    real(dp)     :: cp, cc, c0, accept_rate, diff

    call timer%start(TOD_BP, tod%band)
    if (tod%myid == 0) then
       ndelta  = size(chisq_S,2)
       current = 1
       c0      = sum(chisq_S(:,current))
       do k = 2, ndelta
          cp          = sum(chisq_S(:,k))
          cc          = sum(chisq_S(:,current))
          diff        = max(cp-cc,0.d0)
          accept_rate = exp(-0.5d0*diff)  
          if (trim(tod%operation) == 'optimize') then
             accept = cp <= cc
          else
             accept = (rand_uni(handle) < accept_rate)
          end if
          !write(*,*) k, cp, cc, accept
          if (accept) then
             cc = cp
             current = k
          end if
       end do
       do k = 1, ndelta
          write(*,*) '|    BP delta = ', real(delta(0:1,1,k),sp), sum(chisq_S(:,k))
       end do
       write(*,fmt='(a,f16.1,a,l3)') ' |    BP chisq diff = ', cc-c0, ', accept = ', current /= 1
    end if

    ! Broadcast new saved data
    call mpi_bcast(current, 1,  MPI_INTEGER, 0, tod%info%comm, ierr)
    if (current /= 1) then
       ! Set current to proposal
       map_sky(:,:,:,1) = map_sky(:,:,:,current)
       delta(:,:,1) =  delta(:,:,current)
    end if
    call timer%stop(TOD_BP, tod%band)
    
  end subroutine sample_bp


end module comm_tod_bandpass_mod
  
