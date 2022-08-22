

!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the comm_crosstalk_mod.f90License, or
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
module comm_crosstalk_mod

use comm_status_mod
use comm_tod_driver_mod

implicit none

contains

    subroutine add_crosstalk(sd, xt_mat, N)

    !Adds crosstalk to the TODs for a given scan
    !sd: data for given scan
    !W: crosstalk matrix
    !N: number of detectors we consider a crosstalk contribution from,
    !   i.e. for a detector i we will consider crosstalk from detectors i - N to i + N

        implicit none
        type(comm_scandata), intent(inout) :: sd
        !type(comm_scandata) :: sd_wxt

        real(dp), allocatable, dimension(:,:)    :: tod_wxt
        real(dp), dimension(sd%ndet,sd%ndet), intent(in) :: xt_mat
        integer(i4b), intent(in)  :: N
        integer(i4b)        :: i, j, k, l


        allocate(tod_wxt(1:sd%ntod, 1:sd%ndet))


        tod_wxt = 0.

        !write(*,*) 'before', 'tod', shape(tod_wxt), '1 tod', shape(tod_wxt(:, 1)), 'xt mat', shape(xt_mat), '1 el xt mat', shape(xt_mat(1, 2)), 'xt times 1 tod', shape(xt_mat(1, 2)*sd%tod(:,1)), 'val xt mat', xt_mat, 'val 1 el xt mat',  xt_mat(0,1)

        do i = 1, sd%ndet
            tod_wxt(:, i) = tod_wxt(:, i) + sd%tod(:, i)
            if (i .le. N) then
                do j = 1, N
                    if (i + j .le. sd%ndet) then
                        tod_wxt(:, i) = tod_wxt(:, i) + xt_mat(i, i + j)*sd%tod(:, i + j)

                    end if
                end do
                if (i .ne. 1) then
                    do k = 1, N - i
                        if (i - k .ge. 1) then
                            tod_wxt(:, i) = tod_wxt(:, i) + xt_mat(i, i - k)*sd%tod(:, i - k)
                        end if
                    end do
                end if
            else if (sd%ndet - i .le. N) then 
                do j = 1, N
                    if (i - j .ge. 1) then
                        tod_wxt(:, i) = tod_wxt(:, i) + xt_mat(i, i - j)*sd%tod(:, i - j)
                    end if
                end do
                if (i .ne. sd%ndet) then
                    do k = 1, i - N
                        if (i + k .le. sd%ndet) then
                            tod_wxt(:, i) = tod_wxt(:, i) + xt_mat(i, i + k)*sd%tod(:, i + k)
                        end if
                    end do
                end if
            else
                do j = 1, N
                    tod_wxt(:, i) = tod_wxt(:, i) + xt_mat(i, i + j)*sd%tod(:, i + j) + xt_mat(i, i - j)*sd%tod(:, i - j)
                end do
            end if
        end do

        sd%tod = tod_wxt

        !write(*,*) 'after', 'shape xt tod', shape(tod_wxt), 'shape sd tod', shape(sd%tod)


    end subroutine add_crosstalk


    subroutine estimate_crosstalk(sd, N)

    end subroutine estimate_crosstalk


end module comm_crosstalk_mod