

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
module comm_crosstalk_mod

implicit none
use comm_status_mod




call update_status(status, "add_crosstalk_start"//ctext)

subroutine add_crosstalk(sd, W, N)

!Adds crosstalk to the TODs for a given scan
!sd: data for given scan
!W: crosstalk matrix
!N: number of detectors 

    type(comm_scandata), intent(inout) :: sd
    !type(comm_scandata) :: sd_wxt

    real(sp), allocatable, dimension(1:,1:)    :: tod_wxt
    type(dp), dimension(1:,1:), intent(in) :: W 
    integer(i4b), intent(in)  :: N
    integer(i4b)        :: i, j, k

    allocate(tod_wxt(sd%ndet, sd%ntod))

    tod_wxt = 0.


    do i = 1, sd%ndet
        tod_wxt(:, i) = tod_wxt(:, i) + sd%tod(:, i)
        if (i .leq. N) then
            do j = 1, N
                tod_wxt(:, i) = tod_wxt(:, i) + W(i, i + j)*sd%tod(:, i + j)
            end do
            if (i .neq. 1) then
                do k = 1, N - i
                    tod_wxt(:, i) = tod_wxt(:, i) + W(i, i - k)*sd%tod(:, i - k)
                end do
            end if
        else if (sd%ndet - i .leq. N) then 
            do j = 1, N  
                tod_wxt(:, i) = tod_wxt(:, i) + W(i, i - j)*sd%tod(:, i - j)
            end do
            if (i .neq. sd%ndet) then
                do k = 1, i - N
                    tod_wxt(:, i) = tod_wxt(:, i) + W(i, i + k)*sd%tod(:, i + k)
                end do
            end if
        else
            do j = 1, N
                tod_wxt(:, i) = tod_wxt(:, i) + W(i, i + j)*sd%tod(:, i + j) + W(i, i - j)*sd%tod(:, i - j)
            end do
        end if
    end do

    sd%tod = tod_wxt

end subroutine

call update_status(status, "add_crosstalk_end"//ctext)


end module comm_crosstalk_mod