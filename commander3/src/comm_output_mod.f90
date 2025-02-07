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
module comm_output_mod
  use comm_chisq_mod
  implicit none

contains

  subroutine init_chain_file(cpar, iter)
    implicit none
    
    type(comm_params), intent(in)  :: cpar
    integer(i4b),      intent(out) :: iter

    integer(i4b)                 :: i, j, hdferr, ierr
    logical(lgt)                 :: exist
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix, chainfile, hdfpath, fg_file
    type(hdf_file)   :: file
    TYPE(h5o_info_t) :: object_info

    call int2string(cpar%mychain, ctext)
    chainfile = trim(adjustl(cpar%outdir)) // '/chain' // &
         & '_c' // trim(adjustl(ctext)) // '.h5'

    ! Delete existing chain file if necessary; create new file if necessary; open file
    iter = 1
    if (cpar%myid_chain == 0) then
       inquire(file=trim(chainfile), exist=exist)
       if (trim(cpar%chain_status) == 'new' .or. .not. exist) then
          call open_hdf_file(chainfile, file, 'w')

          ! testing hdf parameter output
          call close_hdf_file(file)
          iter = -1
       else if (trim(cpar%chain_status) == 'append') then
          call open_hdf_file(chainfile, file, 'r')
          exist = .true.
          do while (exist)
             call int2string(iter, itext)
             call h5eset_auto_f(0, hdferr)
             call h5oget_info_by_name_f(file%filehandle, itext, object_info, hdferr)
             exist = (hdferr == 0)
             if (exist) iter = iter+1
          end do
          iter = max(1,iter-1)
          write(*,*) '|  Continuing chain '//ctext// ' on iteration ', iter
          call close_hdf_file(file)          
       else
          write(*,*) 'Unsupported chain mode =', trim(cpar%chain_status)
          call mpi_finalize(ierr)
          stop
       end if
       
       !delete fg_ind_mean_cXXXX.dat if it exists
       fg_file=trim(cpar%outdir)//'/fg_ind_mean_c' // trim(adjustl(ctext))//'.dat'
       inquire(file=fg_file, exist=exist)
    end if
    call mpi_bcast(iter, 1, MPI_INTEGER, 0, cpar%comm_chain, ierr)

  end subroutine init_chain_file

  subroutine output_FITS_sample(cpar, iter, output_hdf)
    implicit none
    
    type(comm_params), intent(in) :: cpar
    integer(i4b),      intent(in) :: iter
    logical(lgt),      intent(in) :: output_hdf

    integer(i4b)                 :: i, j, p, hdferr, ierr, unit, p_min, p_max
    real(dp)                     :: chisq, chisq_eff, t1, t2, t3, t4, theta_sum, uscale
    logical(lgt)                 :: exist, init, new_header
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix, chainfile, hdfpath, fg_file, temptxt, fg_txt
    character(len=512)           :: label, label1, label2, label3
    character(len=2048)          :: outline, fg_header
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map),     pointer :: map => null(), chisq_map => null(), chisq_sub => null()
    class(comm_map),     pointer :: rms_map => null(), chisq_map_eff => null()
    class(comm_comp),    pointer :: c => null()
    class(comm_N),      pointer :: N => null()
    type(hdf_file) :: file
    TYPE(h5o_info_t) :: object_info

  end subroutine output_FITS_sample

  subroutine output_inst_params(cpar, chainfile, iter, filename, output_hdf)
    implicit none
    type(comm_params), intent(in) :: cpar
    type(hdf_file),    intent(in) :: chainfile
    character(len=*),  intent(in) :: iter
    character(len=*),  intent(in) :: filename
    logical(lgt),      intent(in) :: output_hdf
    
    integer(i4b) :: unit, i, j
    real(dp), allocatable, dimension(:,:) :: bp_delta

    if (cpar%myid_chain /= 0) return

    if (output_hdf) then
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/gain')
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/bandpass')
    end if
    

  end subroutine output_inst_params

  ! =============== testing hdf parameter output ===============
  ! ============================================================
end module comm_output_mod
