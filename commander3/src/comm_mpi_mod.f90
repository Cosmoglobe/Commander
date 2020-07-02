module comm_mpi_mod
  use healpix_types
  implicit none

  include "mpif.h"

  integer(i4b), parameter :: MPI_BUF_SIZE = 10000

contains

  
  subroutine mpi_allreduce_buffer(oper, dim, comm, ierr, in_place, &
       & in_dp_1d, out_dp_1d, in_dp_2d, out_dp_2d, in_dp_3d, out_dp_3d)
    implicit none
    integer,                        intent(in)              :: oper
    integer(i4b),                   intent(in)              :: dim, comm
    integer(i4b),                   intent(inout)           :: ierr
    logical(lgt),                   intent(in)              :: in_place
    real(dp),     dimension(:),     intent(in),    optional :: in_dp_1d
    real(dp),     dimension(:),     intent(inout), optional :: out_dp_1d
    real(dp),     dimension(:,:),   intent(in),    optional :: in_dp_2d
    real(dp),     dimension(:,:),   intent(inout), optional :: out_dp_2d
    real(dp),     dimension(:,:,:), intent(in),    optional :: in_dp_3d
    real(dp),     dimension(:,:,:), intent(inout), optional :: out_dp_3d

    integer(i4b) :: i, j, n, m

    if (present(out_dp_1d)) then
       n = size(out_dp_1d,dim=1)
       i = 0
       do while (i <= n)
          j = min(i+MPI_BUF_SIZE-1,n)
          if (in_place) then
             call mpi_allreduce(MPI_IN_PLACE, out_dp_1d(i:j), j-i+1, &
                  & MPI_DOUBLE_PRECISION, oper, comm, ierr)
          else
             call mpi_allreduce(in_dp_1d(i:j), out_dp_1d(i:j), j-i+1, &
                  & MPI_DOUBLE_PRECISION, oper, comm, ierr)
          end if
          i = i+MPI_BUF_SIZE
       end do
    end if

    if (present(out_dp_2d)) then
       n = size(out_dp_2d,dim=dim)
       m = size(out_dp_2d)/n
       i = 0
       do while (i <= n)
          j = min(i+MPI_BUF_SIZE-1,n)
          if (in_place) then
             if (dim == 1) then
                call mpi_allreduce(MPI_IN_PLACE, out_dp_2d(i:j,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)

             else
                call mpi_allreduce(MPI_IN_PLACE, out_dp_2d(:,i:j), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             end if
          else
             if (dim == 1) then
                call mpi_allreduce(in_dp_2d(i:j,:), out_dp_2d(i:j,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             else
                call mpi_allreduce(in_dp_2d(:,i:j), out_dp_2d(:,i:j), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             end if
          end if
          i = i+MPI_BUF_SIZE
       end do
    end if

    if (present(out_dp_3d)) then
       n = size(out_dp_3d,dim=dim)
       m = size(out_dp_3d)/n
       i = 0
       do while (i <= n)
          j = min(i+MPI_BUF_SIZE-1,n)
          if (in_place) then
             if (dim == 1) then
                call mpi_allreduce(MPI_IN_PLACE, out_dp_3d(i:j,:,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             else if (dim == 2) then
                call mpi_allreduce(MPI_IN_PLACE, out_dp_3d(:,i:j,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             else
                call mpi_allreduce(MPI_IN_PLACE, out_dp_3d(:,:,i:j), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             end if
          else
             if (dim == 1) then
                call mpi_allreduce(in_dp_3d(i:j,:,:), out_dp_3d(i:j,:,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             else if (dim == 2) then
                call mpi_allreduce(in_dp_3d(:,i:j,:), out_dp_3d(:,i:j,:), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)

             else
                call mpi_allreduce(in_dp_3d(:,:,i:j), out_dp_3d(:,:,i:j), (j-i+1)*m, &
                     & MPI_DOUBLE_PRECISION, oper, comm, ierr)
             end if
          end if
          i = i+MPI_BUF_SIZE
       end do
    end if

  end subroutine mpi_allreduce_buffer



end module comm_mpi_mod
