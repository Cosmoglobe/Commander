module comm_B_mod
  use comm_param_mod
  use comm_map_mod
  use spline_1D_mod
  implicit none

  private
  public comm_B

  type, abstract :: comm_B
     ! Data variables
     character(len=512)           :: type
     real(dp)                     :: r_max
     class(comm_mapinfo), pointer :: info
     real(dp),          allocatable, dimension(:,:) :: b_l
     type(spline_type), allocatable, dimension(:)   :: b_theta  ! {nmaps}
   contains
     ! Data procedures
     procedure(matmulB),     deferred :: conv
     procedure(matmulInvB),  deferred :: deconv
     procedure                     :: getBTheta
     procedure                     :: initBTheta
  end type comm_B

  abstract interface
     subroutine matmulB(self, trans, map)
       import comm_map, comm_B, dp, i4b, lgt
       implicit none
       class(comm_B),   intent(in)    :: self
       logical(lgt),    intent(in)    :: trans
       class(comm_map), intent(inout) :: map
     end subroutine matmulB

     subroutine matmulInvB(self, trans, map)
       import comm_map, comm_B, dp, i4b, lgt
       implicit none
       class(comm_B),   intent(in)    :: self
       logical(lgt),    intent(in)    :: trans
       class(comm_map), intent(inout) :: map
     end subroutine matmulInvB
  end interface

  ! Local variables
  integer(i4b), parameter :: n_beam = 1000  ! Number of sample point in beam spline

contains

  ! Note: Assume b_l(P) = b_l(T) for now.
  subroutine initBTheta(self, b_l, filename)
    implicit none
    class(comm_B),                      intent(inout) :: self
    real(dp),         dimension(0:,1:), intent(in), optional  :: b_l
    character(len=*),                   intent(in), optional  :: filename

    integer(i4b) :: i, j, l, n, unit
    real(dp)     :: norm
    character(len=512) :: line
    real(dp), allocatable, dimension(:) :: x, y, pl

    if (present(filename)) then
       ! Find number of entries
       unit = getlun()
       n    = 0
       open(unit,file=trim(filename))
       do while (.true.)
          read(unit,'(a)',end=1) line
          line = trim(line)
          if (line(1:1) == '#' .or. trim(line) == '') cycle
          n = n+1
       end do
1      close(unit)

       if (n == 0) call report_error("No valid entries in beam profile file " // trim(filename))

       allocate(x(n), y(n))
       open(unit,file=trim(filename))
       n = 0
       do while (.true.)
          read(unit,'(a)',end=2) line
          line = trim(line)
          if (line(1:1) == '#' .or. trim(line) == '') cycle
          n = n+1
          read(line,*) x(n), y(n)
       end do
2      close(unit)
       x = x * pi/180.d0/60.d0
       write(*,*) n
       
    else if (present(b_l)) then
       
       ! Set up radius grid, and compute unnormalized real-space beam profile
       allocate(x(n_beam), y(n_beam), pl(0:self%info%lmax))
       x = 0.d0
       y = 0.d0
       do i = 1, n_beam
          x(i) = self%r_max / (n_beam-1) * (i-1)
          call comp_normalised_Plm(self%info%lmax, 0, x(i), pl)
          do l = 0, self%info%lmax
             !y(i) = y(i) + (2*l+1) * b_l(l,1) * pl(l) * sqrt((2*l+1)/4.d0*pi)
             y(i) = y(i) + (2*l+1) * &
                  & exp(-0.5d0*l*(l+1)*(30.d0*pi/180./60./sqrt(8.d0*log(2.d0)))**2) * pl(l) * sqrt((2*l+1)/4.d0*pi)
          end do
       end do
       
       ! Normalize to unity 2D integral; flat sky approximation for now
       norm = 0.d0
       do i = 2, n_beam
          norm = norm + 2*pi* 0.5d0*(x(i-1)+x(i)) * &  ! 2*pi*r
               & 0.5d0*(y(i-1)-y(i)) * &        ! b(theta)
               & (x(i)-x(i-1))                  ! dr
       end do
       y = y / norm

       deallocate(pl)
    end if

    ! Spline
    allocate(self%b_theta(self%info%nmaps))
    do i = 1, self%info%nmaps
       call spline(self%b_theta(i), x, y, regular=.true.)
    end do

    deallocate(x,y)
    
  end subroutine initBTheta

  function getBTheta(self, r, pol)
    implicit none
    class(comm_B), intent(in) :: self
    real(dp),      intent(in) :: r
    integer(i4b),  intent(in) :: pol
    real(dp)                  :: getBTheta

    getBTheta = splint(self%b_theta(pol), r)

  end function getBTheta
    
  
end module comm_B_mod
