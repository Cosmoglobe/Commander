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
module comm_B_FIRAS_mod
  use comm_B_mod
  implicit none

  private
  public comm_B_FIRAS, comm_B_FIRAS_ptr

  type, extends (comm_B) :: comm_B_FIRAS
     integer(i4b) :: nside_hires, npix_hires
     real(dp)     :: M_ecl(3,3)
     real(dp), allocatable, dimension(:,:) :: vecs
   contains
     ! Data procedures
     procedure :: conv           => matmulB_firas
     procedure :: deconv         => matmulInvB_firas
     procedure :: update         => updateBeam
  end type comm_B_FIRAS

  interface comm_B_FIRAS
     procedure constructor
  end interface comm_B_FIRAS

  type comm_B_FIRAS_ptr
     type(comm_B_FIRAS), pointer :: p => null()
  end type comm_B_FIRAS_ptr

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, info, id, id_abs, fwhm, mb_eff, nside, pixwin, init_realspace)
    implicit none
    type(comm_params),                  intent(in)           :: cpar
    type(comm_mapinfo), target,         intent(in)           :: info
    integer(i4b),                       intent(in)           :: id, id_abs
    real(dp),                           intent(in), optional :: fwhm
    real(dp),                           intent(in), optional :: mb_eff
    character(len=*),                   intent(in), optional :: pixwin
    integer(i4b),                       intent(in), optional :: nside
    logical(lgt),                       intent(in), optional :: init_realspace
    class(comm_B_FIRAS),   pointer                              :: constructor

    integer(i4b)       :: i
    character(len=4)   :: nside_text
    character(len=512) :: dir
    real(dp)           :: phi_ecl, theta_ecl
    
    ! General parameters
    allocate(constructor)
    dir = trim(cpar%datadir) // '/'

    constructor%info        => info
    constructor%type        =  'FIRAS'    
    constructor%almFromConv = .false.
    constructor%nside_hires =  4*constructor%info%nside
    constructor%npix_hires  =  12*constructor%nside_hires**2

    ! Component specific parameters
    if (present(fwhm)) then
       allocate(constructor%b_l(0:constructor%info%lmax,constructor%info%nmaps))
       if (present(nside)) then
          call int2string(nside, nside_text)
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm, &
               & pixwin=trim(dir)//'/pixel_window_n'//nside_text//'.fits')
       else if (present(pixwin)) then
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm, &
               & pixwin=trim(dir)//'/'//trim(pixwin))
       else
          call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, fwhm=fwhm)
       end if
    else
       call read_beam(constructor%info%lmax, constructor%info%nmaps, constructor%b_l, &
            & beamfile=trim(cpar%ds_blfile(id_abs)), &
            & pixwin=trim(cpar%ds_pixwin(id_abs)))
    end if

    ! Multiply with main beam filling factor
    constructor%mb_eff = 1.d0; if (present(mb_eff)) constructor%mb_eff = mb_eff

    ! Initialize real-space profile
    call constructor%initBTheta(filename=cpar%ds_btheta_file(id_abs))
!    call compute_radial_beam(constructor%info%nmaps, constructor%b_l, constructor%b_theta)
    constructor%r_max = maxval(constructor%b_theta(1)%x)

    ! Precompute ecliptic rotation matrix
    phi_ecl   =        96.38d0  * pi/180.d0
    theta_ecl = (90.d0-29.81d0) * pi/180.d0
    call compute_euler_matrix_zyz(-phi_ecl, -theta_ecl, 0.d0, constructor%M_ecl)

    allocate(constructor%vecs(3,0:constructor%npix_hires-1))
    do i = 0, constructor%npix_hires-1
       call pix2vec_ring(constructor%nside_hires, i, constructor%vecs(:,i))
    end do
       
  end function constructor
  
  subroutine matmulB_firas(self, trans, map)
    implicit none
    class(comm_B_FIRAS), intent(in)    :: self
    logical(lgt),        intent(in)    :: trans
    class(comm_map),     intent(inout) :: map

    real(dp)           :: vec(3), theta, vec0(3), weight, w, vec_rot(3), vec1(3)
    real(dp)           :: M_rot(3,3), dalpha0, alpha0, dalpha1, alpha1
    real(dp)           :: phi0, phi1, phi_rot, theta_rot, M_sat(3,3), alpha, phi_ecl, theta_ecl
    integer(i4b)       :: i, j, k, l, n, nstep
    integer(i4b), allocatable, dimension(:) :: pixlist
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map),     pointer :: map_hires => null()
    real(dp), allocatable, dimension(:,:) :: map_full 

    nstep     = 11
    alpha0    =    3.19d0 * pi/180.d0
    dalpha0   =  -0.585d0 * pi/180.d0
    alpha1    =   -2.78d0 * pi/180.d0
    dalpha1   =   0.520d0 * pi/180.d0
    phi0      =    200.d0 * pi/180.d0
    phi1      =    350.d0 * pi/180.d0

    ! Create and distribute high-res map to all cores
    allocate(map_full(0:self%npix_hires-1,1), pixlist(0:self%npix_hires-1))
    info  => comm_mapinfo(self%info%comm, self%nside_hires, &
         & self%info%lmax, self%info%nmaps, self%info%pol)
    map_hires => comm_map(info)
    call map%alm_equal(map_hires)
    call map_hires%Y()
    call map_hires%bcast_fullsky_map(map_full)
    call map_hires%dealloc(); deallocate(map_hires)

    map%map = 0.d0
    do i = 1, self%info%np
       call pix2vec_ring(self%info%nside, self%info%pix(i), vec0)
       vec0 = matmul(transpose(self%M_ecl), vec0)
       call crossproduct(vec0, [0.d0,0.d0,1.d0], vec_rot)
       vec_rot   = vec_rot / sqrt(sum(vec_rot**2))
       call vec2ang(vec0, theta_ecl, phi_ecl)
       call vec2ang(vec_rot, theta_rot, phi_rot)
       call compute_euler_matrix_zyz(phi_rot, theta_rot, 0.d0, M_rot)
       weight = 0.d0
       do k = 1, nstep
          if (phi_ecl < phi0 .or. phi_ecl > phi1) then
             alpha = alpha0 + (k-1)*dalpha0
          else
             alpha = alpha1 + (k-1)*dalpha1
          end if
          call compute_euler_matrix_zyz(0.d0, 0.d0, alpha, M_sat)
          vec1 = matmul(self%M_ecl, matmul(M_rot, matmul(M_sat, &
               & matmul(transpose(M_rot), vec0))))
          call query_disc(self%nside_hires, vec1, self%r_max, pixlist, n)
          do j = 0, n-1
             vec = self%vecs(:,pixlist(j))
             call angdist(vec1, vec, theta)
             w = splint(self%b_theta(1), theta)
             map%map(i-1,1) = map%map(i-1,1) + map_full(pixlist(j),1) * w
             weight = weight + w
          end do
       end do
       map%map(i-1,1) = map%map(i-1,1) / weight
    end do

!!$    call map%writeFITS("chains_hke/test.fits")
!!$    call mpi_finalize(i)
!!$    stop

    deallocate(map_full, pixlist)

  end subroutine matmulB_firas

  subroutine matmulInvB_firas(self, trans, map)
    implicit none
    class(comm_B_FIRAS), intent(in)    :: self
    logical(lgt),     intent(in)    :: trans
    class(comm_map),  intent(inout) :: map

    integer(i4b) :: i, l, j

    write(*,*) 'matmulInvB not implemented in b_FIRAS'
    call mpi_finalize(i)
    stop

  end subroutine matmulInvB_firas

  subroutine updateBeam(self, b_l_norm, mb_eff) 
    implicit none
    class(comm_B_FIRAS),                intent(inout)           :: self
    real(dp),         dimension(0:,1:), intent(in),    optional :: b_l_norm
    real(dp),                           intent(in),    optional :: mb_eff

    if (present(mb_eff)) self%mb_eff = mb_eff
    
  end subroutine updateBeam

  subroutine crossproduct(vector1, vector2, crossvector)
    implicit none 
    real(dp), dimension(3), intent(in)  :: vector1, vector2
    real(dp), dimension(3), intent(out) :: crossvector
    crossvector=[vector1(2)*vector2(3)-vector1(3)*vector2(2), &
               & vector1(3)*vector2(1)-vector1(1)*vector2(3), &
               & vector1(1)*vector2(2)-vector1(2)*vector2(1) ]
  end subroutine crossproduct
  
end module comm_B_FIRAS_mod
