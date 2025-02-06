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
module comm_chisq_mod
  use comm_comp_interface_mod
  implicit none

contains

  subroutine compute_chisq(comm, chisq_map, chisq_fullsky, mask, maskpath, lowres_eval, band_list, evalpol)
    implicit none
    integer(i4b),                   intent(in)              :: comm
    logical(lgt),                   intent(in),    optional :: lowres_eval
    logical(lgt),                   intent(in),    optional :: evalpol
    !character(len=512),             intent(in),    optional :: evalsig
    !logical(lgt),                   intent(in),    optional :: udgrade_chisq
    class(comm_map),                intent(inout), optional :: chisq_map
    real(dp),                       intent(out),   optional :: chisq_fullsky
    type(map_ptr),   dimension(1:), intent(in),    optional :: mask
    character(len=512),             intent(in),    optional :: maskpath
    integer(i4b), dimension(:),     intent(in),    optional :: band_list

    integer(i4b) :: i, j, k, p, ierr, nmaps, nbands
    integer(i4b), dimension(:), allocatable :: bandlist
    real(dp)     :: t1, t2
    logical(lgt) :: lowres
    class(comm_map), pointer :: res, res_lowres => null(), res_lowres_temp, chisq_sub, mask_tmp
    class(comm_mapinfo), pointer :: info, info_lowres


  end subroutine compute_chisq


  subroutine compute_jeffreys_prior(c, df, pol, par, chisq_jeffreys)
    implicit none
    class(comm_diffuse_comp),       intent(in)              :: c
    type(map_ptr),   dimension(1:), intent(in)              :: df
    integer(i4b),                   intent(in)              :: pol, par
    real(dp),                       intent(out)             :: chisq_jeffreys

    integer(i4b) :: i, j, k, p, ierr, nmaps
    class(comm_map),     pointer :: map
    class(comm_mapinfo), pointer :: info 

      chisq_jeffreys = 0d0


  end subroutine compute_jeffreys_prior



  function compute_residual(band, exclude_comps, cg_samp_group) result (res)
    implicit none

    integer(i4b),                     intent(in)           :: band
    character(len=512), dimension(:), intent(in), optional :: exclude_comps
    integer(i4b),                     intent(in), optional :: cg_samp_group
    class(comm_map),    pointer                            :: res

    integer(i4b) :: i
    logical(lgt) :: skip
    real(dp)     :: t1, t2, t3, t4
    class(comm_comp),    pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm
    integer(i4b), allocatable, dimension(:)   :: pix
    integer(i4b) :: ierr
    logical(lgt) :: nonzero
    class(comm_map), pointer :: ptsrc
    
    ! Initialize to full data set
    res   => comm_map(data(band)%info)  ! Diffuse

  end function compute_residual

  subroutine subtract_fiducial_CMB_dipole(band, map)
    implicit none
    integer(i4b),    intent(in)    :: band
    class(comm_map), intent(inout) :: map

    integer(i4b)        :: i, j, l, m
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer     :: dipole
    real(dp),            allocatable, dimension(:,:) :: alm
    class(comm_comp),    pointer :: c


  end subroutine subtract_fiducial_CMB_dipole

  subroutine add_fiducial_CMB_dipole(info, RJ2unit, alm)
    implicit none
    class(comm_mapinfo),                   intent(in)    :: info
    real(dp),                              intent(in)    :: RJ2unit
    real(dp),            dimension(0:,1:), intent(inout) :: alm

    integer(i4b)        :: i, j, l, m


  end subroutine add_fiducial_CMB_dipole

  subroutine output_signals_per_band(outdir, postfix)
    implicit none
    character(len=*), intent(in) :: outdir, postfix
    
    integer(i4b) :: i
    logical(lgt) :: skip
    character(len=1024) :: filename
    class(comm_comp), pointer :: c
    class(comm_map),  pointer :: out
    real(dp),     allocatable, dimension(:,:) :: map, alm
    integer(i4b), allocatable, dimension(:)   :: pix
    

  end subroutine output_signals_per_band

  subroutine get_sky_signal(band, det, map_out, mono, cmbmap, abscal_comps, gainmap)
    implicit none
    integer(i4b),    intent(in)     :: band, det
    class(comm_map), pointer        :: map_out
    logical(lgt),    optional       :: mono 
    class(comm_map), optional       :: cmbmap
    character(len=512), intent(in), optional :: abscal_comps
    class(comm_map), pointer, intent(inout), optional    :: gainmap

    integer(i4b) :: i, j, k, n
    logical(lgt) :: skip, mono_, calmap
    real(dp)     :: rms_EE2_prior
    class(comm_map),  pointer :: map_diff, cmbmap_band, gaindiff
    class(comm_comp), pointer :: c
    real(dp),     allocatable, dimension(:,:) :: map, alm
    real(dp),                  dimension(5)   :: P_quad
    character(len=16),         dimension(100) :: abscal_labels
    

  end subroutine get_sky_signal


  subroutine compute_marginal(mixing, data, invN, marg_map, marg_fullsky)
    implicit none
    
    real(c_double),  intent(in),    dimension(:,:,0:):: mixing   !(nbands,ncomp,npix) mixing matrix
    real(c_double),  intent(in),    dimension(:,0:)  :: invN     !(nbands,npix) inverse noise matrix
    real(c_double),  intent(in),    dimension(:,:)   :: data     !(nbands,npix) data matrix
    class(comm_map), intent(inout), optional         :: marg_map
    real(dp),        intent(out),   optional         :: marg_fullsky
    integer :: i, j, k, l, p, ierr, nb, npix, nc
    logical :: temp_bool
    double precision     :: temp_marg
    double precision, allocatable, dimension(:)    :: MNd    ! (M.T*invN*M)
    double precision, allocatable, dimension(:)    :: M_d    ! (M.T*invN*M)^-1 * (M.T*invN*d)
    double precision, allocatable, dimension(:,:)  :: MN     ! M.T*ivnN
    double precision, allocatable, dimension(:,:)  :: MNM    ! M.T*ivnN*M (and its inverse)
    double precision, allocatable, dimension(:,:)  :: invmat ! matrix to invert MNM
    double precision, allocatable, dimension(:)    :: temp_arr ! array to flip rows/columns in matrices

  end subroutine compute_marginal

end module comm_chisq_mod
