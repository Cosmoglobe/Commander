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
module comm_comp_mod
  use comm_cr_utils
  use comm_cr_precond_mod
  use comm_data_mod
  implicit none

!  private
!  public  :: comm_comp, ncomp, compList, update_mixing_matrices, comp_ptr!, dumpCompMaps
  
  !**************************************************
  !        Generic component class definition - top level class
  !**************************************************
  type, abstract :: comm_comp !commander components
     ! Linked list variables
     class(comm_comp), pointer :: nextLink => null()
     class(comm_comp), pointer :: prevLink => null()

     ! Data variables
     logical(lgt)       :: active
     integer(i4b)       :: npar, ncr, id, nmaps, myid, comm, numprocs, cg_unique_sampgroup
     character(len=512) :: label, class, type, unit, operation, init_from_HDF
     logical(lgt)       :: output
     real(dp)           :: nu_ref(3), RJ2unit_(3), nu_min, nu_max
     character(len=512), allocatable, dimension(:)   :: indlabel
     integer(i4b),       allocatable, dimension(:)   :: poltype
     real(dp),           allocatable, dimension(:)   :: theta_def
     real(dp),        allocatable, dimension(:,:) :: theta_steplen
     real(dp),        allocatable, dimension(:)   :: scale_sigma

     real(dp),           allocatable, dimension(:,:) :: p_gauss
     real(dp),           allocatable, dimension(:,:) :: p_uni
     integer(i4b),       allocatable, dimension(:)   :: smooth_scale
     real(dp),           allocatable, dimension(:)   :: nu_min_ind
     real(dp),           allocatable, dimension(:)   :: nu_max_ind
     logical(lgt),       allocatable, dimension(:)   :: active_samp_group

   contains
     ! Linked list procedures
     procedure :: nextComp    ! get the link after this link
     procedure :: prevComp    ! get the link before this link
     procedure :: setNextComp ! set the link after this link
     procedure :: addComp     ! add new link at the end

     ! Data procedures
     procedure                          :: initComp
     procedure(evalSED),       deferred :: S
     procedure(evalBand),      deferred :: getBand
     procedure(projectBand),   deferred :: projectBand
     procedure                          :: dumpSED
     procedure(dumpFITS),      deferred :: dumpFITS
     procedure(initHDFcomp),   deferred :: initHDFcomp
     procedure                          :: RJ2unit
     procedure(sampleSpecInd), deferred :: sampleSpecInd
     procedure                          :: CG_mask
     procedure(update_F_int),  deferred :: update_F_int
     procedure(updateMixmat),  deferred :: updateMixmat
  end type comm_comp

  abstract interface
     ! Evaluate SED in brightness temperature normalized to reference frequency
     function evalSED(self, nu, band, pol, theta)
       import comm_comp, dp, i4b
       class(comm_comp),        intent(in)           :: self
       real(dp),                intent(in), optional :: nu
       integer(i4b),            intent(in), optional :: band
       integer(i4b),            intent(in), optional :: pol
       real(dp), dimension(1:), intent(in), optional :: theta
       real(dp)                                      :: evalSED
     end function evalSED

     ! Return effective signal at given frequency band
     function evalBand(self, band, amp_in, pix, alm_out, det)
       import i4b, dp, comm_comp, lgt
       class(comm_comp),                             intent(in)            :: self
       integer(i4b),                                 intent(in)            :: band
       integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
       real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
       logical(lgt),                                 intent(in),  optional :: alm_out
       integer(i4b),                                 intent(in),  optional :: det
       real(dp),        dimension(:,:), allocatable                        :: evalBand
     end function evalBand

     ! Return component projected from map
     function projectBand(self, band, map, alm_in, det)
       import i4b, dp, comm_comp, comm_map, lgt
       class(comm_comp),                             intent(in)            :: self
       integer(i4b),                                 intent(in)            :: band
       class(comm_map),                              intent(in)            :: map
       logical(lgt),                                 intent(in), optional  :: alm_in
       integer(i4b),                                 intent(in), optional  :: det
       real(dp),        dimension(:,:), allocatable                        :: projectBand
     end function projectBand

!!$     ! Evaluate amplitude map in brightness temperature at reference frequency
!!$     function evalAmp(self, nside, nmaps, pix, x_1D, x_2D)
!!$       import comm_comp, dp, i4b
!!$       class(comm_comp),                          intent(in)           :: self
!!$       integer(i4b),                              intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:),   intent(in), optional :: pix
!!$       real(dp),                  dimension(:),   intent(in), optional :: x_1D
!!$       real(dp),                  dimension(:,:), intent(in), optional :: x_2D
!!$       real(dp),     allocatable, dimension(:,:)                       :: evalAmp
!!$     end function evalAmp
!!$
!!$     ! Evaluate mixing matrix in brightness temperature at reference frequency
!!$     function evalMixmat(self, nside, nmaps, pix)
!!$       import comm_comp, dp, i4b
!!$       class(comm_comp),                        intent(in)           :: self
!!$       integer(i4b),                            intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:), intent(in), optional :: pix
!!$       real(dp),     allocatable, dimension(:,:)                     :: evalMixmat
!!$     end function evalMixmat
!!$
!!$     ! Generate simulated component
!!$     function simComp(self, handle, nside, nmaps, pix)
!!$       import planck_rng, comm_comp, dp, i4b
!!$       class(comm_comp),                        intent(in)           :: self
!!$       type(planck_rng),                        intent(inout)        :: handle
!!$       integer(i4b),                            intent(in)           :: nside, nmaps
!!$       integer(i4b),              dimension(:), intent(in), optional :: pix
!!$       real(dp),     allocatable, dimension(:,:)                     :: simComp
!!$     end function simComp
!!$
!!$     ! Dump current sample to HDF chain files
!!$     subroutine dumpHDF(self, filename)
!!$       import comm_comp
!!$       class(comm_comp),                        intent(in)           :: self
!!$       character(len=*),                        intent(in)           :: filename
!!$     end subroutine dumpHDF
!!$
     ! Dump current sample to HEALPix FITS file
     subroutine dumpFITS(self, iter, chainfile, output_hdf, postfix, dir)
       import comm_comp, i4b, hdf_file, lgt
       class(comm_comp),                        intent(inout)        :: self
       integer(i4b),                            intent(in)           :: iter
       type(hdf_file),                          intent(in)           :: chainfile
       logical(lgt),                            intent(in)           :: output_hdf
       character(len=*),                        intent(in)           :: postfix
       character(len=*),                        intent(in)           :: dir
     end subroutine dumpFITS

     ! Initialize from HDF chain file
     subroutine initHDFcomp(self, cpar, hdffile, hdfpath)
       import comm_comp, comm_params, hdf_file
       class(comm_comp),                        intent(inout)        :: self
       type(comm_params),                       intent(in)           :: cpar
       type(hdf_file),                          intent(in)           :: hdffile
       character(len=*),                        intent(in)           :: hdfpath
     end subroutine initHDFcomp

     ! Sample spectral parameters
     subroutine sampleSpecInd(self, cpar, handle, id, iter)
       import comm_comp, comm_params, planck_rng, i4b
       class(comm_comp),                        intent(inout)        :: self
       type(comm_params),                       intent(in)           :: cpar
       type(planck_rng),                        intent(inout)        :: handle
       integer(i4b),                            intent(in)           :: id
       integer(i4b),                            intent(in)           :: iter
     end subroutine sampleSpecInd

     ! Update mixing matrices
     subroutine updateMixmat(self, theta, beta, band, df, par)
       import comm_comp, comm_map, dp, i4b, map_ptr
       class(comm_comp),                        intent(inout)           :: self
       class(comm_map), dimension(:),           intent(in),    optional :: theta
       real(dp),        dimension(:,:,:),       intent(in),    optional :: beta
       integer(i4b),                            intent(in),    optional :: band
       class(map_ptr), dimension(:),            intent(inout), optional :: df
       integer(i4b),                            intent(in),    optional :: par
     end subroutine updateMixmat

     ! Update band integration lookup tables
     subroutine update_F_int(self, band)
       import comm_comp, i4b
       class(comm_comp),                        intent(inout)        :: self
       integer(i4b),                            intent(in), optional :: band
     end subroutine update_F_int
       
  end interface


  type comp_ptr
     class(comm_comp), pointer :: p => null()
  end type comp_ptr

  !**************************************************
  !             Auxiliary variables
  !**************************************************

  integer(i4b)              :: n_dump = 10000
  real(dp),    dimension(2) :: nu_dump = [0.1d9, 300000.d9] 

  !**************************************************
  !             Internal module variables
  !**************************************************
  integer(i4b)              :: ncomp
  class(comm_comp), pointer :: compList => null(), currComp => null()
  
contains

  subroutine initComp(self, cpar, id, id_abs)
    implicit none
    class(comm_comp)               :: self
    type(comm_params),  intent(in) :: cpar
    integer(i4b),       intent(in) :: id, id_abs

    integer(i4b) :: i, j, n
    character(len=16), dimension(1000) :: comp_label


  end subroutine initComp

  subroutine dumpSED(self, unit)
    implicit none
    class(comm_comp), intent(in) :: self
    integer(i4b),     intent(in) :: unit

    integer(i4b) :: i
    real(dp)     :: nu, dnu, S

    
  end subroutine dumpSED


  function RJ2unit(self, pol)
    implicit none

    class(comm_comp), intent(in)           :: self
    integer(i4b),     intent(in)           :: pol
    real(dp)                               :: RJ2unit


    RJ2unit = 0d0

    
  end function RJ2unit

  subroutine CG_mask(self, samp_group, mask)
    implicit none

    class(comm_comp),               intent(in)    :: self
    integer(i4b),                   intent(in)    :: samp_group
    real(dp),         dimension(:), intent(inout) :: mask

    call cr_mask(self%id, self%active_samp_group(samp_group), mask)

  end subroutine CG_mask
  
  function nextComp(self)
    class(comm_comp) :: self
    class(comm_comp), pointer :: nextComp
    nextComp => self%nextLink
  end function nextComp

  function prevComp(self)
    class(comm_comp) :: self
    class(comm_comp), pointer :: prevComp
    prevComp => self%prevLink
  end function prevComp
  
  subroutine setNextComp(self,nextComp)
    class(comm_comp) :: self
    class(comm_comp), pointer :: nextComp
    self%nextLink => nextComp
  end subroutine setNextComp

  subroutine addComp(self,link)
    class(comm_comp), target  :: self
    class(comm_comp), pointer :: link

    class(comm_comp), pointer :: c => null()
    
    c => self
    do while (associated(c%nextLink))
       c => c%nextLink
    end do
    link%prevLink => c
    c%nextLink    => link
  end subroutine addComp

  subroutine update_mixing_matrices(band, update_F_int)
    implicit none
    integer(i4b), intent(in), optional :: band
    logical(lgt), intent(in), optional :: update_F_int
    integer(i4b) :: i, j, k
    logical(lgt) :: update_F

    class(comm_comp), pointer :: c => null()
    update_F =.false.; if (present(update_F_int)) update_F = update_F_int 

    c => compList
    do while (associated(c))
       if (present(band)) then
          if (update_F) call c%update_F_int(band)
          call c%updateMixmat(band=band)
       else
          if (update_F) call c%update_F_int
          call c%updateMixmat()
       end if
       c => c%nextLink
    end do
  end subroutine update_mixing_matrices

end module comm_comp_mod
