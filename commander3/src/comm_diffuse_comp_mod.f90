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
module comm_diffuse_comp_mod ! only interfaces in this file, accompanying smod.f90 file
  use comm_comp_mod
  use comm_Cl_mod
  use comm_F_mod
  implicit none

  private
  public comm_diffuse_comp, add_to_npre, updateDiffPrecond, initDiffPrecond, applyDiffPrecond, &
       & res_smooth, res_lowres, dust_lowres, hotpah_lowres, rms_smooth, print_precond_mat, nullify_monopole_amp, &
       & get_monopole_amp, set_monopole_amp, recompute_diffuse_precond, precond_type, diff_ptr
  
  !**************************************************
  !            Diffuse component class - subclass under component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nx, x0, ndet
     logical(lgt)       :: pol, output_mixmat, output_EB, apply_jeffreys, almsamp_pixreg, priorsamp_local
     logical(lgt)       :: output_localsamp_maps
     integer(i4b)       :: lmin_amp, lmax_amp, lmax_ind, lmax_prior, lpiv, l_apod, lmax_pre_lowl
     integer(i4b)       :: lmax_def, nside_def, ndef, nalm_tot, sample_first_niter

     real(dp),     allocatable, dimension(:,:)   :: sigma_priors, steplen, chisq_min
     real(dp),     allocatable, dimension(:,:,:,:)   :: L
     integer(i4b), allocatable, dimension(:,:)   :: corrlen     
     logical(lgt),    dimension(:), allocatable :: L_read, L_calculated

     integer(i4b)       :: cg_samp_group_md         ! in the case we prior sample the monopole of the component amplitude
                                                    ! we will need to re-estimate the monopoles afterwards. This parameter
                                                    ! stores the CG sampling group number for the CG group with "md" and 
                                                    ! nothing more. (= -1 if no pure "md" CG group)
     real(dp)           :: cg_scale(1:3)                                  ! make the cg_scale vary between T and P
     real(dp)           :: latmask, fwhm_def, test
     real(dp),           allocatable, dimension(:,:)   :: cls
     real(dp),           allocatable, dimension(:,:,:) :: F_mean
     character(len=128), allocatable, dimension(:,:)   :: pol_lnLtype     ! {'chisq', 'ridge', 'marginal', 'prior'}
     integer(i4b),       allocatable, dimension(:,:)   :: lmax_ind_pol    ! lmax per poltype per spec. ind  
     integer(i4b),       allocatable, dimension(:,:)   :: lmax_ind_mix    ! equal to lmax_ind_pol, but 0 where lmax=-1 and fullsky pixreg
     integer(i4b),       allocatable, dimension(:,:)   :: pol_pixreg_type ! {1=fullsky, 2=single_pix, 3=pixel_regions}
     integer(i4b),       allocatable, dimension(:,:)   :: nprop_uni       ! limits on nprop
     integer(i4b),       allocatable, dimension(:,:)   :: npixreg          ! number of pixel regions
     integer(i4b),       allocatable, dimension(:,:,:) :: ind_pixreg_arr  ! number of pixel regions
     real(dp),           allocatable, dimension(:,:,:) :: theta_pixreg    ! thetas for pixregs, per poltype, per ind.
     real(dp),           allocatable, dimension(:,:,:) :: theta_pixreg_buff    ! thetas for pixregs, per poltype, per ind.
     real(dp),           allocatable, dimension(:,:,:) :: prior_pixreg    ! thetas for pixregs, per poltype, per ind.
     real(dp),           allocatable, dimension(:,:,:) :: proplen_pixreg  ! proposal length for pixregs
     real(dp),           allocatable, dimension(:,:,:) :: pixreg_priors   ! individual priors for pixel regions
     integer(i4b),       allocatable, dimension(:,:,:) :: nprop_pixreg    ! number of proposals for pixregs
     integer(i4b),       allocatable, dimension(:,:,:) :: npix_pixreg     ! number of pixels per pixel region
     logical(lgt),       allocatable, dimension(:,:)   :: pol_sample_nprop   ! sample the corr. length in first iteration
     logical(lgt),       allocatable, dimension(:,:)   :: pol_sample_proplen ! sample the prop. length in first iteration
     logical(lgt),       allocatable, dimension(:)     :: spec_corr_convergence ! local sampling pixel regions until a minimum running correlation is achieved
     real(dp),           allocatable, dimension(:)     :: spec_corr_limit ! local sampling running correlation limit

     logical(lgt),       allocatable, dimension(:,:)   :: first_ind_sample
     logical(lgt),       allocatable, dimension(:,:,:) :: fix_pixreg
     real(dp),           allocatable, dimension(:,:,:) :: theta_prior ! prior if pixel region (local sampler) == prior
     class(map_ptr),     allocatable, dimension(:)     :: pol_ind_mask
     class(map_ptr),     allocatable, dimension(:)     :: pol_proplen
     class(map_ptr),     allocatable, dimension(:)     :: ind_pixreg_map   !map with defined pixelregions
     class(map_ptr),     allocatable, dimension(:)     :: pol_nprop
     class(map_ptr),     allocatable, dimension(:)     :: spec_mono_mask
     logical(lgt),       allocatable, dimension(:)     :: spec_mono_combined
     character(len=512), allocatable, dimension(:)     :: spec_mono_type
     character(len=512), allocatable, dimension(:)     :: spec_mono_freeze
     class(comm_B_bl_ptr), allocatable, dimension(:)   :: B_pp_fr
     class(comm_B_bl_ptr), allocatable, dimension(:)   :: B_smooth_amp, B_smooth_specpar

     character(len=512) :: mono_prior_type, mono_prior_band
     real(dp)           :: mono_prior_gaussian_mean, mono_prior_gaussian_rms, mono_prior_fwhm
     integer(i4b)       :: mono_prior_nside, mono_prior_Nthresh
     real(dp), allocatable, dimension(:) :: mono_prior_threshold
     class(comm_map),               pointer     :: mono_prior_map => null()
     class(comm_map),               pointer     :: mask => null()
     class(comm_map),               pointer     :: procmask => null()
     class(map_ptr),  dimension(:), allocatable :: indmask
     class(comm_map),               pointer     :: defmask => null()
     class(comm_map),               pointer     :: priormask => null()
     class(comm_map),               pointer     :: x => null()           ! Spatial parameters
     real(dp)                                   :: x_scale !overall scaling parameter for component
     class(comm_map),               pointer     :: x_smooth => null()    ! Spatial parameters
     class(comm_map),               pointer     :: mu => null()          ! Spatial prior mean
     class(comm_B),                 pointer     :: B_out => null()       ! Output beam
     class(comm_B),                 pointer     :: B_mono_prior => null() ! monopole prior beam
     class(comm_Cl),                pointer     :: Cl => null()          ! Power spectrum
     class(map_ptr),  dimension(:), allocatable :: theta        ! Spectral parameters
     class(map_ptr),  dimension(:), allocatable :: theta_smooth ! Spectral parameters
     type(map_ptr),   dimension(:,:), allocatable :: F          ! Mixing matrix
     type(map_ptr),   dimension(:,:), allocatable :: dF         ! Derivative of mixing matrix
     real(dp),        dimension(:,:), allocatable :: invM_lowl  ! (0:nalm-1,0:nalm-1)
     real(dp),        dimension(:,:), allocatable :: Z_def      ! (0:nalm-1,ndef)
     real(dp),        dimension(:,:), allocatable :: invM_def   ! (0:nalm-1,0:nalm-1)
     logical(lgt),    dimension(:,:), allocatable :: F_null     ! Don't allocate space for null mixmat's
     type(F_int_ptr), dimension(:,:,:), allocatable :: F_int        ! SED integrator
     integer(i4b) :: ntab
     real(dp), allocatable, dimension(:,:) :: SEDtab        ! (2+npar_tab, nbin)
     real(dp), allocatable, dimension(:,:) :: SEDtab_buff   ! 
     real(dp)                              :: SEDtab_prior  ! (npar_tab), Single value for MH proposals, per comp
   contains
     procedure :: initDiffuse
     procedure :: initPixregSampling
     procedure :: initSpecindProp
     procedure :: initLmaxSpecind
     procedure :: updateMixmat  => updateDiffuseMixmat
     procedure :: update_F_int  => updateDiffuseFInt
!!$     procedure :: dumpHDF    => dumpDiffuseToHDF
     procedure :: getBand       => evalDiffuseBand
     procedure :: projectBand   => projectDiffuseBand
     procedure :: dumpFITS      => dumpDiffuseToFITS
     procedure :: initHDFComp   => initDiffuseHDF
     procedure :: sampleSpecInd => sampleDiffuseSpecInd
     procedure :: updateLowlPrecond
     procedure :: applyLowlPrecond
     procedure :: updateDeflatePrecond
     procedure :: applyDeflatePrecond
     procedure :: applyMonoDipolePrior
  end type comm_diffuse_comp

  type diff_ptr
     class(comm_diffuse_comp), pointer :: p => null()
  end type diff_ptr
  
  integer(i4b) :: npre      =  0
  integer(i4b) :: lmax_pre  = -1
  integer(i4b) :: nside_pre = 1000000
  integer(i4b) :: nmaps_pre = -1
  logical(lgt) :: recompute_diffuse_precond = .true.
  logical(lgt) :: output_cg_eigenvals
  logical(lgt), private :: only_pol, only_I
  character(len=512) :: outdir, precond_type
  integer(i4b),        allocatable, dimension(:) :: ind_pre
  class(comm_mapinfo), pointer                   :: info_pre => null()
  class(diff_ptr),     allocatable, dimension(:) :: diffComps

  character(len=24), private :: operation

  ! Variables for non-linear search
  class(comm_diffuse_comp), pointer,       private :: c_lnL => null()
  integer(i4b),                            private :: k_lnL, p_lnL, id_lnL
  real(dp), allocatable, dimension(:),     private :: a_lnL
  real(dp), allocatable, dimension(:),     private :: theta_lnL        
  real(dp), allocatable, dimension(:,:),   private :: buffer_lnL        
  logical(lgt),                            private :: apply_mixmat = .true.
  type(map_ptr),               allocatable, dimension(:) :: res_smooth
  type(map_ptr),               allocatable, dimension(:) :: res_lowres
  type(map_ptr),               allocatable, dimension(:) :: dust_lowres
  type(map_ptr),               allocatable, dimension(:) :: hotpah_lowres
  class(comm_N_ptr),           allocatable, dimension(:) :: rms_smooth
  
interface

  module subroutine initDiffuse(self, cpar, id, id_abs)
    !
    ! Routine that initializes a diffuse type component. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! cpar: Commander parameter type
    !       Incudes all information from the parameter file
    !
    ! id: integer
    !       Integer ID of the diffuse component with respect to the active components
    !
    ! id_abs: integer
    !       Integer ID of the diffuse component with respect to all components defined in the parameter file
    !       (and also the id in the 'cpar' parameter)
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    
  end subroutine initDiffuse

  module subroutine initLmaxSpecind(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs


  end subroutine initLmaxSpecind

  module subroutine initPixregSampling(self, cpar, id, id_abs)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs


  end subroutine initPixregSampling

  module subroutine initSpecindProp(self,cpar, id, id_abs)
    !
    !  Subroutine that initializes the spectral index proposal matrix
    !
    !  Arguments:
    !  ----------
    !  self: comm_diffuse_component
    !       Object that contains all of the diffuse component properties
    !  cpar: comm_params
    !       Object that contains parameters input to Commander from parameter
    !       file.
    !  id: int
    !       index of the components, out of only the ones being used
    !  id_abs: int
    !       index of the component, out of all, whether they are used or not. 
    !  Returns:
    !  --------
    !  None, but modifies self by changing the proposal matrices.
    !
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id, id_abs

    
  end subroutine initSpecindProp


  module subroutine initDiffPrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm


  end subroutine initDiffPrecond

  module subroutine initDiffPrecond_diagonal(comm)
    implicit none
    integer(i4b),                intent(in) :: comm


  end subroutine initDiffPrecond_diagonal


  module subroutine initDiffPrecond_pseudoinv(comm)
    implicit none
    integer(i4b),                intent(in) :: comm


  end subroutine initDiffPrecond_pseudoinv


  module subroutine updateDiffPrecond(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update


  end subroutine updateDiffPrecond

  module subroutine updateDiffPrecond_diagonal(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update

       
  end subroutine updateDiffPrecond_diagonal


  module subroutine updateDiffPrecond_pseudoinv(samp_group, force_update)
    implicit none
    integer(i4b), intent(in) :: samp_group
    logical(lgt), intent(in) :: force_update


  end subroutine updateDiffPrecond_pseudoinv
    
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  module subroutine updateDiffuseMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta
    real(dp),  dimension(:,:,:),               intent(in),    optional :: beta  ! Not used here
    integer(i4b),                              intent(in),    optional :: band
    class(map_ptr), dimension(:),              intent(inout), optional :: df    ! Derivative of mixmat with respect to parameter par; for Jeffreys prior
    integer(i4b),                              intent(in),    optional :: par   ! Parameter ID for derivative


  end subroutine updateDiffuseMixmat



  module function evalDiffuseBand(self, band, amp_in, pix, alm_out, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: res


  end function evalDiffuseBand

  ! Return component projected from map
  module function projectDiffuseBand(self, band, map, alm_in, det) result(res)
    implicit none
    class(comm_diffuse_comp),                     intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: res


  end function projectDiffuseBand


  module subroutine applyDiffPrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x


  end subroutine applyDiffPrecond


  module subroutine applyDiffPrecond_diagonal(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x


  end subroutine applyDiffPrecond_diagonal


  module subroutine applyDiffPrecond_pseudoinv(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x


  end subroutine applyDiffPrecond_pseudoinv


  
  ! Dump current sample to HEALPix FITS file
  module subroutine dumpDiffuseToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    !
    ! Routine that writes a diffuce component to FITS (and HDF) files. 
    !
    ! Arguments:
    ! self: comm_diffuse_comp 
    !       Diffuse type component
    !
    ! iter: integer
    !       Sample number in the Gibb's chain.
    !
    ! chainfile: hdf_file
    !       HDF file to write the component to
    !
    ! output_hdf: logical
    !       Logical parameter to tell whether or not to write the component to the specified HDF file
    !
    ! postfix: string
    !       A string label to be added to the end of FITS-files.
    !       (default format: cXXXX_kYYYYYY; XXXX = chain number, YYYYYY = sample number)
    !
    ! dir: string
    !       Output directory to which output is written
    !
    ! Returns:
    !       The diffuse component parameter is returned (self).
    !       Any other changes are done internally
    !
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir


  end subroutine dumpDiffuseToFITS

  ! Dump current sample to HEALPix FITS file
  module subroutine initDiffuseHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_diffuse_comp),  intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath
  end subroutine initDiffuseHDF
  
  module subroutine add_to_npre(n, nside, lmax, nmaps)
    implicit none
    integer(i4b), intent(in) :: n, nside, lmax, nmaps
  end subroutine add_to_npre

  ! Sample spectral parameters. This subroutine is obsolete, it is defined and used in comm_nonlin_mod instead
  module subroutine sampleDiffuseSpecInd(self, cpar, handle, id, iter)
    implicit none
    class(comm_diffuse_comp),                intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter
  end subroutine sampleDiffuseSpecInd

  module function lnL_diffuse_multi(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: lnL_diffuse_multi
  end function lnL_diffuse_multi


  
  module subroutine print_precond_mat
    implicit none
  end subroutine print_precond_mat

  module subroutine updateDiffuseFInt(self, band)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: band
  end subroutine updateDiffuseFInt

  module subroutine updateLowlPrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
  end subroutine updateLowlPrecond

  module subroutine applyLowlPrecond(self, alm, alm0)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(inout)  :: alm
    real(dp),                 dimension(0:,1:), intent(in)     :: alm0
  end subroutine applyLowlPrecond


  module subroutine updateDeflatePrecond(self)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
  end subroutine updateDeflatePrecond



  module subroutine applyDeflatePrecond(self, alm, Qalm)
    implicit none
    class(comm_diffuse_comp),                   intent(in)     :: self
    real(dp),                 dimension(0:,1:), intent(in)     :: alm
    real(dp),                 dimension(0:,1:), intent(out)    :: Qalm
  end subroutine applyDeflatePrecond

  module subroutine setup_needlets(info, nside_def, defmask, Z, ndef)
    implicit none
    class(comm_mapinfo), pointer,          intent(in)  :: info
    integer(i4b),                          intent(in)  :: nside_def
    class(comm_map),                       intent(in)  :: defmask
    real(dp),            dimension(0:,1:), intent(out) :: Z
    integer(i4b),                          intent(out) :: ndef
  end subroutine setup_needlets

  module subroutine applyMonoDipolePrior(self, handle)
    implicit none
    class(comm_diffuse_comp), intent(inout)          :: self
    type(planck_rng),         intent(inout)          :: handle
  end subroutine applyMonoDipolePrior

  module subroutine nullify_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band
  end subroutine nullify_monopole_amp

  module function get_monopole_amp(band)
    implicit none
    character(len=*), intent(in) :: band
    real(dp)                     :: get_monopole_amp
  end function get_monopole_amp

  module subroutine set_monopole_amp(band, mono)
    implicit none
    character(len=*), intent(in) :: band
    real(dp),         intent(in) :: mono
  end subroutine set_monopole_amp

end interface

end module comm_diffuse_comp_mod
