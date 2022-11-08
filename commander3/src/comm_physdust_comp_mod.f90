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
module comm_physdust_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_diffuse_comp_mod
  use comm_map_mod
  use comm_F_int_1D_mod
  use comm_data_mod
  use spline_2D_mod
  implicit none

  private
  public comm_physdust_comp

  !**************************************************
  !           Power-law component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_physdust_comp
     integer(i4b) :: num_nu, num_u, num_comp, num_u_int
     real(dp)     :: log_umax, gamma, alpha
     real(dp)     :: u_arr_min, u_arr_max, du_arr
     real(dp), allocatable, dimension(:)         :: log_dust_wav, dust_u, extcrv, extpol, scat, amps
     real(dp), allocatable, dimension(:,:,:,:,:) :: coeff_arr
   contains
     procedure :: S    => evalSED
  end type comm_physdust_comp

  interface comm_physdust_comp
     procedure constructor
  end interface comm_physdust_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor(cpar, id, id_abs)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_physdust_comp), pointer   :: constructor

    real(dp)     :: dummy
    real(dp), allocatable, dimension(:,:,:) :: comp_mat
    integer(i4b) ::  unit, num_nu
    integer(i4b) :: i, j, k, l, m, n, p, ierr
    type(comm_mapinfo), pointer :: info => null()
    real(dp)           :: par_dp
    integer(i4b), allocatable, dimension(:) :: sum_pix
    real(dp),    allocatable, dimension(:) :: sum_theta, sum_proplen, sum_nprop
    character(len=512) :: temptxt, partxt
    integer(i4b) :: smooth_scale, p_min, p_max
    class(comm_mapinfo), pointer :: info2 => null()
    class(comm_map),     pointer :: tp => null() 
    class(comm_map),     pointer :: tp_smooth => null() 


    ! General parameters
    allocate(constructor)


    constructor%npar         = 2
    allocate(constructor%poltype(constructor%npar))
    do i = 1, constructor%npar
       constructor%poltype(i)   = cpar%cs_poltype(i,id_abs)
    end do
    call constructor%initLmaxSpecind(cpar, id, id_abs)

    call constructor%initDiffuse(cpar, id, id_abs)

    ! Component specific parameters
    allocate(constructor%theta_def(1), constructor%p_gauss(2,1), constructor%p_uni(2,1))
    allocate(constructor%indlabel(1))
    constructor%theta_def(1) = cpar%cs_theta_def(1,id_abs)
    constructor%p_uni(:,1)   = cpar%cs_p_uni(id_abs,:,1)
    constructor%p_gauss(:,1) = cpar%cs_p_gauss(id_abs,:,1)
    constructor%indlabel(1)  = 'Umin'

    ! Initialize spectral index map
    info => comm_mapinfo(cpar%comm_chain, constructor%nside, constructor%lmax_ind, &
         & constructor%nmaps, constructor%pol)

    allocate(constructor%theta(constructor%npar))
    do i = 1, constructor%npar
       if (trim(cpar%cs_input_ind(i,id_abs)) == 'default' .or. trim(cpar%cs_input_ind(i,id_abs)) == 'none') then
          constructor%theta(i)%p => comm_map(info)
          constructor%theta(i)%p%map = constructor%theta_def(1)
       else
          ! Read map from FITS file, and convert to alms
          constructor%theta(i)%p => comm_map(info, trim(cpar%cs_input_ind(i,id_abs)))
       end if

       !convert spec. ind. pixel map to alms if lmax_ind >= 0
       if (constructor%lmax_ind >= 0) then
          ! if lmax >= 0 we can get alm values for the theta map
          call constructor%theta(i)%p%YtW_scalar
       end if
    end do

    ! Initialize SED templates
    unit                  = getlun()
    num_nu                = 1000    ! Number of frequencies in dust files
    constructor%num_nu    = num_nu
    constructor%num_u     = 11      ! Number of radiation field strengths
    constructor%u_arr_min = -0.5d0
    constructor%u_arr_max =  0.5d0
    constructor%du_arr    = (constructor%u_arr_max - constructor%u_arr_min)/(constructor%num_u - 1.d0)
    constructor%num_comp  = 4 ! Number of dust components

    allocate(constructor%coeff_arr(4,4,num_nu,constructor%num_u,constructor%num_comp))
    allocate(constructor%log_dust_wav(num_nu), constructor%dust_u(constructor%num_u))
    allocate(constructor%extcrv(num_nu), constructor%extpol(num_nu), constructor%scat(num_nu))
    allocate(constructor%amps(constructor%num_comp))

    constructor%log_umax  = cpar%cs_auxpar(1,id_abs)
    constructor%gamma     = cpar%cs_auxpar(2,id_abs)
    constructor%alpha     = cpar%cs_auxpar(3,id_abs)
    constructor%amps      = cpar%cs_auxpar(4:3+constructor%num_comp,id_abs)

    ! Make array of log U
    do i = 1, constructor%num_u
       constructor%dust_u(i) = -0.5d0 + constructor%du_arr*(i-1)
    enddo

    allocate(comp_mat(num_nu,constructor%num_u,constructor%num_comp))    
    do k = 1, constructor%num_comp
       open(unit,file=trim(cpar%cs_SED_template(k,id_abs)),recl=4096)
       do i = 1, num_nu
          read(unit,fmt='(1P26E11.3)') constructor%log_dust_wav(i), constructor%extcrv(i),&
               & constructor%extpol(i), constructor%scat(i), (comp_mat(i,j,k),j=1,constructor%num_u), &
                (dummy,j=1,constructor%num_u)
          constructor%log_dust_wav(i) = log(constructor%log_dust_wav(i))
       enddo
       call splie2_full_precomp(constructor%log_dust_wav,constructor%dust_u,&
            & log(comp_mat(:,:,k)), constructor%coeff_arr(:,:,:,:,k))
        close(59)
     enddo
     deallocate(comp_mat)


    ! Precompute mixmat integrator for each band
    allocate(constructor%F_int(3,numband,0:constructor%ndet))
    do k = 1, 3
       do i = 1, numband
          do j = 0, data(i)%ndet
             if (k > 1) then
                if (constructor%nu_ref(k) == constructor%nu_ref(k-1)) then
                   constructor%F_int(k,i,j)%p => constructor%F_int(k-1,i,j)%p
                   cycle
                end if
             end if
             constructor%F_int(k,i,j)%p => comm_F_int_1D(constructor, data(i)%bp(j)%p, k)
          end do
       end do
    end do

    call constructor%initPixregSampling(cpar, id, id_abs)
    ! Init alm 
    if (constructor%lmax_ind >= 0) call constructor%initSpecindProp(cpar, id, id_abs)

    ! Initialize mixing matrix
    call constructor%updateMixmat

  end function constructor

  ! Definition:
  !    SED  = (nu/nu_ref)**beta
  ! where 
  !    beta = theta(1)
  function evalSED(self, nu, band, pol, theta)
    implicit none
    class(comm_physdust_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED

    integer(i4b) :: i, j, num_u_int
    real(dp) :: SED_physdust, SED_norm, wav_in, wav_ref, c
    real(dp) :: log_umin, umin, log_umax, umax, uval, du, fdu

    if (nu < 2d9) then
       evalSED = 0.d0
       return
    end if
     
    num_u_int = 100 ! How many dU steps to perform
    c         = 2.99792458d8
     
    log_umin     = theta(1)
    log_umax     = self%log_umax
    SED_physdust = 0.d0
    SED_norm     = 0.d0
    wav_in       = c/nu               * 1d6 ! In um
    wav_ref      = c/self%nu_ref(pol) * 1d6 ! In um
    umin         = 10**log_umin
    umax         = 10**log_umax
     do i = 1, self%num_comp
        ! Delta function component
        SED_physdust = SED_physdust + (1.d0-self%gamma)* &
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,&
             & self%coeff_arr(:,:,:,:,i), log(wav_in), log_umin)))
        SED_norm = SED_norm + (1.d0-self%gamma)*&
             & (self%amps(i)*exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
             log(wav_ref),log_umin)))

        ! Distribution of U values
        ! See, e.g., Aniano et al 2012, Eqs. 9 & 10
        if (self%gamma /= 0.d0) then
           du = umin*((umax/umin)**(1.d0/(num_u_int-1.d0)) - 1d0) ! Formally dU/U
           do j = 1, num_u_int ! integrate over U distribution
              uval = umin*(umax/umin)**((j-1.d0)/(num_u_int-1.d0))
              if (self%alpha /= 1.d0) then
                 fdu = uval**(1.d0-self%alpha)*du*self%gamma*&
                      & (self%alpha-1.d0)/(umin**(1.d0-self%alpha) - umax**(1.d0-self%alpha))
              else
                 ! Special case of alpha = 1 to ensure integrability
                 fdu = du*self%gamma/log(umax/umin)
              endif
              SED_physdust = SED_physdust + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   & log(wav_in),log10(uval)))*fdu
              SED_norm = SED_norm + self%amps(i)* &
                   & exp(splin2_full_precomp_irreg(self%log_dust_wav,self%dust_u,self%coeff_arr(:,:,:,:,i), &
                   log(wav_ref),log10(uval)))*fdu
           enddo
        endif
     enddo
     evalSED = (SED_physdust / SED_norm) * (self%nu_ref(pol)/nu)**3 ! Normalize to reference in T units

  end function evalSED
  
end module comm_physdust_comp_mod
