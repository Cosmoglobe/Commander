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
module comm_line_comp_mod
  use comm_comp_interface_mod
  implicit none

  private
  public comm_line_comp

  !**************************************************
  !           Line emission component
  !**************************************************
  type, extends (comm_diffuse_comp) :: comm_line_comp
     integer(i4b)                            :: ref_band
     real(dp)                                :: line2RJ_ref
     integer(i4b), allocatable, dimension(:) :: ind2band
     real(dp),     allocatable, dimension(:) :: line2RJ
   contains
     procedure :: S    => evalSED_line
     procedure :: sampleSpecInd => sampleLineRatios
  end type comm_line_comp

  interface comm_line_comp
     procedure constructor_line
  end interface comm_line_comp

contains

  !**************************************************
  !             Routine definitions
  !**************************************************
  function constructor_line(cpar, id, id_abs) result(c)
    implicit none
    type(comm_params),   intent(in) :: cpar
    integer(i4b),        intent(in) :: id, id_abs
    class(comm_line_comp), pointer   :: c

    integer(i4b) :: i, j, k, l, m, nline, b, n, ierr
    real(dp)     :: f
    logical(lgt) :: ref_exist
    character(len=512), allocatable, dimension(:) :: label
    real(dp),           allocatable, dimension(:) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:) :: poltype
    type(comm_mapinfo), pointer :: info
    
    ! General parameters
    allocate(c)
    c%npar = 0 !temporary value so that lmax_ind is correcty set (to 0) in initDiffuse
    call c%initDiffuse(cpar, id, id_abs)

    ! Read line template file
    call read_line_template(trim(cpar%cs_SED_template(1,id_abs)), &
         & nline, label, mu, sigma, line2RJ, poltype)

    ! Check how many lines are included in current run
    n         = 0
    ref_exist = .false.
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) n = n+1
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) ref_exist = .true.
    end do
    if (.not. ref_exist) call report_error("Line component reference band does not exist, need "//trim(cpar%cs_band_ref(id_abs)))

    allocate(c%ind2band(n))
    c%npar = n

    allocate(c%lmax_ind_pol(3,c%npar))
    c%lmax_ind_pol = 0 !always fullsky (lmax=0) for line component
    if (allocated(c%lmax_ind_mix)) deallocate(c%lmax_ind_mix)
    allocate(c%lmax_ind_mix(3,c%npar))
    c%lmax_ind_mix = 0 !always fullsky (lmax=0) for line component
    allocate(c%pol_pixreg_type(3,c%npar))
    c%pol_pixreg_type = 0

    allocate(c%theta_def(n), c%p_gauss(2,n), c%p_uni(2,n))
    allocate(c%poltype(n), c%indlabel(n), c%line2RJ(n))
    n         = 0
    do i = 1, numband
       do j = 1, nline
          if (trim(label(j)) == trim(data(i)%label)) then
             n = n+1
             c%ind2band(n)  = i
             c%theta_def(n) = mu(j)
             c%p_gauss(1,n) = mu(j)
             c%p_gauss(2,n) = sigma(j)
             c%p_uni(1,n)   = -100.d30   !mu(j)-5*sigma(j)
             c%p_uni(2,n)   =  100.d30   !mu(j)+5*sigma(j)
             c%poltype(n)   = poltype(j)
             c%indlabel(n)  = label(j)
             c%line2RJ(n)   = line2RJ(j)
             exit
          end if
       end do
       if (trim(data(i)%label) == trim(cpar%cs_band_ref(id_abs))) then
          c%ref_band    = i
          c%line2RJ_ref = c%line2RJ(n)
       end if
    end do

    ! Update reference band unit conversion
    do i = 1, c%x%info%nmaps
       c%RJ2unit_(i) = 1.d0 / c%line2RJ_ref
       c%x%map(:,i)  = c%x%map(:,i) / c%RJ2unit_(i)
       c%x%alm(:,i)  = c%x%alm(:,i) / c%RJ2unit_(i)
    end do

    ! Initialize spectral index maps
    info => comm_mapinfo(cpar%comm_chain, c%nside, c%lmax_ind, &
         & c%nmaps, c%pol)

    allocate(c%theta(n))
    do i = 1, n
       c%theta(i)%p     => comm_map(info)
       c%theta(i)%p%map = c%theta_def(i)
       if (c%lmax_ind >= 0) call c%theta(i)%p%YtW_scalar
    end do


    ! Precompute mixmat integrator for each band
    allocate(c%F_int(3,numband,0:c%ndet))
    j = 1
    do l = 1, 3
       do i = 1, numband
          if (l > 1) then
             do m = 0,c%ndet
                c%F_int(l,i,m)%p => c%F_int(l-1,i,m)%p
             end do
             cycle
          end if
          if (any(c%ind2band == i)) then
             do k = 0, data(i)%ndet
!                write(*,*) 'line disabled'
                c%F_int(l,i,k)%p => comm_F_line(c, data(i)%bp(k)%p, .true., &
                     & c%line2RJ(j) / c%line2RJ_ref * data(i)%RJ2data(k), j)
             end do
             j = j+1
          else
             do k = 0, data(i)%ndet
!                write(*,*) 'line disabled'
                c%F_int(l,i,k)%p => comm_F_line(c, data(i)%bp(k)%p, .false., 0.d0, j)
             end do
          end if
       end do
    end do
    
    ! Initialize mixing matrix
    if (trim(cpar%init_chain_prefix) == 'none' &
         & .or. trim(c%init_from_HDF) == 'none') &
         & call c%updateMixmat

    deallocate(label, mu, sigma, line2RJ, poltype)

  end function constructor_line

  ! Definition:
  !    SED  = delta_{band,
  function evalSED_line(self, nu, band, pol, theta)
    class(comm_line_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSED_line

    integer(i4b) :: i, ind

    if (band == self%ref_band) then
       evalSED_line = 1.d0
    else 
       do i = 1, self%npar
          if (band == self%ind2band(i)) exit
       end do
       if (i > self%npar) then
          evalSED_line = 0.d0
       else
          evalSED_line = theta(i) * self%line2RJ(i) / self%line2RJ_ref
       end if
    end if

  end function evalSED_line

  subroutine read_line_template(filename, nline, label, mu, sigma, line2RJ, poltype)
    implicit none
    character(len=*),                              intent(in)  :: filename
    integer(i4b),                                  intent(out) :: nline
    character(len=512), allocatable, dimension(:), intent(out) :: label
    real(dp),           allocatable, dimension(:), intent(out) :: mu, sigma, line2RJ
    integer(i4b),       allocatable, dimension(:), intent(out) :: poltype

    integer(i4b)        :: i, unit
    character(len=1024) :: line

    unit  = getlun()

    ! Find number of lines
    nline = 0
    open(unit, file=trim(filename), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
    end do
1   close(unit)
    
    allocate(label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline))
    open(unit, file=trim(filename), recl=1024)
    nline = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       nline = nline+1
       read(line,*) label(nline), mu(nline), sigma(nline), line2RJ(nline), poltype(nline)
    end do
2   close(unit)

  end subroutine read_line_template

  ! Sample line ratios
  subroutine sampleLineRatios(self, cpar, handle, id, iter)
    implicit none
    class(comm_line_comp),                   intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration

    integer(i4b)    :: i, j, l, n, m, band, ierr
    real(dp)        :: A, b, mu, sigma, par, sigma_p, scale, w
    class(comm_map), pointer :: invN_amp, amp
    character(len=2) :: id_text
    
    band = self%ind2band(id)
    !if (band == self%ref_band) return

    
    ! Compute likelihood term
    w            = self%theta(id)%p%map(1,1)
    amp          => comm_map(data(band)%info)
    invN_amp     => comm_map(data(band)%info)
    amp%map      =  self%getBand(band)/w
    invN_amp%map = amp%map
    call data(band)%N%invN(invN_amp)     ! Inverse noise variance weighted amplitude map
    
!!$    call int2string(id, id_text)
!!$    call mask%writeFITS('co_mask'//id_text//'.fits')
!!$    call amp%writeFITS('co_amp'//id_text//'.fits')
!!$    call data(band)%res%writeFITS('co_res'//id_text//'.fits')
!!$    call data(band)%N%invN_diag%writeFITS('co_invN'//id_text//'.fits')

    ! Reduce across processors
    if (associated(self%indmask(band)%p)) then
       A = sum(invN_amp%map * self%indmask(band)%p%map * amp%map)
       b = sum(invN_amp%map * self%indmask(band)%p%map * data(band)%res%map)
    else
       A = sum(invN_amp%map * amp%map)
       b = sum(invN_amp%map * data(band)%res%map)
    end if
    call mpi_allreduce(MPI_IN_PLACE, A, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
    
    call amp%dealloc(); deallocate(amp)
    call invN_amp%dealloc(); deallocate(invN_amp)
    
    ! Compute new line ratio; just root processor
    if (self%x%info%myid == 0) then

       if (A > 0.d0) then
          mu    = b / A
          sigma = sqrt(1.d0 / A)
       else if (self%p_gauss(2,id) > 0.d0) then
          mu    = 0.d0
          sigma = 1.d30
       else
          mu    = self%p_uni(1,id) + (self%p_uni(2,id)-self%p_uni(1,id))*rand_uni(handle)
          sigma = 0.d0
       end if

       ! Add prior
       if (self%p_gauss(2,id) > 0.d0) then
          sigma_p = self%p_gauss(2,id) !/ sqrt(real(npix_reg,dp))
          mu      = (mu*sigma_p**2 + self%p_gauss(1,id) * sigma**2) / (sigma_p**2 + sigma**2)
          sigma   = sqrt(sigma**2 * sigma_p**2 / (sigma**2 + sigma_p**2))
       end if

       ! Draw sample
       par = -1.d300
       if (trim(self%operation) == 'optimize') then
          if (mu < self%p_uni(1,id)) then
             par = self%p_uni(1,id)
          else if (mu > self%p_uni(2,id)) then
             par = self%p_uni(2,id)
          else
             par = mu
          end if
       else
          do while (par < self%p_uni(1,id) .or. par > self%p_uni(2,id))
             if (mu < self%p_uni(1,id)) then
                par = rand_trunc_gauss(handle, mu, self%p_uni(1,id), sigma)
             else if (mu > self%p_uni(2,id)) then
                par = 2.d0*mu-rand_trunc_gauss(handle, mu, 2.d0*mu-self%p_uni(2,id), sigma)
             else
                par = mu + sigma * rand_gauss(handle)
             end if
          end do
       end if
       
       if (band == self%ref_band) then
          write(*,*) '  Line ratio i = ', id, ' = ', par, ' (at refband)'
       else
          write(*,*) '  Line ratio i = ', id, ' = ', par
       end if
    end if
    
    ! Distribute new relative line ratio, and update
    call mpi_bcast(par, 1, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)

    if (band == self%ref_band) then
       self%x%map = self%x%map * par  ! Rescale amplitude map, but leave mixing matrix
       self%x%alm = self%x%alm * par  
       do i = 1, self%npar            ! Rescale line ratios at other frequencies
          if (self%ind2band(i) == self%ref_band) cycle
          self%theta(i)%p%map = self%theta(i)%p%map / par
          if (self%lmax_ind >= 0) then
             self%theta(i)%p%alm(:,1) = self%theta(i)%p%alm(:,1) / par
          end if
       end do
    else
       self%theta(id)%p%map = par
       if (self%lmax_ind >= 0) then
          self%theta(id)%p%alm(:,1) = par * sqrt(4.d0*pi)
       end if
    end if
    call self%updateMixmat()

  end subroutine sampleLineRatios
  
end module comm_line_comp_mod
