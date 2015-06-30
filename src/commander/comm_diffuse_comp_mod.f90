module comm_diffuse_comp_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_map_mod
  use comm_data_mod
  use comm_F_int_mod
  implicit none

  private
  public comm_diffuse_comp
  
  !**************************************************
  !            Diffuse component class
  !**************************************************
  type, abstract, extends (comm_comp) :: comm_diffuse_comp
     character(len=512) :: cltype
     integer(i4b)       :: nside, nmaps, nx, x0
     logical(lgt)       :: pol
     integer(i4b)       :: lmax_amp, lmax_ind, lpiv
     real(dp), allocatable, dimension(:,:) :: cls

     class(comm_map),               pointer     :: mask
     class(comm_map),               pointer     :: x      ! Spatial parameters
     class(map_ptr),  dimension(:), allocatable :: theta  ! Spectral parameters
     type(map_ptr),   dimension(:), allocatable :: F      ! Mixing matrix
     type(F_int_ptr), dimension(:), allocatable :: F_int  ! SED integrator
   contains
     procedure :: initDiffuse
     procedure :: updateMixmat
!!$     procedure :: a        => evalDiffuseAmp
!!$     procedure :: F        => evalDiffuseMixmat
!!$     procedure :: sim      => simDiffuseComp
!!$     procedure :: dumpHDF  => dumpDiffuseToHDF
     procedure :: dumpFITS => dumpDiffuseToFITS
  end type comm_diffuse_comp

contains

  subroutine initDiffuse(self, cpar, id)
    implicit none
    class(comm_diffuse_comp)            :: self
    type(comm_params),       intent(in) :: cpar
    integer(i4b),            intent(in) :: id

    integer(i4b) :: i
    type(comm_mapinfo), pointer :: info
    
    call self%initComp(cpar, id)

    ! Initialize variables specific to diffuse source type
    self%pol      = cpar%cs_polarization(id)
    self%nside    = cpar%cs_nside(id)
    self%lmax_amp = cpar%cs_lmax_amp(id)
    self%lmax_ind = cpar%cs_lmax_ind(id)
    self%nmaps    = 1; if (self%pol) self%nmaps = 3
    info          => comm_mapinfo(cpar%comm_chain, self%nside, self%lmax_amp, self%nmaps, self%pol)

    ! Initialize amplitude map
    if (trim(cpar%cs_input_amp(id)) == 'zero' .or. trim(cpar%cs_input_amp(id)) == 'none') then
       self%x => comm_map(info)
    else
       ! Read map from FITS file, and convert to alms
       write(*,*) 'hei'
       self%x => comm_map(info, cpar%cs_input_amp(id))
       !call self%x%YtW
    end if
    self%ncr = size(self%x%alm)

    ! Allocate mixing matrix
    allocate(self%F(numband))
    do i = 1, numband
       info      => comm_mapinfo(cpar%comm_chain, data(i)%info%nside, &
            & self%lmax_ind, data(i)%info%nmaps, data(i)%info%pol)
       self%F(i)%p => comm_map(info)
    end do

  end subroutine initDiffuse

  ! Evaluate amplitude map in brightness temperature at reference frequency
  subroutine updateMixmat(self, theta)
    implicit none
    class(comm_diffuse_comp),                  intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta

    integer(i4b) :: i, j, k, n
    real(dp),       allocatable, dimension(:) :: theta_p
    real(dp),       allocatable, dimension(:) :: nu, s
    class(map_ptr), allocatable, dimension(:) :: t
    
    ! Copy over alms from input structure, and compute pixel-space parameter maps
    if (present(theta)) then
       do i = 1, self%npar
          self%theta(i)%p%alm = theta(i)%alm
       end do
    end if

    ! Compute mixing matrix
    allocate(t(self%npar), theta_p(self%npar))
    do i = 1, numband

       write(*,*) trim(self%label), self%npar
       
       ! Compute spectral parameters at the correct resolution for this channel
       do j = 1, self%npar
          t(j)%p     => comm_map(self%F(j)%p%info)
          t(j)%p%alm =  self%theta(j)%p%alm
          !call t(j)%p%Y
       end do

       write(*,*) shape(t(1)%p%map), lbound(t(1)%p%map), ubound(t(1)%p%map)
       write(*,*) t(1)%p%map(700000,1)
       call mpi_finalize(j)
       stop

       ! Loop over all pixels, computing mixing matrix for each
       do j = 0, self%F(i)%p%info%np-1
          ! Temperature
          do k = 1, self%npar
             write(*,*) j, shape(t(k)%p%map)
             theta_p(k) = t(k)%p%map(j,1)
          end do
          self%F(i)%p%map(j,1) = self%F_int(i)%p%eval(theta_p)

          ! Polarization
          if (self%nmaps == 3) then
             ! Stokes Q
             if (all(self%poltype > 1)) then
                do k = 1, self%npar
                   if (self%poltype(k) > 1) then
                      theta_p(k) = t(k)%p%map(j,2)
                   else
                      theta_p(k) = t(k)%p%map(j,1)
                   end if
                end do
                self%F(i)%p%map(j,2) = self%F_int(i)%p%eval(theta_p)
             else
                self%F(i)%p%map(j,2) = self%F(i)%p%map(j,1)
             end if
       
             ! Stokes U
             if (all(self%poltype > 2)) then
                do k = 1, self%npar
                   if (self%poltype(k) > 2) then
                      theta_p(k) = t(k)%p%map(j,3)
                   else
                      theta_p(k) = t(k)%p%map(j,2)
                   end if
                end do
                self%F(i)%p%map(j,3) = self%F_int(i)%p%eval(theta_p)
             else
                self%F(i)%p%map(j,3) = self%F(i)%p%map(j,2)
             end if
          end if
          
       end do

    end do
    deallocate(t)

  end subroutine updateMixmat

  
!!$  ! Evaluate amplitude map in brightness temperature at reference frequency
!!$  function evalDiffuseAmp(self, nside, nmaps, pix, x_1D, x_2D)
!!$    class(comm_diffuse_comp),                  intent(in)           :: self
!!$    integer(i4b),                              intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:),   intent(in), optional :: pix
!!$    real(dp),                  dimension(:),   intent(in), optional :: x_1D
!!$    real(dp),                  dimension(:,:), intent(in), optional :: x_2D
!!$    real(dp),     allocatable, dimension(:,:)                       :: evalDiffuseAmp
!!$
!!$    if (present(x_1D)) then
!!$       ! Input array is complete packed data vector
!!$       
!!$    else if (present(x_2D)) then
!!$       ! Input array is alms(:,:)
!!$       
!!$    else
!!$       ! Use internal array
!!$       
!!$    end if
!!$    
!!$  end function evalDiffuseAmp
!!$
!!$  ! Evaluate amplitude map in brightness temperature at reference frequency
!!$  function evalDiffuseMixmat(self, nside, nmaps, pix)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    integer(i4b),                            intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:), intent(in), optional :: pix
!!$    real(dp),     allocatable, dimension(:,:)                     :: evalDiffuseMixmat
!!$  end function evalDiffuseMixmat
!!$
!!$  ! Generate simulated component
!!$  function simDiffuseComp(self, handle, nside, nmaps, pix)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    type(planck_rng),                        intent(inout)        :: handle
!!$    integer(i4b),                            intent(in)           :: nside, nmaps
!!$    integer(i4b),              dimension(:), intent(in), optional :: pix
!!$    real(dp),     allocatable, dimension(:,:)                     :: simDiffuseComp
!!$  end function simDiffuseComp
!!$  
!!$  ! Dump current sample to HDF chain files
!!$  subroutine dumpDiffuseToHDF(self, filename)
!!$    class(comm_diffuse_comp),                intent(in)           :: self
!!$    character(len=*),                        intent(in)           :: filename
!!$  end subroutine dumpDiffuseToHDF
!!$  

  ! Dump current sample to HEALPix FITS file
  subroutine dumpDiffuseToFITS(self, postfix, dir)
    implicit none
    character(len=*),                        intent(in)           :: postfix
    class(comm_diffuse_comp),                intent(in)           :: self
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i
    character(len=512) :: filename

    write(*,*) 'q1'
    ! Write amplitude
    filename = trim(self%label) // '_' // trim(postfix) // '.fits'
    call self%x%Y
    write(*,*) 'q2'
    call self%x%writeFITS(trim(dir)//'/'//trim(filename))
    write(*,*) 'q3'

    ! Write spectral index maps
    write(*,*) 'q4'
    do i = 1, self%npar
       filename = trim(self%label) // '_' // trim(self%indlabel(i)) // '_' // &
            & trim(postfix) // '.fits'
       call self%theta(i)%p%Y
       call self%theta(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
    end do
    write(*,*) 'q4'

    ! Write mixing matrices
    write(*,*) numband
    do i = 1, numband
       filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
            & trim(postfix) // '.fits'
       write(*,*) i, trim(filename)
       call self%F(i)%p%writeFITS(trim(dir)//'/'//trim(filename))
    end do
    write(*,*) 'q5'
        
    
  end subroutine dumpDiffuseToFITS
  
end module comm_diffuse_comp_mod
