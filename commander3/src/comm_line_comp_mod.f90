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
     real(dp),           allocatable, dimension(:) :: line2RJ      ! This is the actual free parameter
     integer(i4b),       allocatable, dimension(:) :: ind2band
     integer(i4b),       allocatable, dimension(:) :: ind2det
     character(len=128), allocatable, dimension(:) :: freq_label
     character(len=128), allocatable, dimension(:) :: det_label
   contains
     procedure :: S                 => evalSED_line
     procedure :: sampleSpecInd     => sampleLineRatios
     procedure :: updateMixmat      => updateLineMixmat
     procedure :: update_F_int      => updateLineFInt
     procedure :: dumpFITS          => dumpLineToFITS
     procedure :: read_line_template
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
    type(comm_params),     intent(in) :: cpar
    integer(i4b),          intent(in) :: id, id_abs
    class(comm_line_comp), pointer   :: c

    integer(i4b) :: i, j, k
    
    ! General parameters
    allocate(c)
    c%npar = 0 !temporary value so that lmax_ind is correcty set (to 0) in initDiffuse
    allocate(c%poltype(1))
    c%poltype  = 1

    ! Read line template file
    call c%read_line_template(trim(cpar%cs_SED_template(1,id_abs)))

    allocate(c%lmax_ind_pol(3,c%npar))
    c%lmax_ind_pol = 0 !always fullsky (lmax=0) for line component
    if (allocated(c%lmax_ind_mix)) deallocate(c%lmax_ind_mix)
    allocate(c%lmax_ind_mix(3,c%npar))
    c%lmax_ind_mix = 0 !always fullsky (lmax=0) for line component
    allocate(c%pol_pixreg_type(3,c%npar))
    c%pol_pixreg_type = 0

    call c%initDiffuse(cpar, id, id_abs)

    ! Identify active parameters
    c%F_null = .true.
    do i = 1, c%npar
       j = c%ind2band(i)
       k = c%ind2det(i)
       if (j > 0) c%F_null(j,k) = .false.
    end do
    
    ! Update mixing matrix
    call c%updateMixmat
    
  end function constructor_line

  subroutine read_line_template(self, filename)
    implicit none
    class(comm_line_comp),                         intent(inout) :: self
    character(len=*),                              intent(in)    :: filename
    
    integer(i4b)        :: i, j, k, n, unit
    character(len=1024) :: line

    unit  = getlun()

    ! Find number of entries
    n = 0
    open(unit, file=trim(filename), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       n = n+1
    end do
1   close(unit)
    self%npar = n

    ! Allocate data structures
    allocate(self%theta_def(n), self%p_gauss(2,n), self%p_uni(2,n))
    if (allocated(self%poltype)) deallocate(self%poltype)
    allocate(self%poltype(n), self%indlabel(n), self%line2RJ(n), self%ind2band(n))
    allocate(self%freq_label(n), self%det_label(n), self%ind2det(n))

    ! Read actual data
    open(unit, file=trim(filename), recl=1024)
    j = 0
    do while (.true.)
       read(unit,'(a)', end=2) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       j = j+1
       read(line,*) self%freq_label(j), self%det_label(j), self%theta_def(j), self%p_gauss(1,j), self%p_gauss(2,j), self%poltype(j)
       self%p_uni(1,j)  = -1.d30
       self%p_uni(2,j)  =  1.d30
       self%indlabel(j) = 'na'
       self%ind2band(j) = 0
       self%ind2det(j)  = 0
       self%line2RJ(j)  = self%theta_def(j)

       ! Identify correct frequency and detector ID; default to 0 if they doesn't exist
       do i = 1, numband
          if (trim(self%freq_label(j)) == trim(data(i)%label)) then
             self%ind2band(j)  = i
             if (trim(data(i)%tod_type) == 'none' .or. trim(self%freq_label(j)) == 'all') cycle
             do k = 1, data(i)%tod%ndet
                if (trim(self%det_label(j)) == trim(data(i)%tod%label(k))) then
                   self%ind2det(j) = k
                   exit
                end if
             end do
             exit
          end if
       end do
       
    end do
2   close(unit)

  end subroutine read_line_template

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

    write(*,*) 'Warning: Should never be in evalSED_line!'
    evalSED_line = 0.d0

!!$    if (band == self%ref_band) then
!!$       evalSED_line = 1.d0
!!$    else 
!!$       do i = 1, self%npar
!!$          if (band == self%ind2band(i)) exit
!!$       end do
!!$       if (i > self%npar) then
!!$          evalSED_line = 0.d0
!!$       else
!!$          evalSED_line = theta(i) * self%line2RJ(i) / self%line2RJ_ref
!!$       end if
!!$    end if

  end function evalSED_line
  
  ! Evaluate amplitude map in brightness temperature at reference frequency
  subroutine updateLineMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_line_comp),                     intent(inout)           :: self
    class(comm_map),           dimension(:),   intent(in),    optional :: theta
    real(dp),  dimension(:,:,:),               intent(in),    optional :: beta  ! Not used here
    integer(i4b),                              intent(in),    optional :: band
    class(map_ptr), dimension(:),              intent(inout), optional :: df    ! Derivative of mixmat with respect to parameter par; for Jeffreys prior
    integer(i4b),                              intent(in),    optional :: par   ! Parameter ID for derivative

    integer(i4b) :: i, j, k, l

    ! Ccompute F_mean (= fullsky mixing matrix) for each band
    self%F_mean = 0.d0
    do i = 1, self%npar
       j = self%ind2band(i)
       k = self%ind2det(i)
       if (j > 0) then
          do l = 1, self%nmaps
             self%F_mean(j,k,l) = self%line2RJ(i) * data(j)%gain * self%cg_scale(l)
          end do
       end if
    end do

  end subroutine updateLineMixmat

  subroutine updateLineFInt(self, band)
    implicit none
    class(comm_line_comp), intent(inout)          :: self
    integer(i4b),             intent(in),   optional :: band
    ! Do nothing
  end subroutine updateLineFInt
  
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

    if (self%ind2det(id) > 0) then
       write(*,*) 'Warning: Line ratio sampling not yet supported: Freq, det = ', trim(self%freq_label(id)), trim(self%det_label(id))
       return
    end if
        
!!$    band = self%ind2band(id)
!!$    
!!$    ! Compute likelihood term
!!$    w            = self%theta(id)%p%map(1,1)
!!$    amp          => comm_map(data(band)%info)
!!$    invN_amp     => comm_map(data(band)%info)
!!$    amp%map      =  self%getBand(band)/w
!!$    invN_amp%map = amp%map
!!$    call data(band)%N%invN(invN_amp)     ! Inverse noise variance weighted amplitude map
!!$    
!!$    call int2string(id, id_text)
!!$    call mask%writeFITS('co_mask'//id_text//'.fits')
!!$    call amp%writeFITS('co_amp'//id_text//'.fits')
!!$    call data(band)%res%writeFITS('co_res'//id_text//'.fits')
!!$    call data(band)%N%invN_diag%writeFITS('co_invN'//id_text//'.fits')
!!$
!!$    ! Reduce across processors
!!$    if (associated(self%indmask(band)%p)) then
!!$       A = sum(invN_amp%map * self%indmask(band)%p%map * amp%map)
!!$       b = sum(invN_amp%map * self%indmask(band)%p%map * data(band)%res%map)
!!$    else
!!$       A = sum(invN_amp%map * amp%map)
!!$       b = sum(invN_amp%map * data(band)%res%map)
!!$    end if
!!$    call mpi_allreduce(MPI_IN_PLACE, A, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
!!$    call mpi_allreduce(MPI_IN_PLACE, b, 1, MPI_DOUBLE_PRECISION, MPI_SUM, self%x%info%comm, ierr)
!!$    
!!$    call amp%dealloc(); deallocate(amp)
!!$    call invN_amp%dealloc(); deallocate(invN_amp)
!!$    
!!$    ! Compute new line ratio; just root processor
!!$    if (self%x%info%myid == 0) then
!!$
!!$       if (A > 0.d0) then
!!$          mu    = b / A
!!$          sigma = sqrt(1.d0 / A)
!!$       else if (self%p_gauss(2,id) > 0.d0) then
!!$          mu    = 0.d0
!!$          sigma = 1.d30
!!$       else
!!$          mu    = self%p_uni(1,id) + (self%p_uni(2,id)-self%p_uni(1,id))*rand_uni(handle)
!!$          sigma = 0.d0
!!$       end if
!!$
!!$       ! Add prior
!!$       if (self%p_gauss(2,id) > 0.d0) then
!!$          sigma_p = self%p_gauss(2,id) !/ sqrt(real(npix_reg,dp))
!!$          mu      = (mu*sigma_p**2 + self%p_gauss(1,id) * sigma**2) / (sigma_p**2 + sigma**2)
!!$          sigma   = sqrt(sigma**2 * sigma_p**2 / (sigma**2 + sigma_p**2))
!!$       end if
!!$
!!$       ! Draw sample
!!$       par = -1.d300
!!$       if (trim(self%operation) == 'optimize') then
!!$          if (mu < self%p_uni(1,id)) then
!!$             par = self%p_uni(1,id)
!!$          else if (mu > self%p_uni(2,id)) then
!!$             par = self%p_uni(2,id)
!!$          else
!!$             par = mu
!!$          end if
!!$       else
!!$          do while (par < self%p_uni(1,id) .or. par > self%p_uni(2,id))
!!$             if (mu < self%p_uni(1,id)) then
!!$                par = rand_trunc_gauss(handle, mu, self%p_uni(1,id), sigma)
!!$             else if (mu > self%p_uni(2,id)) then
!!$                par = 2.d0*mu-rand_trunc_gauss(handle, mu, 2.d0*mu-self%p_uni(2,id), sigma)
!!$             else
!!$                par = mu + sigma * rand_gauss(handle)
!!$             end if
!!$          end do
!!$       end if
!!$       
!!$       if (band == self%ref_band) then
!!$          write(*,*) '  Line ratio i = ', id, ' = ', par, ' (at refband)'
!!$       else
!!$          write(*,*) '  Line ratio i = ', id, ' = ', par
!!$       end if
!!$    end if
!!$    
!!$    ! Distribute new relative line ratio, and update
!!$    call mpi_bcast(par, 1, MPI_DOUBLE_PRECISION, 0, self%x%info%comm, ierr)
!!$
!!$    if (band == self%ref_band) then
!!$       self%x%map = self%x%map * par  ! Rescale amplitude map, but leave mixing matrix
!!$       self%x%alm = self%x%alm * par  
!!$       do i = 1, self%npar            ! Rescale line ratios at other frequencies
!!$          if (self%ind2band(i) == self%ref_band) cycle
!!$          self%theta(i)%p%map = self%theta(i)%p%map / par
!!$          if (self%lmax_ind >= 0) then
!!$             self%theta(i)%p%alm(:,1) = self%theta(i)%p%alm(:,1) / par
!!$          end if
!!$       end do
!!$    else
!!$       self%theta(id)%p%map = par
!!$       if (self%lmax_ind >= 0) then
!!$          self%theta(id)%p%alm(:,1) = par * sqrt(4.d0*pi)
!!$       end if
!!$    end if
!!$    call self%updateMixmat()

  end subroutine sampleLineRatios

  ! Dump current sample to HEALPix FITS file
  subroutine dumpLineToFITS(self, iter, chainfile, output_hdf, postfix, dir)
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
    class(comm_line_comp),                   intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    integer(i4b)       :: i, l, j, k, m, ierr, unit, nnu, nuc
    integer(i4b)       :: p, p_min, p_max, npr, npol
    real(dp)           :: vals(10),theta(2)
    real(dp)           :: nu1, nu2, dlognu, nu, sed
    logical(lgt)       :: exist, first_call = .true.
    character(len=6)   :: itext
    character(len=512) :: filename, path
    class(comm_mapinfo), pointer :: info => null()
    class(comm_map), pointer :: map => null(), tp => null()
    real(dp), allocatable, dimension(:,:) :: sigma_l
    real(dp),     allocatable, dimension(:,:) :: dp_pixreg
    integer(i4b), allocatable, dimension(:,:) :: int_pixreg

    if (.not. self%output) return

    call update_status(status, "writeLine_1")

    ! Write amplitude
    map => comm_map(self%x)
    do i = 1, map%info%nmaps
       map%alm(:,i) = map%alm(:,i) * self%RJ2unit_(i) * self%cg_scale(i)  ! Output in requested units
    end do

    if (output_hdf) then
       call int2string(iter, itext)
       path = '/'//trim(adjustl(itext))//'/'//trim(adjustl(self%label))
       if (self%x%info%myid == 0) call create_hdf_group(chainfile, trim(adjustl(path)))
    end if
    
    filename = trim(self%label) // '_' // trim(postfix) // '.fits'
    call self%B_out%conv(trans=.false., map=map)
    call map%Y
    do i = 1, map%info%nmaps
       map%alm(:,i) = self%x%alm(:,i) * self%RJ2unit_(i) * self%cg_scale(i)  ! Replace convolved with original alms
    end do

    if (output_hdf) then
       call map%writeFITS(trim(dir)//'/'//trim(filename), &
            & hdffile=chainfile, hdfpath=trim(path)//'/amp_', output_hdf_map=.false.)
       !if we have set the overall scale parameter
       if (self%x_scale /= 1.d0 .and. self%myid == 0) then
          call write_hdf(chainfile, trim(path)//'/x_scale', self%x_scale)
       end if
    else
       call map%writeFITS(trim(dir)//'/'//trim(filename))
    end if
    call map%dealloc(); deallocate(map)
    call update_status(status, "writeLine_5")


!!$       ! Output Sampled SED's
!!$       if (output_hdf .and. allocated(self%SEDtab) .and. self%x%info%myid == 0) then
!!$         call write_hdf(chainfile, trim(path)//'/SED', self%SEDtab)
!!$         !!write the mbbTab SED for a range of frequencies from nu1 to nu2
!!$         ! this could maybe be updated for a more custom 'range' of frequencies in the future,
!!$         ! currently runs from 30GHz to the higher frequency in the table, and for 
!!$         ! 500 logarithmically spaced samples between those two frequencies (this could
!!$         ! also maybe be done more cleanly)
!!$         filename = trim(dir)// '/mbbTab_SED_' // trim(self%label) //'_'  // trim(postfix) // '.dat'
!!$         unit = getlun()
!!$         open(unit, file=trim(filename), status='replace')
!!$         write(unit,'(a)') '# nu[Hz]    SED[muK_RJ]'
!!$         nu1=30d0*1e9
!!$         nu2=self%SEDtab(2,self%ntab)
!!$         dlognu = (log(nu2) - log(nu1)) / 500d0
!!$         theta(1)=self%theta(1)%p%map(1,1)
!!$         theta(2)=self%theta(2)%p%map(1,1)
!!$         do nuc = 0, 500
!!$            nu  = exp(log(nu1) + dlognu*nuc)
!!$            sed = self%S(nu=nu, pol=1, theta=theta)
!!$            write(unit,'(2E20.10)') nu, sed
!!$         end do
!!$         close(unit)
!!$       end if
!!$       
!!$       ! Write mixing matrices
!!$       if (self%output_mixmat) then
!!$          do i = 1, numband
!!$             if (self%F_null(i,0)) cycle
!!$             filename = 'mixmat_' // trim(self%label) // '_' // trim(data(i)%label) // '_' // &
!!$                  & trim(postfix) // '.fits'
!!$             call self%F(i,0)%p%writeFITS(trim(dir)//'/'//trim(filename))
!!$          end do
!!$       end if
!!$       call update_status(status, "writeFITS_9")
!!$    end if

  end subroutine dumpLineToFITS
  
end module comm_line_comp_mod
 
