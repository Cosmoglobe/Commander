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
module comm_template_comp_mod
  use math_tools
  use comm_param_mod
  use comm_comp_mod
  use comm_F_int_mod
  use comm_data_mod
  use pix_tools
  use comm_hdf_mod
  use comm_cr_utils
  use comm_cr_precond_mod
  use locate_mod
  use spline_1D_mod
  implicit none

  private
  public comm_template_comp, initTemplatePrecond, updateTemplatePrecond, applyTemplatePrecond, &
       & initialize_template_comps

  !**************************************************
  !            Compact object class
  !**************************************************
  type template

  end type template
  
  type, extends (comm_comp) :: comm_template_comp
     character(len=512)         :: outprefix
     real(dp)                   :: cg_scale
     integer(i4b)               :: band
     real(dp),  dimension(1,1)  :: x         ! Amplitude
     real(dp),  dimension(1,1)  :: x_buff    ! Amplitude
     real(dp),  dimension(2)    :: P         ! Gaussian prior on amplitude (mean,sigma)
     real(dp),  dimension(2)    :: P_cg      ! Gaussian prior on amplitude for CG (mean,sigma)
     class(comm_map), pointer   :: T    => null()   ! Template
     real(dp)                   :: T_scale          ! overall scaling amplitude for the template
     class(comm_map), pointer   :: mask => null()   ! Template mask
   contains
     procedure :: dumpFITS      => dumpTemplateToFITS
     procedure :: getBand       => evalTemplateBand
     procedure :: projectBand   => projectTemplateBand
     procedure :: S             => evalSED_template
     procedure :: initHDFComp   => initTemplateHDF
     procedure :: sampleSpecInd => sampleTempSpecInd
     procedure :: updateMixmat  => updateTempMixmat
     procedure :: update_F_int  => updateTempFInt
  end type comm_template_comp

  interface comm_template_comp
     procedure constructor_template
  end interface comm_template_comp

  type template_ptr
     class(comm_template_comp), pointer :: p => null()
  end type template_ptr
  
  integer(i4b) :: npre         =   0
  integer(i4b) :: comm_pre     =  -1
  integer(i4b) :: myid_pre     =  -1
  integer(i4b) :: numprocs_pre =  -1
  logical(lgt) :: recompute_template_precond = .true.
  
contains

  function constructor_template(cpar, id, id_abs, mu, rms, def, band, label, mapfile, maskfile) result(c)
    implicit none
    class(comm_params),        intent(in)           :: cpar
    integer(i4b),              intent(in)           :: id, id_abs
    real(dp),                  intent(in)           :: mu, rms, def
    integer(i4b),              intent(in), optional :: band
    character(len=*),          intent(in), optional :: label
    character(len=*),          intent(in), optional :: mapfile, maskfile
    class(comm_template_comp), pointer              :: c

    integer(i4b) :: i, j, n
    character(len=24), dimension(1000) :: comp_label
    character(len=512) :: dir

    ! General parameters
    allocate(c)


    ! Initialize general parameters
    c%class     = cpar%cs_class(id_abs)
    c%type      = cpar%cs_type(id_abs)
    c%label     = label !cpar%cs_label(id_abs)
    c%id        = id
    c%nmaps     = 1    ! Only used for CR book-keeping; must be 1 for templates
    c%outprefix = trim(cpar%cs_label(id_abs))
    c%init_from_HDF   = cpar%cs_initHDF(id_abs)
    c%cg_scale  = 1.d0
    c%output    = .true.
    c%myid      = cpar%myid_chain
    c%comm      = cpar%comm_chain
    c%numprocs  = cpar%numprocs_chain
    c%P         = [mu,rms]
    c%T_scale   = 1.d0
    npre                  = npre + 1
    comm_pre              = cpar%comm_chain
    myid_pre              = cpar%myid_chain
    numprocs_pre          = cpar%numprocs_chain

    if (c%myid == 0) then
       c%ncr = 1
       c%x   = def
    else
       c%ncr = 0
    end if

    if (present(mapfile)) then
       c%band      = band

       ! Read template and mask
       c%T => comm_map(data(band)%info, trim(mapfile))
       if (trim(maskfile) /= 'fullsky') then
          c%mask  => comm_map(data(band)%info, trim(maskfile))
          c%P_cg  =  c%P      ![mu,1.d-6]
       else
          c%P_cg  =  c%P      
       end if
    end if

    ! Set up CG sampling groups
    allocate(c%active_samp_group(cpar%cg_num_samp_groups))
    c%active_samp_group = .false.
    if (mu > 0.d0) then
       do i = 1, cpar%cg_num_samp_groups
          call get_tokens(cpar%cg_samp_group(i), ",", comp_label, n)
          do j = 1, n
             if (trim(c%label) == trim(comp_label(j))) then
                c%active_samp_group(i) = .true.
                if (n == 1) c%cg_unique_sampgroup = i ! Dedicated sampling group for this component
                exit
             end if
          end do
       end do
    end if

  end function constructor_template
  


  function initialize_template_comps(cpar, id, id_abs, n)
    implicit none
    type(comm_params),   intent(in)  :: cpar
    integer(i4b),        intent(in)  :: id, id_abs
    integer(i4b),        intent(out) :: n
    class(comm_template_comp), pointer   :: initialize_template_comps

    integer(i4b)        :: i, unit
    real(dp)            :: mu, rms, def
    character(len=1024) :: line, label, mapfile, maskfile
    class(comm_comp), pointer :: c => null()

    unit  = getlun()

    ! Find number of lines
    n = 0
    open(unit, file=trim(cpar%cs_SED_template(1,id_abs)), recl=1024)
    do while (.true.)
       read(unit,'(a)', end=1) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       read(line,*) label, mapfile, maskfile, mu, rms, def
       do i = 1, numband
          if (trim(label) == trim(data(i)%label)) exit
       end do
       if (i > numband) cycle

       if (n == 0) then
          initialize_template_comps => comm_template_comp(cpar, id+n, id_abs, mu, rms, def, &
               & i, trim(cpar%cs_label(id_abs))//'_'//trim(label), mapfile, maskfile)
       else
          c => comm_template_comp(cpar, id+n, id_abs, mu, rms, def, i, trim(cpar%cs_label(id_abs))//'_'//trim(label), mapfile, maskfile)
          call initialize_template_comps%addComp(c)
       end if
       n = n+1
    end do
1   close(unit)
  
  end function initialize_template_comps

  function evalSED_template(self, nu, band, pol, theta)
    class(comm_template_comp),  intent(in)           :: self
    real(dp),                   intent(in), optional :: nu
    integer(i4b),               intent(in), optional :: band
    integer(i4b),               intent(in), optional :: pol
    real(dp), dimension(1:),    intent(in), optional :: theta
    real(dp)                                        :: evalSED_template

    if (band == self%band) then
       evalSED_template = 1.d0
    else
       evalSED_template = 0.d0
    end if
    
  end function evalSED_template

  function evalTemplateBand(self, band, amp_in, pix, alm_out, det)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    integer(i4b),    dimension(:),   allocatable, intent(out), optional :: pix
    real(dp),        dimension(:,:),              intent(in),  optional :: amp_in
    logical(lgt),                                 intent(in),  optional :: alm_out
    integer(i4b),                                 intent(in),  optional :: det
    real(dp),        dimension(:,:), allocatable                        :: evalTemplateBand

    integer(i4b) :: i, j, p, q, ierr
    real(dp)     :: a

    if (.not. allocated(evalTemplateBand)) &
         & allocate(evalTemplateBand(0:data(band)%info%np-1,data(band)%info%nmaps))

    if (band /= self%band) then
       evalTemplateBand = 0.d0
       return
    end if

    if (self%myid == 0) then
       if (present(amp_in)) then
          a = amp_in(1,1)
       else
          a = self%x(1,1)
       end if
    end if
    call mpi_bcast(a, 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)

    evalTemplateBand = a * self%T%map
    
  end function evalTemplateBand
  
  ! Return component projected from map
  function projectTemplateBand(self, band, map, alm_in, det)
    implicit none
    class(comm_template_comp),                    intent(in)            :: self
    integer(i4b),                                 intent(in)            :: band
    class(comm_map),                              intent(in)            :: map
    logical(lgt),                                 intent(in), optional  :: alm_in
    integer(i4b),                                 intent(in), optional  :: det
    real(dp),        dimension(:,:), allocatable                        :: projectTemplateBand

    integer(i4b) :: i, j, q, p, ierr, d
    real(dp)     :: val, val2
    real(dp), allocatable, dimension(:,:) :: amp, amp2
    
    if (.not. allocated(projectTemplateBand)) allocate(projectTemplateBand(1,1))

    if (band /= self%band) then
       projectTemplateBand = 0.d0
       return
    end if

    val = sum(map%map * self%T%map)
    call mpi_reduce(val, val2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    if (self%myid == 0) projectTemplateBand = val2

  end function projectTemplateBand
  
  ! Dump current sample to HEALPix FITS file
  subroutine dumpTemplateToFITS(self, iter, chainfile, output_hdf, postfix, dir)
    class(comm_template_comp),               intent(inout)        :: self
    integer(i4b),                            intent(in)           :: iter
    type(hdf_file),                          intent(in)           :: chainfile
    logical(lgt),                            intent(in)           :: output_hdf
    character(len=*),                        intent(in)           :: postfix
    character(len=*),                        intent(in)           :: dir

    character(len=6) :: itext
    character(len=512) :: path

    if (.not. self%output) return

    if (self%myid == 0) write(*,*) '     Temp amp ', trim(adjustl(data(self%band)%label)), ' = ', self%x

    if (output_hdf .and. self%myid == 0) then
       call int2string(iter, itext)
       path = trim(adjustl(itext))//'/'//trim(self%outprefix)
       call create_hdf_group(chainfile, trim(adjustl(path)))
       path = trim(adjustl(path))//'/'//trim(adjustl(data(self%band)%label))
       call write_hdf(chainfile, trim(adjustl(path)), self%x)
       if(self%T_scale /= 1.d0 .and. self%myid == 0) then
         call write_hdf(chainfile, trim(adjustl(path)) // 'T_scale', self%T_scale)
       end if
    end if
    
  end subroutine dumpTemplateToFITS

  ! Dump current sample to HEALPix FITS file
  subroutine initTemplateHDF(self, cpar, hdffile, hdfpath)
    implicit none
    class(comm_template_comp),    intent(inout) :: self
    type(comm_params),         intent(in)    :: cpar    
    type(hdf_file),            intent(in)    :: hdffile
    character(len=*),          intent(in)    :: hdfpath

    character(len=512) :: path

    path = trim(adjustl(hdfpath))//'/'//trim(self%outprefix)//'/'//trim(adjustl(data(self%band)%label))
    call read_hdf(hdffile, path, self%x)

!!$    path = trim(adjustl(hdfpath))//trim(adjustl(self%label)) // '/'
!!$    if (self%myid == 0) then
!!$       call read_hdf(hdffile, trim(adjustl(path))//'/amp', self%x)
!!$       self%x = self%x/self%cg_scale
!!$       
!!$       allocate(theta(self%nsrc,self%nmaps,self%npar))
!!$       call read_hdf(hdffile, trim(adjustl(path))//'/specind', theta)
!!$       do i = 1, self%nsrc
!!$          do j = 1, self%nmaps
!!$             self%src(i)%theta(:,j) = theta(i,j,:) 
!!$          end do
!!$       end do
!!$       deallocate(theta)
!!$    end if
  end subroutine initTemplateHDF


  subroutine initTemplatePrecond(comm)
    implicit none
    integer(i4b),                intent(in) :: comm

    integer(i4b) :: i, i1, i2, j, j1, j2, k1, k2, q, l, m, n, p, p1, p2, n1, n2, myid, ierr, cnt
    real(dp)     :: t1, t2
    logical(lgt) :: skip
    class(comm_comp),          pointer :: c => null(), c1 => null(), c2 => null()
    class(comm_template_comp), pointer :: pt1 => null(), pt2 => null()
    class(comm_map),           pointer :: invN_T => null()
    real(dp),     allocatable, dimension(:,:) :: mat, mat2

    if (npre == 0) return
    if (allocated(P_cr%invM_temp)) return

    call mpi_comm_rank(comm, myid, ierr)
        
    ! Build frequency-dependent part of preconditioner
    call wall_time(t1)
    allocate(mat(npre,npre), mat2(npre,npre))
    mat = 0.d0
    i1  = 0
    c1 => compList
    do while (associated(c1))
       skip = .true.
       select type (c1)
       class is (comm_template_comp)
          pt1  => c1
          if (c1%band /= 0) skip = .false.
       end select
       if (skip) then
          c1 => c1%nextComp()
          cycle
       end if
       i1 = i1+1

       invN_T => comm_map(pt1%T)
       call data(pt1%band)%N%invN(invN_T) ! Multiply with invN

       i2 = 0
       c2 => compList
       do while (associated(c2))
          skip = .true.
          select type (c2)
          class is (comm_template_comp)
             pt2 => c2
             if (pt2%band == pt1%band) skip = .false.
          end select
          if (skip) then
             c2 => c2%nextComp()
             cycle
          end if
          i2 = i2+1
          !if (i2 < i1) cycle

          mat(i1,i2) = sum(invN_T%map * pt2%T%map)
          !mat(i2,i1) = mat(i1,i2)
          c2 => c2%nextComp()
       end do
       c1 => c1%nextComp()
    end do

    ! Collect contributions from all cores
    call mpi_reduce(mat, mat2, size(mat2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

    ! Invert matrix to finalize preconditioner
    allocate(P_cr%invM_temp(1,1))
    if (myid == 0) then
       ! Multiply with sqrtS from left side
       i1  = 0
       c1 => compList
       do while (associated(c1))
          skip = .true.
          select type (c1)
          class is (comm_template_comp)
             pt1  => c1
             if (c1%band /= 0) skip = .false.
          end select
          if (skip) then
             c1 => c1%nextComp()
             cycle
          end if
          i1         = i1+1
          mat2(i1,:) = mat2(i1,:) * pt1%P_cg(2)
          c1 => c1%nextComp()
       end do
       ! Multiply with sqrtS from right side
       i1  = 0
       c1 => compList
       do while (associated(c1))
          skip = .true.
          select type (c1)
          class is (comm_template_comp)
             pt1  => c1
             if (c1%band /= 0) skip = .false.
          end select
          if (skip) then
             c1 => c1%nextComp()
             cycle
          end if
          i1         = i1+1
          mat2(:,i1) = mat2(:,i1) * pt1%P_cg(2)
          c1 => c1%nextComp()
       end do
       ! Add unity
       do i1 = 1, npre
          mat2(i1,i1) = mat2(i1,i1) + 1.d0
       end do
       ! Invert matrix to finalize preconditioner
       call invert_matrix_with_mask(mat2)
       allocate(P_cr%invM_temp(1,1)%M(npre,npre))
       P_cr%invM_temp(1,1)%M = mat2
    end if

    deallocate(mat,mat2)
    
  end subroutine initTemplatePrecond

  subroutine updateTemplatePrecond(samp_group)
    implicit none
    integer(i4b),          intent(in) :: samp_group

    ! Placeholder for now; already fully initialized
    if (npre == 0) return
       
  end subroutine updateTemplatePrecond


  subroutine applyTemplatePrecond(x)
    implicit none
    real(dp),           dimension(:), intent(inout) :: x

    integer(i4b)              :: i, j, k, l, m, nmaps
    logical(lgt)              :: skip
    real(dp),                  allocatable, dimension(:,:) :: amp
    real(dp),                  allocatable, dimension(:)   :: y
    class(comm_comp),          pointer                     :: c => null() 
    class(comm_template_comp), pointer                     :: pt => null()

    return
    if (npre == 0 .or. myid_pre /= 0) return

    ! Reformat linear array into y(npre) structure
    allocate(y(npre))
    y = 0.d0
    l = 1
    c => compList
    do while (associated(c))
       skip = .true.
       select type (c)
       class is (comm_template_comp)
          pt => c
          if (pt%band /= 0) skip = .false.
       end select
       if (skip) then
          c => c%nextComp()
          cycle
       end if
       call cr_extract_comp(pt%id, x, amp)
       y(l) = amp(0,1)
       l  = l + 1
       c => c%nextComp()
       deallocate(amp)
    end do

    ! Multiply with preconditioner
    y = matmul(P_cr%invM_temp(1,1)%M, y)

    ! Reformat y(npre) structure back into linear array
    l = 1
    c => compList
    do while (associated(c))
       skip = .true.
       select type (c)
       class is (comm_template_comp)
          pt => c
          if (pt%band /= 0) skip = .false.
       end select
       if (skip) then
          c => c%nextComp()
          cycle
       end if
       allocate(amp(0:0,1:1))
       amp(0,1) = y(l)
       call cr_insert_comp(pt%id, .false., amp, x)
       l = l + 1
       c => c%nextComp()
       deallocate(amp)
    end do

    deallocate(y)

  end subroutine applyTemplatePrecond

  ! Sample spectral parameters
  subroutine sampleTempSpecInd(self, cpar, handle, id, iter)
    implicit none
    class(comm_template_comp),               intent(inout)        :: self
    type(comm_params),                       intent(in)           :: cpar
    type(planck_rng),                        intent(inout)        :: handle
    integer(i4b),                            intent(in)           :: id
    integer(i4b),                            intent(in)           :: iter    !Gibbs iteration
  end subroutine sampleTempSpecInd


  subroutine updateTempMixmat(self, theta, beta, band, df, par)
    implicit none
    class(comm_template_comp),         intent(inout)           :: self
    class(comm_map), dimension(:),     intent(in),    optional :: theta
    real(dp),        dimension(:,:,:), intent(in),    optional :: beta  ! (npar,nmaps,nsrc)
    integer(i4b),                      intent(in),    optional :: band
    class(map_ptr), dimension(:),      intent(inout), optional :: df
    integer(i4b),                      intent(in),    optional :: par


  end subroutine updateTempMixmat

  subroutine updateTempFInt(self, band)
    implicit none
    class(comm_template_comp), intent(inout)          :: self
    integer(i4b),              intent(in),   optional :: band

  end subroutine updateTempFInt

  
end module comm_template_comp_mod
