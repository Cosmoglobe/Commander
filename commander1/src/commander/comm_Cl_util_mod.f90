module comm_Cl_util_mod
  use healpix_types
  use rngmod
  use comm_utils
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************


  integer(i4b)                                  :: num_Cl_bin
  integer(i4b),     allocatable, dimension(:,:) :: Cl_bins
  character(len=1), allocatable, dimension(:,:) :: Cl_bin_stat
 
  integer(i4b),               private :: nspec, lmax
  integer(i4b), parameter,    private :: MAX_NUM_BIN = 10000

contains

  subroutine initialize_Cl_util_mod(paramfile)
    implicit none

    character(len=128),                   intent(in) :: paramfile

    integer(i4b)       :: i, j, k, l, ind, upper, lower, unit, mcmc_lmin(6)
    logical(lgt)       :: polarization, cl_binning, ok
    character(len=2)   :: label
    character(len=128) :: filename, cl_prior_filename, t1
    character(len=512) :: line
    character(len=1), allocatable, dimension(:) :: stat
    real(dp), allocatable, dimension(:) :: input_vals

    unit = getlun()
    call get_parameter(paramfile, 'POLARIZATION',       par_lgt=polarization)
    call get_parameter(paramfile, 'LMAX',               par_int=lmax)

    if (polarization) then
       nspec = 6
    else
       nspec = 1
    end if

    ! Read binning information
    call get_parameter(paramfile, 'CL_BINNING', par_lgt=cl_binning)     
    allocate(Cl_bins(MAX_NUM_BIN,2), Cl_bin_stat(MAX_NUM_BIN,nspec), stat(nspec))
    Cl_bin_stat = '0'

    ! Add temperature monopole and dipole by default
    num_Cl_bin = 2
    Cl_bins(1,:) = [0,0]
    Cl_bins(2,:) = [1,1]

    if (cl_binning) then
       
       call get_parameter(paramfile, 'BINFILE',   par_string=filename)     

       open(unit,file=trim(filename))
       do while (.true.)
          read(unit,'(a)',end=813) line
          line = trim(line)
          if (line(1:1) == '#') cycle
          read(line,*) lower, upper, stat
          if (lower <= lmax .and. upper <= lmax .and. lower > 1 .and. upper > 1) then
             where (stat == 'H')
                stat = 'M'
             end where
             num_Cl_bin                  = num_Cl_bin + 1
             Cl_bins(num_Cl_bin,:)       = [lower,upper]
             Cl_bin_stat(num_Cl_bin,:)   = stat
          end if
       end do
813    close(unit)

       if (num_Cl_bin == 0) then
          write(*,*) 'Error: Cl bin file does not contain any valid entries'
          stop
       end if

    else
       
       do l = 2, lmax
          num_Cl_bin = num_Cl_bin + 1
          Cl_bins(num_Cl_bin,:) = l
          Cl_bin_stat(num_Cl_bin,:) = 'S'
       end do

    end if


    ! Read in overrides from the parameter file
    call get_parameter(paramfile, 'SAMPLE_TT_MODES', par_lgt=ok)
    if (.not. ok) Cl_bin_stat(:,1) = '0'
    if (polarization) then
       call get_parameter(paramfile, 'SAMPLE_TE_MODES', par_lgt=ok)
       if (.not. ok) Cl_bin_stat(:,2) = '0'
       call get_parameter(paramfile, 'SAMPLE_TB_MODES', par_lgt=ok)
       if (.not. ok) Cl_bin_stat(:,3) = '0'
       call get_parameter(paramfile, 'SAMPLE_EE_MODES', par_lgt=ok)
       if (.not. ok) Cl_bin_stat(:,4) = '0'
       call get_parameter(paramfile, 'SAMPLE_EB_MODES', par_lgt=ok)
       if (.not. ok) Cl_bin_stat(:,5) = '0'
       call get_parameter(paramfile, 'SAMPLE_BB_MODES', par_lgt=ok)
       if (.not. ok) Cl_bin_stat(:,6) = '0'
    end if

    do i = 1, num_Cl_bin
       do j = 1, nspec
          if (Cl_bin_stat(i,j) == 'M' .and. mcmc_lmin(j) > Cl_bins(i,2)) Cl_bin_stat(i,j) = 'S'
       end do
    end do

  end subroutine initialize_Cl_util_mod

  subroutine cleanup_Cl_util_mod
    implicit none

    if (allocated(Cl_bins))     deallocate(Cl_bins)
    if (allocated(Cl_bin_stat)) deallocate(Cl_bin_stat)

  end subroutine cleanup_Cl_util_mod
  
  
  subroutine initialize_ran_powspec(paramfile, rng_handle, cls_out)
    implicit none

    character(len=128),                 intent(in)    :: paramfile
    type(planck_rng),                   intent(inout) :: rng_handle
    real(dp),         dimension(0:,1:), intent(out)   :: cls_out

    real(dp)           :: val, disp_fac
    integer(i4b)       :: i, j, k, l, attempts, unit
    logical(lgt)       :: enforce_zero_cl, exist, condition_on_cl
    character(len=256) :: filename, clfile
    character(len=2048) :: line
    real(dp),     allocatable, dimension(:)     :: buffer, cl
    real(dp),     allocatable, dimension(:,:,:) :: cl_in

    unit = getlun()
    call get_parameter(paramfile, 'DISPERSION_FACTOR',  par_dp=disp_fac)
    call get_parameter(paramfile, 'ENFORCE_ZERO_CL',    par_lgt=enforce_zero_cl)
    call get_parameter(paramfile, 'IN_POWSPECFILE',     par_string=clfile)
    call get_parameter(paramfile, 'CONDITION_ON_CL',    par_lgt=condition_on_cl)

    cls_out = 0.d0
    if (enforce_zero_cl) return

    allocate(cl(nspec), cl_in(0:lmax, nspec, 2), buffer(2*nspec))
    
    ! Read input file
    cl_in = 0.d0
    inquire(file=trim(clfile),exist=exist)
    if (exist) then
       open(unit, file=trim(clfile), recl=1024)
       do while (.true.)
          read(unit,'(a)',end=1) line
          line = trim(line)
          if (line(1:1)=='#') cycle
          read(line,*) l, buffer
          if (l > lmax)       goto 1
          cl_in(l,:,1) = buffer(1::2)
          cl_in(l,:,2) = buffer(2::2)
       end do
1      close(unit)
    else
       write(*,*) 'Input power spectrum file does not exist. Exiting.'
       call mpi_finalize(i)
       stop
    end if
    if (condition_on_cl) cl_in(:,:,2) = 0.d0

    ! Set up Cl's
    do i = 1, num_Cl_bin
       if (all(Cl_bin_stat(i,:) == '0')) cycle
       attempts = 0
       do while (attempts < 1000)
          cl = 0.d0
          do l = Cl_bins(i,1), Cl_bins(i,2)
             do j = 1, nspec
                if (Cl_bin_stat(i,j) == 'C') then
                   cl(j) = cl(j) +  cl_in(l,j,1) 
                else if (Cl_bin_stat(i,j) == 'S' .or. Cl_bin_stat(i,j) == 'M') then
                   cl(j) = cl(j) + (cl_in(l,j,1) + disp_fac * rand_gauss(rng_handle) * cl_in(l,j,2))
                end if
             end do
          end do
          cl = cl / (Cl_bins(i,2)-Cl_bins(i,1)+1)
          if (is_posdef(cl)) then
             do l = Cl_bins(i,1), Cl_bins(i,2)
                cls_out(l,:) = cl
             end do
             exit
          else
             attempts = attempts+1
          end if
       end do

       if (attempts >= 1000) then
          write(*,*) 'Too many failed initialization attempts. Check input spectrum.'
          write(*,*) '   bin           = ', i
          write(*,*) '   l_low         = ', Cl_bins(i,1)
          write(*,*) '   l_high        = ', Cl_bins(i,2)
          call mpi_finalize(i)
          stop
       end if

    end do

    if (any(Cl_bin_stat(:,1) /= '0')) then
       ! Set monopole and dipole to a large value
       cls_out(0:1,1) = 1.d3
    end if

    deallocate(cl, cl_in, buffer)

  end subroutine initialize_ran_powspec


  subroutine bin_spectrum(cls)
    implicit none

    real(dp), dimension(0:,1:), intent(inout) :: cls 

    real(dp)     :: val, n
    integer(i4b) :: i, j, k, l

    do k = 1, num_Cl_bin
       do j = 1, nspec
          if (Cl_bin_stat(k,j) == '0' .or. Cl_bin_stat(k,j) == 'C') cycle
          val = 0.d0
          n   = 0.d0
          do l = Cl_bins(l,1), Cl_bins(l,2)
             val = val + (2*l+1) * cls(l,j)
             n   = n   + 2*l+1
          end do
          cls(Cl_bins(k,1):Cl_bins(k,2),j) = val/n
       end do
    end do
    
  end subroutine bin_spectrum

end module comm_Cl_util_mod
