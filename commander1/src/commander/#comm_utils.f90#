module comm_utils
  use healpix_types
  use alm_tools
  use rngmod
  use fitstools
  use head_fits
  use pix_tools
  use udgrade_nr
  use mpi_alm_tools
  use math_tools
  use iso_c_binding
  use comm_system_mod
  use sort_utils
  implicit none

  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2006, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

  !integer(i8b), external :: get_mem_use, get_max_mem_use
  !integer(i4b), external :: get_pid, open_atomic_file, nfork, popen, svn_revision, ishift, count_set_bits
  !logical(lgt), external :: is_dir, mkdir_c, mv_c, rm_c
  !integer(i4b), parameter :: SIGINT = 2, SIGSEGV = 11, SIGTERM = 15, SIGBUS = 7
  !real(c_double),     bind(C, name="nan")          :: nan
  !real(c_double),     bind(C, name="snan")         :: snan
  !real(c_double),     bind(C, name="infinity")     :: infinity
  !integer(c_int64_t), bind(C, name="fe_divbyzero") :: fe_divbyzero
  !integer(c_int64_t), bind(C, name="fe_inexact")   :: fe_inexact
  !integer(c_int64_t), bind(C, name="fe_nan")       :: fe_nan
  !integer(c_int64_t), bind(C, name="fe_overflow")  :: fe_overflow
  !integer(c_int64_t), bind(C, name="fe_underflow") :: fe_underflow

  interface compute_sigma_l
     module procedure compute_sigma_l_vec, compute_sigma_l_mat
  end interface

  real(dp), allocatable, dimension(:,:) :: plms

contains

  function getlun()
    implicit none
    integer(i4b) :: getlun
    logical(lgt) :: exists, isopen
    getlun = 9
    do
       getlun = getlun+1
       inquire(unit=getlun,exist=exists)
       if(exists) then
          inquire(unit=getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function getlun

  subroutine read_beam(beamfile, beam)
    implicit none

    character(len=128),                   intent(in)  :: beamfile
    real(dp),           dimension(0:,1:), intent(out) :: beam

    integer(i4b) :: l, lmax, nmaps
    real(dp)     :: sigma_sq
    real(dp),          allocatable, dimension(:,:)   :: inbeam
    character(len=80),              dimension(1:180) :: header
    
    lmax  = size(beam(:,1))-1
    nmaps = size(beam(0,:))

    ! Seem to remember there is something weird going on with the WMAP beams when reading only 
    ! one component at a time. Remove this wrapper if you feel more comfortable with that...
    allocate(inbeam(0:lmax,4))

!    call read_asctab(beamfile, inbeam, lmax, 4, header)
    call fits2cl(beamfile, inbeam, lmax, 4, header)

    beam(0:lmax,1:nmaps) = inbeam(0:lmax,1:nmaps)

    if (nmaps > 1) then
       if (sum(beam(:,2)) < 1.d0) beam(:,2) = beam(:,1)
       if (sum(beam(:,3)) < 1.d0) beam(:,3) = beam(:,2)
    end if

    deallocate(inbeam)

  end subroutine read_beam


!!$  subroutine multiply_with_beam(beam, alms)
!!$    implicit none
!!$
!!$    real(dp), dimension(0:,1:), intent(in)    :: beam
!!$    real(dp), dimension(1:,1:), intent(inout) :: alms
!!$
!!$    integer(i4b) :: j, l, m, lmax
!!$
!!$    lmax = size(beam(:,1))-1
!!$
!!$    j = 1
!!$    do l = 0, lmax
!!$       do m = -l, l
!!$          alms(j,:) = alms(j,:) * beam(l,:)
!!$          j = j+1
!!$       end do
!!$    end do
!!$
!!$  end subroutine multiply_with_beam


  subroutine convert_RMS_to_invN(mask, rmsmap, inv_N_covar)
    implicit none

    logical(lgt), dimension(0:,1:), intent(in)  :: mask
    real(dp),     dimension(0:,1:), intent(in)  :: rmsmap
    real(dp),     dimension(0:,1:), intent(out) :: inv_N_covar

    integer(i4b) :: i, j, map_size, nmaps

    map_size = size(mask(:,1))
    nmaps    = size(rmsmap(0,:))

    inv_N_covar = 0.d0
    do j = 1, nmaps
       do i = 0, map_size-1
          if (mask(i,j)) then
             inv_N_covar(i,j) = 1.d0 / rmsmap(i,j)**2
          end if
       end do
    end do

  end subroutine convert_RMS_to_invN



  subroutine convert_real_to_harmonic_space(map, alms)
    implicit none

    real(dp), dimension(0:,1:), intent(in)            :: map
    real(dp), dimension(1:,1:), intent(out), optional :: alms

    real(dp) :: t1, t2, zbounds(2)
    integer(i4b) :: lmax, nmaps, nside
    logical(lgt) :: single_core
    real(dp), allocatable, dimension(:,:) :: weights
    complex(dpc), allocatable, dimension(:,:,:)        :: alms_cmplx

    nside       = get_mpi_nside()
    single_core = (12*nside**2 == size(map,1))

    if (single_core) then

       nmaps = size(alms(1,:))
       lmax  = numcomp2lmax(size(alms(:,1)))
       allocate(alms_cmplx(nmaps,0:lmax,0:lmax), weights(1:2*nside, nmaps))
       weights = 1.d0
       zbounds = 0.d0
       if (nmaps == 1) then
          if (allocated(plms)) then
             call map2alm(nside, lmax, lmax, map(:,1), alms_cmplx, zbounds, weights, plms(:,1))
          else
             call map2alm(nside, lmax, lmax, map(:,1), alms_cmplx, zbounds, weights)
          end if
       else
          if (allocated(plms)) then
             call map2alm(nside, lmax, lmax, map,      alms_cmplx, zbounds, weights, plms)
          else
             call map2alm(nside, lmax, lmax, map,      alms_cmplx, zbounds, weights)
          end if
       end if
       call convert_complex_to_real_alms(alms_cmplx, alms)
       deallocate(alms_cmplx, weights)

    else

       if (present(alms)) then
          
          nmaps = size(alms(1,:))
          lmax  = numcomp2lmax(size(alms(:,1)))
          allocate(alms_cmplx(nmaps,0:lmax,0:lmax))
          
          call mpi_map2alm_dist(map, alms_cmplx)
          
          call convert_complex_to_real_alms(alms_cmplx, alms)
          
          deallocate(alms_cmplx)
          
       else 
          
          call mpi_map2alm_dist_slave(map)
          
       end if
       
    end if

  end subroutine convert_real_to_harmonic_space


  subroutine convert_harmonic_to_real_space(map, alms)
    implicit none

    real(dp), dimension(0:,1:), intent(out)           :: map
    real(dp), dimension(1:,1:), intent(in),  optional :: alms

    real(dp)     :: t1, t2
    integer(i4b) :: npix, lmax, nmaps, nside
    logical(lgt) :: single_core
    complex(dpc), allocatable, dimension(:,:,:)  :: alms_cmplx

    nside       = get_mpi_nside()
    single_core = (12*nside**2 == size(map,1))

    if (single_core) then

       lmax  = numcomp2lmax(size(alms(:,1)))
       nmaps = size(alms(1,:))
          
       allocate(alms_cmplx(nmaps,0:lmax,0:lmax))
       call convert_real_to_complex_alms(alms, alms_cmplx)
       if (nmaps == 1) then
          if (allocated(plms)) then
             call alm2map(nside, lmax, lmax, alms_cmplx, map(:,1), plms(:,1))
          else
             call alm2map(nside, lmax, lmax, alms_cmplx, map(:,1))
          end if
       else
          if (allocated(plms)) then
             call alm2map(nside, lmax, lmax, alms_cmplx, map, plms)
          else
             call alm2map(nside, lmax, lmax, alms_cmplx, map)
          end if
       end if
       deallocate(alms_cmplx)

    else

       if (present(alms)) then
          
          lmax  = numcomp2lmax(size(alms(:,1)))
          nmaps = size(alms(1,:))
          
          allocate(alms_cmplx(nmaps,0:lmax,0:lmax))
          
          call convert_real_to_complex_alms(alms, alms_cmplx)
          
          call cpu_time(t1)
          call mpi_alm2map_dist(alms_cmplx, map)
          call cpu_time(t2)
       
          deallocate(alms_cmplx)
          
       else
          
          call cpu_time(t1)
          call mpi_alm2map_dist_slave(map)
          call cpu_time(t2)
          
       end if
       
    end if

  end subroutine convert_harmonic_to_real_space


  subroutine convert_real_to_complex_alms(alms_real, alms_cmplx)
    implicit none

    real(dp),     dimension(1:,1:),    intent(in)  :: alms_real
    complex(dpc), dimension(1:,0:,0:), intent(out) :: alms_cmplx

    integer(i4b) :: ind, l, m, lmax
    real(dp)     :: one_over_sqrt_two 

    lmax              = size(alms_cmplx(1,:,0))-1
    one_over_sqrt_two = 1.d0 / sqrt(2.d0)

    alms_cmplx = cmplx(0.d0, 0.d0)
    do l = 0, lmax
       ind = lm2ind(l,0)

       alms_cmplx(:,l,0) = cmplx(alms_real(ind,:),0.d0)
       do m = 1, l
          alms_cmplx(:,l,m) = one_over_sqrt_two * cmplx(alms_real(ind+m,:), alms_real(ind-m,:))
       end do
    end do

  end subroutine convert_real_to_complex_alms

  subroutine convert_complex_to_real_alms(alms_cmplx, alms_real)
    implicit none

    complex(dpc), dimension(1:,0:,0:), intent(in)  :: alms_cmplx
    real(dp),     dimension(1:,1:),    intent(out) :: alms_real

    integer(i4b) :: ind, l, m, lmax
    real(dp)     :: sqrt_two

    lmax     = size(alms_cmplx(1,:,0))-1
    sqrt_two = sqrt(2.d0)

    alms_real = cmplx(0.d0, 0.d0)
    do l = 0, lmax
       ind = lm2ind(l,0)

       alms_real(ind,:) = real(alms_cmplx(:,l,0),dp)
       do m = 1, l
          alms_real(ind-m,:) = sqrt_two * imag(alms_cmplx(:,l,m))
          alms_real(ind+m,:) = sqrt_two * real(alms_cmplx(:,l,m))
       end do
    end do

  end subroutine convert_complex_to_real_alms



  function numcomp2lmax(numcomp)
    implicit none

    integer(i4b), intent(in)  :: numcomp
    integer(i4b)              :: numcomp2lmax

!    numcomp2lmax = (nint(sqrt(real(8*numcomp+1,dp)))-3)/2
    numcomp2lmax = nint(sqrt(real(numcomp,dp)))-1

  end function numcomp2lmax


  function lmax2numcomp(lmax)
    implicit none

    integer(i4b), intent(in)  :: lmax
    integer(i4b)              :: lmax2numcomp

!    lmax2numcomp = (lmax+1) * (lmax+2) / 2
    lmax2numcomp = (lmax+1)**2

  end function lmax2numcomp

  function lm2ind(l,m)
    implicit none

    integer(i4b), intent(in)  :: l, m
    integer(i4b)              :: lm2ind

!    lm2ind = (l-1)*l/2 + l + m + 1
    lm2ind = l**2 + l + m + 1

  end function lm2ind




  subroutine output_spectra_from_iteration(chain, iteration, chain_dir, cls, sigmas)
    implicit none

    integer(i4b),                   intent(in) :: chain, iteration
    character(len=128),             intent(in) :: chain_dir
    real(dp),     dimension(0:,1:), intent(in) :: cls, sigmas

    integer(i4b)       :: l, i, j, ind, lmax, nmaps, nspec, unit
    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename

    unit  = getlun()
    lmax  = size(cls(:,1))-1
    nspec = size(cls(0,:))

    call int2string(chain,     chain_text)
    call int2string(iteration, i_text)
    
    filename = trim(chain_dir) // '/cl_c' // chain_text // '_k' // i_text // '.dat'
    
    open(unit,file=trim(filename), recl=1024)
    do l = 2, lmax
       write(unit,*) l, real(cls(l,:),sp)
    end do
    write(unit,*) ''
    do l = 2, lmax
       write(unit,*) l, real(sigmas(l,:),sp)
    end do
    close(unit)

  end subroutine output_spectra_from_iteration


  subroutine output_maps_from_iteration(nside, lmax, chain, iteration, chain_dir, s, fwhm, &
       & pix, cmb)
    implicit none

    integer(i4b),               intent(in)             :: nside, lmax, chain, iteration
    character(len=*),           intent(in)             :: chain_dir
    real(dp), dimension(1:,1:), intent(in)             :: s
    real(dp),                   intent(in)             :: fwhm
    integer(i4b), dimension(:), intent(in),  optional  :: pix
    real(dp), dimension(:,:),   intent(out), optional  :: cmb


    integer(i4b)       :: i, j, l, m, npix, nmaps
    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename
    real(dp),     allocatable, dimension(:,:)    :: map, gb
    complex(dpc), allocatable, dimension(:,:,:)  :: alms

    npix  = 12*nside**2
    nmaps = size(s(1,:))

    allocate(alms(nmaps,0:lmax,0:lmax))
    allocate(map(0:npix-1,nmaps))
    allocate(gb(0:lmax,nmaps))
       
    call convert_real_to_complex_alms(s, alms)
    call gaussbeam(fwhm, lmax, gb)

    do j = 1, nmaps
       do l = 0, lmax
          alms(j,l,0:l) = alms(j,l,0:l) * gb(l,j)
       end do
    end do
    deallocate(gb)

    ! Write standard Healpix map
    !alms(1,0:1,:) = 0.d0
    if (nmaps == 1) then
       call alm2map(nside, lmax, lmax, alms, map(:,1))
    else
       call alm2map(nside, lmax, lmax, alms, map)
    end if
       
    call int2string(chain,     chain_text)
    call int2string(iteration, i_text)
    filename = trim(chain_dir) // '/cmb_Cl_c' // chain_text // '_k' // i_text // '.fits'
    call write_map(filename, map, unit='uK_cmb', comptype='CMB')

    if (present(pix)) then
       if (size(pix) > 0) then
          do i = 1, size(pix)
             cmb(i,:) = map(pix(i),:)
          end do
       end if
    end if

    deallocate(alms)
    deallocate(map)

  end subroutine output_maps_from_iteration

  subroutine output_alms_from_iteration(nside, lmax, chain, iteration, chain_dir, s)
    implicit none

    integer(i4b),               intent(in) :: nside, lmax, chain, iteration
    character(len=128),         intent(in) :: chain_dir
    real(dp), dimension(1:,1:), intent(in) :: s

    integer(i4b)       :: i, l, m, nmaps
    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename
    complex(dpc), allocatable, dimension(:,:,:)  :: alms

    nmaps = size(s(1,:))

    allocate(alms(nmaps,0:lmax,0:lmax))

    call convert_real_to_complex_alms(s, alms)

    ! Write alms
    call int2string(chain,     chain_text)
    call int2string(iteration, i_text)
    filename = trim(chain_dir) // '/alms_c' // chain_text // '_k' // i_text // '.fits'
    call write_s_lm(filename, alms)

    deallocate(alms)

  end subroutine output_alms_from_iteration

  subroutine output_chisq_from_iteration(chain, iteration, chain_dir, chisq)
    implicit none

    integer(i4b),               intent(in) :: chain, iteration
    character(len=128),         intent(in) :: chain_dir
    real(dp), dimension(0:,1:), intent(in) :: chisq

    character(len=4)   :: chain_text
    character(len=5)   :: i_text
    character(len=128) :: filename
       
    call int2string(chain,     chain_text)
    call int2string(iteration, i_text)
    filename = trim(chain_dir) // '/chisq_c' // chain_text // '_k' // i_text // '.fits'
    call write_map(filename, chisq, comptype='Chi^2')

  end subroutine output_chisq_from_iteration


  subroutine smooth_par_map(map, fwhm, prior)
    implicit none

    real(dp), dimension(0:,1:), intent(inout) :: map
    real(dp),                   intent(in)    :: fwhm
    real(dp), dimension(2),     intent(in)    :: prior

    integer(i4b) :: i, j, l, m, nside, npix, nmaps, lmax
    real(dp),     allocatable, dimension(:,:)   :: weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms

    if (fwhm <= 0.d0) return

    npix  = size(map,1)
    nmaps = size(map,2)
    nside = nint(sqrt(real(npix,dp)/12.d0))
    lmax  = 3*nside
    allocate(alms(1,0:lmax,0:lmax), weights(2*nside,1))
    weights = 1.d0
    
    do i = 1, nmaps
       if (all(map(:,i) == map(0,1))) cycle
       call map2alm(nside, lmax, lmax, map(:,i), alms, [0.d0, 0.d0], weights)
       do l = 0, lmax
          alms(1,l,0:l) = alms(1,l,0:l) * exp(-0.5d0 * l*(l+1)*(fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
       end do
       call alm2map(nside, lmax, lmax, alms, map(:,i))
    end do

    where (map < prior(1))
       map = prior(1)
    end where
    where (map > prior(2))
       map = prior(2)
    end where

    deallocate(alms, weights)

  end subroutine smooth_par_map


  subroutine compute_sigma_and_cls(s_i, cl_i, sigma_out, cls_out)
    implicit none

    real(dp), dimension(1:,1:),    intent(in)  :: s_i
    real(dp), dimension(0:,1:,1:), intent(in)  :: cl_i
    real(dp), dimension(0:,1:),    intent(out) :: sigma_out, cls_out

    integer(i4b) :: i, j, k, l, m, ind, nmaps, lmax, nspec
    real(dp) :: sigma

    lmax  = size(cl_i(:,1,1))-1
    nmaps = size(cl_i(0,:,1))
    nspec = nmaps*(nmaps+1)/2

    sigma_out = 0.d0
    cls_out   = 0.d0
    do i = 1, nmaps
       do k = i, nmaps
          ind = i*(1-i)/2 + (i-1)*nmaps + k

          j = 1
          do l = 0, lmax
             sigma = 0.d0
             do m = -l, l
                sigma = sigma + s_i(j,k) * s_i(j,i)
                j = j+1
             end do
             sigma_out(l,ind) = sigma / real(2*l+1,dp) * real(l*(l+1),dp) / (2.d0*pi)
             cls_out(l,ind)   = cl_i(l,i,k)            * real(l*(l+1),dp) / (2.d0*pi)
          end do
       end do
    end do

  end subroutine compute_sigma_and_cls


  subroutine compute_sigma_l_vec(s_i, cl_i)
    implicit none

    real(dp), dimension(1:,1:), intent(in)   :: s_i
    real(dp), dimension(0:,1:), intent(out)  :: cl_i

    integer(i4b) :: i, j, k, l, m, s, ind, nmaps, lmax, nspec
    real(dp) :: sigma

    lmax  = size(cl_i,1)-1
    nmaps = size(s_i,2)
    nspec = nmaps*(nmaps+1)/2

    cl_i = 0.d0
    s = 1
    do i = 1, nmaps
       do k = i, nmaps
          j = 1
          do l = 0, lmax
             sigma = 0.d0
             do m = -l, l
                sigma = sigma + s_i(j,k) * s_i(j,i)
                j = j+1
             end do
             cl_i(l,s) = sigma / real(2*l+1,dp) 
          end do
          s = s+1
       end do
    end do

  end subroutine compute_sigma_l_vec

  subroutine compute_sigma_l_mat(s_i, cl_i)
    implicit none

    real(dp), dimension(1:,1:),    intent(in)   :: s_i
    real(dp), dimension(0:,1:,1:), intent(out)  :: cl_i

    integer(i4b) :: i, j, k, l, m, s, ind, nmaps, lmax, nspec
    real(dp) :: sigma

    lmax  = size(cl_i,1)-1
    nmaps = size(s_i,2)
    nspec = nmaps*(nmaps+1)/2

    cl_i = 0.d0
    do i = 1, nmaps
       do k = i, nmaps
          j = 1
          do l = 0, lmax
             sigma = 0.d0
             do m = -l, l
                sigma = sigma + s_i(j,k) * s_i(j,i)
                j = j+1
             end do
             cl_i(l,i,k) = sigma / real(2*l+1,dp) 
             cl_i(l,k,i) = cl_i(l,i,k)
          end do
       end do
    end do

  end subroutine compute_sigma_l_mat


  ! Small utility for converting an integer to a string
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string



  ! *****************************************************************************************
  !
  ! Routine for reading one parameter from a parameter file
  !     Example usage:   call get_parameter(21, "parfile.txt", "NSIDE", par_int=nside)
  !
  ! *****************************************************************************************

!!$  subroutine get_parameter(parfile, parname, par_int, par_char, &
!!$       & par_string, par_sp, par_dp, par_lgt)
!!$    implicit none
!!$
!!$    character(len=*),  intent(in)  :: parfile, parname
!!$    integer(i4b),      intent(out), optional :: par_int
!!$    logical(lgt),      intent(out), optional :: par_lgt
!!$    character(len=1),  intent(out), optional :: par_char
!!$    character(len=*),  intent(out), optional :: par_string
!!$    real(sp),          intent(out), optional :: par_sp
!!$    real(dp),          intent(out), optional :: par_dp
!!$
!!$
!!$    integer(i4b)        :: i, unit
!!$    character(len=128)  :: string, variable, value
!!$    character(len=1)    :: equals
!!$
!!$    unit = getlun()
!!$    open(unit, file=trim(parfile))
!!$
!!$    do while (.true.)
!!$       read(unit,*,end=1) string
!!$
!!$       if (string(1:1)=='#') cycle
!!$
!!$       backspace(unit)
!!$       read(unit,*) variable, equals, value
!!$
!!$       if (trim(variable) == trim(parname)) then
!!$
!!$          if (present(par_int)) then
!!$             read(value,*) par_int
!!$          else if (present(par_char)) then
!!$             read(value,*) par_char
!!$          else if (present(par_string)) then
!!$             read(value,'(a)') par_string
!!$          else if (present(par_sp)) then
!!$             read(value,*) par_sp
!!$          else if (present(par_dp)) then
!!$             read(value,*) par_dp
!!$          else if (present(par_lgt)) then
!!$             read(value,*) par_lgt
!!$          end if
!!$
!!$          close(unit)
!!$          return
!!$
!!$       end if
!!$
!!$    end do
!!$
!!$1   write(*,*) 'GET_PARAMETER:    Critical error -- parameter not found'
!!$    write(*,*) 'GET_PARAMETER:       Parameter file = ', trim(parfile)
!!$    write(*,*) 'GET_PARAMETER:       Parameter name = ', trim(parname)
!!$
!!$    close(unit)
!!$    stop
!!$
!!$  end subroutine get_parameter


  subroutine get_parameter(parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
    integer(i4b)               :: unit

    unit = getlun()
    found = .false.
    call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
    if(found) then
       if(present(par_present)) par_present = .true.
    else
       call get_parameter_parfile(unit, parfile, parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
    end if
  end subroutine

  subroutine get_parameter_parfile(unit, parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    integer(i4b)               :: unit
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i
    character(len=512)         :: key, value, filenames(maxdepth), line
    character(len=1)           :: equals

    depth = 1
    units(depth) = getlun()
    !write(*,*) "Entering file " // trim(parfile)
    filenames(depth) = parfile
    open(units(depth),file=trim(parfile),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),*,end=1) key, value
             !write(*,*) "Entering file " // trim(value)
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          call parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
          if(found) then
             ! Match found, so clean up and return.
             do i = depth, 1, -1; close(units(i)); end do
             if(present(par_present)) par_present = .true.
             return
          end if
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do

    ! ===== Error handling section ======
    ! Case 1: Failed to find matching parameter in file
    if(present(par_present)) then
       par_present = .false.
       return
    else
       write(*,*) "get_parameter: Fatal error: Cannot find parameter " // trim(parname) // &
            & " in parameter file " // trim(parfile) // " or included files."
       if(present(desc)) write(*,*) trim(desc)
       call mpi_finalize(i)
       stop
    end if

    ! Case 2: Include file error
2   write(*,*) "get_parameter: Fatal error: Cannot open include file '" // trim(value) // "'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    call mpi_finalize(i)
    stop

    ! Case 3: Directive error
3   write(*,*) "get_parameter: Fatal error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    call mpi_finalize(i)
    stop

    ! Case 4: Top level parameter file unreadable
4   write(*,*) "get_parameter: Fatal error: Cannot open parameter file '" // trim(parfile) // "'"
    call mpi_finalize(i)
    stop
  end subroutine

  subroutine get_parameter_arg(parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    character(len=512) :: line
    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, iargc()
       call getarg(i, line)
       if(line(1:2) /= "--") cycle
       call parse_parameter(line(3:), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arg: Fatal error: Cannot find parameter " // trim(parname) // " in argument list!"
       if(present(desc)) write(*,*) trim(desc)
       call mpi_finalize(i)
       stop
    end if
  end subroutine

  subroutine get_parameter_arr(arr, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present)
    implicit none
    character(len=*)           :: arr(:)
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present

    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, size(arr)
       call parse_parameter(arr(i), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arr: Fatal error: Cannot find parameter " // trim(parname) // " in argument list!"
       stop
    end if
  end subroutine


  subroutine parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
    implicit none
    character(len=*)           :: line, parname
    character(len=256)         :: toks(2), key, value, par
    logical(lgt)               :: found
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt

    integer(i4b) :: i, n

    call get_tokens(trim(line), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
    if(n < 2) then
       found = .false.
       return
    end if
    key = get_token(toks(1), " ", 1, group="''" // '""')
    value = get_token(toks(2), " ", 1, group="''" // '""')
    par = parname
    call tolower(key)
    call tolower(par)
    if (trim(key) == trim(par)) then
       if (present(par_int)) then
          read(value,*) par_int
       elseif (present(par_char)) then
          read(value,*) par_char
       elseif (present(par_string)) then
          read(value,*) par_string
       elseif (present(par_sp)) then
          read(value,*) par_sp
       elseif (present(par_dp)) then
          read(value,*) par_dp
       elseif (present(par_lgt)) then
          read(value,*) par_lgt
       else
          write(*,*) "get_parameter: Reached unreachable point!"
       end if
       found = .true.
    else
       found = .false.
    end if
  end subroutine


  ! Loops through the parameter files and children, counting lines.
  ! No error reporting.
  subroutine dump_expanded_paramfile(parfile, outfile)
    implicit none
    character(len=*)           :: parfile, outfile
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i, num, ounit
    character(len=1024)        :: key, value, arg
    character(len=1)           :: equals

    num = 0
    depth = 1
    ounit = getlun()
    open(ounit,file=outfile,action="write")
    write(ounit,fmt="(a)",advance="no") '# Arguments:'
    do i = 1, iargc()
       call getarg(i, arg)
       write(ounit,fmt="(a)",advance="no") " '" // trim(arg) // "'"
    end do
    write(ounit,*)

    units(depth) = getlun()
    open(units(depth),file=parfile,status="old",action="read")
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       backspace(units(depth))
       if (key=='@INCLUDE') then
          ! Recurse into the new file
          read(units(depth),*,end=1) key, value
          write(ounit,fmt='(a)') pad("",depth-1," ") // "# File: " // trim(value)
          depth=depth+1
          units(depth) = getlun()
          open(units(depth),file=value,status="old")
       else
          read(units(depth),fmt="(a)") value
          write(ounit,fmt='(a)') pad("",depth-1," ") // trim(value)
       end if
       cycle
1      close(units(depth))
       depth = depth-1
    end do
    close(ounit)
  end subroutine


  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num; call tokenize(string, sep, ext, group, allow_empty); end do
    res = string(ext(1):ext(2))
  end function

  ! Fill all tokens into toks, and the num filled into num
  subroutine get_tokens(string, sep, toks, num, group, maxnum, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*) :: toks(:)
    character(len=*), optional :: group
    integer(i4b),     optional :: num, maxnum
    logical(lgt),     optional :: allow_empty
    integer(i4b) :: n, ext(2), nmax
    ext = -1
    n = 0
    nmax = size(toks); if(present(maxnum)) nmax = maxnum
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0 .and. n < nmax)
       n = n+1
       toks(n) = string(ext(1):ext(2))
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    if(present(num)) num = n
  end subroutine

  function has_token(token, string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: token, string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    logical(lgt)     :: res
    integer(i4b)     :: ext(2)
    res = .true.
    ext = -1
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       if(string(ext(1):ext(2)) == trim(token)) return
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = .false.
  end function

  function num_tokens(string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)     :: res, ext(2)
    ext = -1
    res = 0
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       res = res+1
       call tokenize(string, sep, ext, group, allow_empty)
    end do
  end function

  subroutine tokenize(string, sep, ext, group, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional   :: group
    character(len=256)  :: op, cl
    integer(i4b), save           :: level(256), nl
    integer(i4b), intent(inout)  :: ext(2)
    logical(lgt), optional       :: allow_empty

    integer(i4b) :: i, j, k, n, o, c, ng
    logical(lgt) :: intok, hit, empty

    empty = .false.; if(present(allow_empty)) empty = allow_empty

    if(ext(2) >= len(string)) then
       ext = (/ 0, -1 /)
       return
    end if
    ng = 0
    if(present(group)) then
       ng = len_trim(group)/2
       do i = 1, ng
          op(i:i) = group(2*i-1:2*i-1)
          cl(i:i) = group(2*i:2*i)
       end do
    end if
    if(ext(2) <= 0) then
       level = 0
       nl = 0
    end if
    intok = .false.
    do i = ext(2)+2, len(string)
       hit = .false.
       c = index(cl(1:ng), string(i:i))
       if(c /= 0) then; if(level(c) > 0) then
          level(c) = level(c) - 1
          if(level(c) == 0) nl = nl - 1
          hit = .true.
       end if; end if
       if(nl == 0) then
          ! Are we in a separator or not
          if(index(sep, string(i:i)) == 0) then
             ! Nope, so we must be in a token. Register start of token.
             if(.not. intok) then
                j = i
                intok = .true.
             end if
          else
             ! Yes. This either means that a token is done, and we should
             ! return it, or that we are waiting for a new token, in
             ! which case do nothing.
             if(intok) then
                ext = (/ j, i-1 /)
                return
             elseif(empty) then
                ext = (/ i, i-1 /)
                return
             end if
          end if
       end if
       o = index(op(1:ng), string(i:i))
       if(o /= 0 .and. .not. hit) then
          if(level(o) == 0) nl = nl + 1
          level(o) = level(o) + 1
       end if
    end do
    ! Handle last token
    if(intok) then
       ext = (/ j, i-1 /)
    elseif(empty) then
       ext = (/ i, i-1 /)
    else
       ext = (/ 0, -1 /)
    end if
  end subroutine



  function pad(a, n, c) result(b)
    implicit none
    character(len=*) :: a
    character(len=1), optional :: c
    character(len=1):: k
    integer(i4b)    :: n, i, m
    character(len=max(len_trim(a),n)) :: b
    if(len_trim(a) > n) then
       b = a
    else
       k = '0'; if(present(c)) k = c
       m = len_trim(a)
       do i = 1, n-m; b(i:i) = k; end do
       b(n-m+1:n) = a(1:m)
    end if
  end function



  subroutine read_map(filename, map)
    implicit none

    character(len=128),                 intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(out) :: map

    integer(i4b)   :: nside, nmaps, ordering, int_npix, temp_i, i, npix, nmaps_in, nside_in
    real(dp)       :: nullval
    logical(lgt)   :: anynull

    npix  = size(map(:,1))
    nmaps = size(map(0,:))
    nside = nint(sqrt(real(npix,sp)/12.))

    temp_i = getsize_fits(trim(filename), ordering=ordering, nside=nside_in, nmaps=nmaps_in)
    if ((nmaps_in < nmaps) .or. (nside_in /= nside)) then
       write(*,*) 'Incorrect nside or nmaps for file called ', trim(filename)
    end if

    !call read_bintab(trim(filename), map, npix, nmaps, nullval, anynull)
    call input_map(filename, map, npix, nmaps)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if

  end subroutine read_map



  subroutine read_ringweights(nside, weights)
    implicit none

    integer(i4b),             intent(in)  :: nside
    real(dp), dimension(:,:), intent(out) :: weights

    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    logical(lgt)        :: exist, anynull
    integer(i4b)        :: nmaps
    real(dp)            :: nullval

    nmaps = size(weights(1,:))

    call int2string(nside, nside_text)
    weight_file = 'weight_ring_n' // nside_text // '.fits'
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(weight_file, weights, 2*nside, nmaps, nullval, anynull)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if

  end subroutine read_ringweights



  ! Routine for initializing the monopole and dipoles
  subroutine initialize_mono_and_dipole_temp(nside, templates, pixels)
    implicit none

    integer(i4b),                   intent(in)            :: nside
    real(dp),     dimension(0:,1:), intent(out)           :: templates
    integer(i4b), dimension(0:),    intent(in), optional  :: pixels


    integer(i4b) :: i, j, k, l, numpix
    real(dp), dimension(3) :: vector

    if (present(pixels)) then

       numpix = size(pixels)
       
       do i = 0, numpix-1
          call pix2vec_ring(nside, pixels(i), vector)
          
          templates(i,1) = 1.d0
          templates(i,2) = vector(1)
          templates(i,3) = vector(2)
          templates(i,4) = vector(3)
       end do

    else

       do i = 0, 12*nside**2-1
          call pix2vec_ring(nside, i, vector)
          
          templates(i,1) = 1.d0
          templates(i,2) = vector(1)
          templates(i,3) = vector(2)
          templates(i,4) = vector(3)
       end do

    end if

!    templates(:,2:4) = sqrt(3.d0) * templates(:,2:4)

  end subroutine initialize_mono_and_dipole_temp

!!$  ! Routine for estimating the best-fit monopole and dipole. Used for initializing
!!$  ! data maps, and the monopole and dipole sampling
!!$  subroutine remove_best_fit_mono_dipole(map)
!!$    implicit none
!!$
!!$    real(dp), dimension(0:), intent(inout) :: map
!!$
!!$    integer(i4b) :: i, j, k, nside, npix
!!$    real(dp), dimension(3)     :: vector
!!$
!!$    real(dp), dimension(4)   :: x, b, Y
!!$    real(dp), dimension(4,4) :: A
!!$
!!$    npix    = size(map)
!!$    nside   = nint(sqrt(real(npix,sp)/12.))
!!$
!!$    b = 0.d0
!!$    A = 0.d0
!!$    do i = 0, npix-1
!!$       call pix2vec_ring(nside, i, vector)
!!$
!!$       Y(1) =  1.d0
!!$       Y(2) =  vector(1)
!!$       Y(3) =  vector(2)
!!$       Y(4) =  vector(3)
!!$
!!$       do j = 1, 4
!!$          b(j) = b(j) + map(i) * Y(j)
!!$       end do
!!$       
!!$       do j = 1, 4
!!$          do k = 1, j
!!$             A(j,k) = A(j,k) + Y(j) * Y(k)
!!$          end do
!!$       end do
!!$       
!!$    end do
!!$
!!$    do j = 1, 4
!!$       do k = 1, j
!!$          A(k,j) = A(j,k)
!!$       end do
!!$    end do
!!$
!!$    ! Solve the linear system
!!$    call solve_system_real(A, x, b)
!!$
!!$    ! Remove monopole and dipole
!!$    do i = 0, npix-1
!!$       call pix2vec_ring(nside, i, vector)
!!$       map(i) = map(i) - x(1)
!!$       map(i) = map(i) - x(2) * vector(1)
!!$       map(i) = map(i) - x(3) * vector(2)
!!$       map(i) = map(i) - x(4) * vector(3)
!!$    end do
!!$
!!$  end subroutine remove_best_fit_mono_dipole
!!$
!!$  subroutine remove_best_fit_mono_dipole2(map, mask)
!!$    implicit none
!!$
!!$    real(dp), dimension(0:), intent(inout)  :: map
!!$    real(dp), dimension(0:), intent(in)     :: mask
!!$
!!$    integer(i4b) :: i, j, k, nside, npix
!!$    real(dp), dimension(3)     :: vector
!!$
!!$    real(dp), dimension(4)   :: x, b, Y
!!$    real(dp), dimension(4,4) :: A
!!$
!!$    npix    = size(map)
!!$    nside   = nint(sqrt(real(npix,sp)/12.))
!!$
!!$    b = 0.d0
!!$    A = 0.d0
!!$    do i = 0, npix-1
!!$       if (mask(i) > 0.) then
!!$
!!$          call pix2vec_ring(nside, i, vector)
!!$
!!$          Y(1) =  1.d0
!!$          Y(2) =  vector(1)
!!$          Y(3) =  vector(2)
!!$          Y(4) =  vector(3)
!!$
!!$          do j = 1, 4
!!$             b(j) = b(j) + mask(i) * map(i) * Y(j)
!!$          end do
!!$
!!$          do j = 1, 4
!!$             do k = 1, j
!!$                A(j,k) = A(j,k) + mask(i) * Y(j) * Y(k)
!!$             end do
!!$          end do
!!$
!!$       end if
!!$    end do
!!$
!!$    do j = 1, 4
!!$       do k = 1, j
!!$          A(k,j) = A(j,k)
!!$       end do
!!$    end do
!!$
!!$    ! Solve the linear system
!!$    call solve_system_real(A, x, b)
!!$
!!$    write(*,*) 'md = ', real(x,sp)
!!$
!!$    ! Remove monopole and dipole
!!$    do i = 0, npix-1
!!$       if (mask(i) > 0.) then
!!$          call pix2vec_ring(nside, i, vector)
!!$          map(i) = map(i) - x(1)
!!$          map(i) = map(i) - x(2) * vector(1)
!!$          map(i) = map(i) - x(3) * vector(2)
!!$          map(i) = map(i) - x(4) * vector(3)
!!$       else
!!$          map(i) = 0.d0
!!$       end if
!!$    end do
!!$
!!$  end subroutine remove_best_fit_mono_dipole2


  subroutine recover_session_data(chain, first_iteration, s, cl, foreground_coeff)
    implicit none

    integer(i4b),                      intent(in)            :: chain
    integer(i4b),                      intent(out)           :: first_iteration
    complex(dpc), dimension(1:,1:),    intent(out)           :: s
    real(dp),     dimension(0:,1:,1:), intent(out)           :: cl
    real(dp),     dimension(1:,1:),    intent(out), optional :: foreground_coeff

    integer(i4b)       :: i, j, l, m, numalms, lmax, numband, numtemp, iteration, nlheader, unit
    logical(lgt)       :: exist
    character(len=4)   :: chain_text
    character(len=5)   :: k_text
    character(len=128) :: filename
    real(dp),     allocatable, dimension(:,:,:) :: alms_in
    complex(dpc), allocatable, dimension(:,:,:) :: alms_cmplx
    character(len=128), dimension(1:80,1)       :: header

    unit            = getlun()
    first_iteration = 1
    s               = cmplx(0.d0, 0.d0)
    cl              = 0.d0
    if (present(foreground_coeff)) foreground_coeff = 0.d0

!!$    call int2string(chain, chain_text)
!!$
!!$    iteration = 1
!!$    exist     = .true.
!!$    do while (exist)
!!$       call int2string(iteration, k_text)
!!$       filename = 'alms_c' // chain_text // '_k' // k_text // '.fits'
!!$       inquire(file=trim(filename), exist=exist)
!!$       if (.not. exist .and. iteration == 1) then
!!$          write(*,*) 'Chain no. ', chain, ' -- No recovery data found.'
!!$          stop
!!$       end if
!!$       iteration = iteration + 1
!!$    end do
!!$    iteration = iteration - 2
!!$
!!$    write(*,*) 'Chain no. ', chain, ' -- Recovering from sample no. ', iteration
!!$
!!$    first_iteration = iteration+1
!!$    
!!$    call int2string(iteration, k_text)
!!$    
!!$    ! Read power spectrum
!!$    filename = 'cl_c' // chain_text // '_k' // k_text // '.dat'       
!!$    open(unit, file=trim(filename))
!!$    lmax = size(cl)-1
!!$    do l = 0, lmax
!!$       read(unit,*) i, cl(l)
!!$       if (l > 1) then
!!$          cl(l) = cl(l) * (2.d0*pi) / real(l*(l+1),dp)
!!$       else
!!$          cl(l) = 0.d0
!!$       end if
!!$    end do
!!$    close(unit)
!!$    
!!$    ! Read alms
!!$    filename = 'alms_c' // chain_text // '_k' // k_text // '.fits'       
!!$    
!!$    numalms = (lmax+1)*(lmax+2)/2
!!$    allocate(alms_in(numalms,4,1))
!!$    allocate(alms_cmplx(1,0:lmax,0:lmax))
!!$    nlheader = size(header)
!!$    call fits2alms(filename, numalms, alms_in, 3, header, nlheader, 1)
!!$    
!!$    alms_cmplx = cmplx(0.d0,0.d0)
!!$    do i = 1, numalms
!!$       l = nint(alms_in(i,1,1))
!!$       m = nint(alms_in(i,2,1))
!!$       alms_cmplx(1,l,m) = cmplx(alms_in(i,3,1), alms_in(i,4,1))
!!$    end do
!!$    
!!$    call convert_complex_to_real_alms(alms_cmplx, s)       
!!$    
!!$    deallocate(alms_in)
!!$    deallocate(alms_cmplx)
!!$    
!!$    
!!$    if (present(foreground_coeff)) then
!!$       ! Read foreground coefficients
!!$       filename = 'foreground_c' // chain_text // '_k' // k_text // '.dat'       
!!$       open(unit, file=trim(filename))
!!$       numtemp = size(foreground_coeff(:,1))
!!$       numband = size(foreground_coeff(1,:))
!!$       do i = 1, numband
!!$          do j = 1, numtemp
!!$             read(unit,*) foreground_coeff(j,i)
!!$          end do
!!$       end do
!!$       close(unit)
!!$    end if

  end subroutine recover_session_data

!  subroutine compute_lowres_fg_amp(fg_amp, lmax_lowres, fwhm_lowres, fg_amp_lowres)
  subroutine compute_lowres_fg_amp(fg_amp, fg_amp_lowres)
    implicit none

!    integer(i4b),                      intent(in)  :: lmax_lowres
!    real(dp),                          intent(in)  :: fwhm_lowres
    real(dp),     dimension(0:,1:,1:), intent(in)  :: fg_amp
    real(dp),     dimension(0:,1:,1:), intent(out) :: fg_amp_lowres

    integer(i4b) :: i, npix_in, npix_out, nmaps
!    real(dp), allocatable, dimension(:,:) :: beam

    npix_in  = size(fg_amp(:,1,1))
    npix_out = size(fg_amp_lowres(:,1,1))
    nmaps    = size(fg_amp(0,:,1))
    
    if (npix_in == npix_out) then
       fg_amp_lowres = fg_amp
    else

       write(*,*) 'Different high- and low-resolution analysis not supported yet.'
       stop

!       allocate(beam(0:lmax_lowres,nmaps))
!       beam = 1.d0

!       do i = 1, size(fg_amp_lowres(0,1,:))
!          call degrade_map(fg_amp(:,:,i), beam, fwhm_lowres, fg_amp_lowres(:,:,i))
!       end do

!       deallocate(beam)
    end if

  end subroutine compute_lowres_fg_amp


  subroutine write_map(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map
    character(len=*),                   intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),                           intent(in), optional :: nu_ref

    integer(i4b)   :: npix, nlheader, nmaps, i, nside
    logical(lgt)   :: exist, polarization

    character(len=80), dimension(1:120)    :: header
    character(len=16) :: unit_, ttype_

    npix         = size(map(:,1))
    nside        = nint(sqrt(real(npix,sp)/12.))
    nmaps        = size(map(0,:))
    polarization = (nmaps == 3)
    unit_        = '';       if (present(unit)) unit_  = unit
    ttype_       = 'Stokes'; if (present(unit)) ttype_ = ttype


    !-----------------------------------------------------------------------
    !                      write the map to FITS file
    !  This is copied from the synfast.f90 file in the Healpix package
    !-----------------------------------------------------------------------
    
    nlheader = SIZE(header)
    do i=1,nlheader
       header(i) = ""
    enddo

    ! start putting information relative to this code and run
    call add_card(header)
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
    call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
    call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
    call add_card(header) ! blank line
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"POLCCONV","COSMO"," Coord. convention for polarisation (COSMO/IAU)")
    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
    call add_card(header) ! blank line
    call add_card(header,"POLAR",polarization," Polarisation included (True/False)")

    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "I_"//ttype_,"Stokes I")
    call add_card(header,"TUNIT1", unit_,"Map unit")
    call add_card(header)

    if (polarization) then
       call add_card(header,"TTYPE2", "Q_"//ttype_,"Stokes Q")
       call add_card(header,"TUNIT2", unit_,"Map unit")
       call add_card(header)
       
       call add_card(header,"TTYPE3", "U_"//ttype_,"Stokes U")
       call add_card(header,"TUNIT3", unit_,"Map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Commander Keywords                        ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","Commander is a code for global CMB analysis    ")
    call add_card(header,"COMMENT","developed in collaboration between the University")
    call add_card(header,"COMMENT","of Oslo and Jet Propulsion Laboratory (NASA).  ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    if (present(comptype)) call add_card(header,"COMPTYPE",trim(comptype), "Component type")
    if (present(nu_ref))   call add_card(header,"NU_REF",  nu_ref,         "Reference frequency")
    if (present(spectrumfile)) call add_card(header,"SPECFILE",  trim(spectrumfile), &
         & "Reference spectrum")
    call add_card(header,"COMMENT","-----------------------------------------------")


    !call write_bintab(map, npix, nmaps, header, nlheader, filename)
    call output_map(map, header, "!"//trim(filename))

  end subroutine write_map


  subroutine write_s_lm(filename, alms)
    implicit none
    
    character(len=128),                      intent(in) :: filename
    complex(dpc),       dimension(1:,0:,0:), intent(in) :: alms

    integer(i4b) :: i, l, m, lmax, nmaps, nlheader
    logical(lgt) :: exist, polarization
    character(len=80), dimension(1:120) :: header

    lmax  = size(alms(1,:,0))-1
    nmaps = size(alms(:,0,0))

    inquire(file=trim(filename),exist=exist)
    if (exist) call system('rm ' // trim(filename))
    
    ! write each component (T,G,C) in a different extension of the same file
    do i = 1, nmaps
       header = ""
       call add_card(header,"EXTNAME","SIMULATION")
       call add_card(header,"COMMENT","-----------------------------------------------")
       call add_card(header,"COMMENT","     Map Simulation Specific Keywords      ")
       call add_card(header,"COMMENT","-----------------------------------------------")
       if (i == 1) then
          call add_card(header,"EXTNAME","'Commander Gibbs sample a_lms (TEMPERATURE)'")
       elseif (i == 2) then
          call add_card(header,"EXTNAME","'Commander Gibbs sample a_lms (GRAD component)'")
       elseif (i == 3) then
          call add_card(header,"EXTNAME","'Commander Gibbs sample a_lms (CURL component)'")
       endif
       call add_card(header)
       call add_card(header,"MAX-LPOL",lmax      ,"Maximum multipole l")
       call add_card(header,"MAX-MPOL",lmax      ,"Maximum m")
       call add_card(header,"HISTORY")
       call add_card(header,"COMMENT","-----------------------------------------------")
       if (i == 1) then
          call add_card(header,"COMMENT","Temperature a_lms")
       else
          call add_card(header,"COMMENT","Polarisation a_lms")
       endif
       call add_card(header,"COMMENT"," The real and imaginary part of the a_lm with m>=0")
       call add_card(header)
       call add_card(header,"TTYPE1", "INDEX"," i = l^2 + l + m + 1")
       call add_card(header,"TUNIT1", "   "," index")
       call add_card(header)
       call add_card(header,"TTYPE2", "REAL"," REAL a_lm")
       call add_card(header,"TUNIT2", 'muK'," alm units")
       call add_card(header)
       !
       call add_card(header,"TTYPE3", "IMAG"," IMAGINARY a_lm")
       call add_card(header,"TUNIT3", 'muK'," alm units")
       call add_card(header)
       !
       call add_card(header)
       call add_card(header)
       call add_card(header)
       !
       nlheader = SIZE(header)
       call dump_alms (filename,alms(i,0:lmax,0:lmax),lmax,header,nlheader,i-1)
    enddo
    
  end subroutine write_s_lm




  ! Routine for writing the results into a FITS-file in HKE's standard result file format
!!$  subroutine write_resfile(filename, lmax, nspec, numsamp, numchain, funcname, data)
!!$    implicit none
!!$
!!$    character(len=*),             intent(in) :: filename, funcname
!!$    integer(i4b),                 intent(in) :: lmax, nspec, numsamp, numchain
!!$    real(sp), dimension(:,:,:,:), intent(in) :: data
!!$
!!$    integer(i4b)         :: status, blocksize, readwrite, unit
!!$    integer(i4b)         :: nelements, fpixel, group, bitpix, naxis
!!$    logical(lgt)         :: simple, extend, anyf, exist
!!$    real(sp)             :: nullval
!!$    character(len=80)    :: comment, errorline
!!$    integer(i4b), dimension(4) :: naxes
!!$
!!$    unit = getlun()
!!$
!!$    ! Output to a fits-file
!!$    status = 0
!!$    bitpix = -32
!!$    simple = .true.
!!$    naxis = 4
!!$    extend = .true.
!!$
!!$    naxes(1) = size(data(:,1,1,1))
!!$    naxes(2) = size(data(1,:,1,1))
!!$    naxes(3) = size(data(1,1,:,1))
!!$    naxes(4) = size(data(1,1,1,:))
!!$
!!$    inquire(file=trim(filename),exist=exist)
!!$    if (exist) call system('mv ' // trim(filename) // ' ' // trim(filename) // '_old')
!!$
!!$    ! Open a new fitsfile
!!$    blocksize = 1
!!$    call ftinit(unit,trim(filename),blocksize,status)
!!$
!!$    ! Write the required headers
!!$    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!!$
!!$    ! Write a few keywords
!!$    call ftpkys(unit,'FUNCNAME',trim(funcname),'Full function name',status)
!!$    call ftpkyj(unit,'LMAX',lmax,'Maximum multipole moment',status)
!!$    call ftpkyj(unit,'NUMSAMP',numsamp,'Number of samples',status)
!!$    call ftpkyj(unit,'NUMCHAIN',numchain,'Number of independent chains',status)
!!$    call ftpkyj(unit,'NUMSPEC',nspec,'Number of power spectra',status)
!!$
!!$    ! Output the array
!!$    group=1
!!$    fpixel=1
!!$    nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
!!$    call ftppre(unit,group,fpixel,nelements,data,status)
!!$
!!$    ! Write the current date
!!$    call ftpdat(unit,status)
!!$
!!$    ! close the file and free the unit number
!!$    call ftclos(unit, status)
!!$    call ftfiou(unit, status)
!!$
!!$    if (status /= 0) then
!!$       write(*,*) 'Error in FITS operations -- status = ', status
!!$
!!$       call ftgmsg(errorline)
!!$       do while (trim(errorline) /= '')
!!$          write(*,*) errorline
!!$          call ftgmsg(errorline)
!!$       end do
!!$
!!$    end if
!!$
!!$  end subroutine write_resfile


!!$  subroutine check_mpi_ierr(ierr)
!!$    implicit none
!!$
!!$    integer(i4b), intent(in) :: ierr
!!$
!!$    if (ierr /= 0) then
!!$       write(*,*) 'Error in MPI routine. ierr = ', ierr
!!$       stop
!!$    end if
!!$
!!$  end subroutine check_mpi_ierr


  subroutine read_pixwin(nside, nmaps, pixwin)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: pixwin

    integer(i4b)        :: nc
    character(len=128)  :: pixwin_file
    character(len=4)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(pixwin(0:4*nside,nmaps))
    
    if (nmaps == 3) then
       nc = 2
    else
       nc = 1
    end if

    call int2string(nside, nside_text)
    pixwin_file = 'pixel_window_n' // nside_text // '.fits'
    inquire(file=pixwin_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside+1, nc, nullval, anynull)
       if (nmaps == 3) pixwin(:,3) = pixwin(:,2)
    else
       pixwin = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Pixel window file ', trim(pixwin_file), ' not found. '
       write(*,*) 'Using unity weights.'
       write(*,*) ''
    end if

  end subroutine read_pixwin


  function matrix_is_positive_definite(nonzero_component, cl)
    implicit none

    logical(lgt), dimension(1:),    intent(in) :: nonzero_component
    real(dp),     dimension(1:,1:), intent(in) :: cl
    logical(lgt)                               :: matrix_is_positive_definite

    integer(i4b) :: i, j, nmaps

    nmaps = size(nonzero_component)

    if (nmaps == 1) then

       matrix_is_positive_definite = (cl(1,1) > 0.d0)
       return

    else if (nmaps == 3) then

       matrix_is_positive_definite = .true.
       do i = 1, nmaps
          if (nonzero_component(i)) then
             if (cl(i,i) < 0.d0) then
                matrix_is_positive_definite = .false.
                return
             end if
             do j = i+1, nmaps
                if (cl(i,i)*cl(j,j) < cl(i,j)**2) then
                   matrix_is_positive_definite = .false.
                   return
                end if
             end do
          end if
       end do

    else
       write(*,*) 'Invalid nmaps = ', nmaps
       stop
    end if
    

  end function matrix_is_positive_definite


  function is_nan(a) result(res)
    implicit none
    real(dp) :: a
    logical(lgt) :: res
    res = (a .ne. a) .or. (a .eq. NaN)
  end function

  subroutine cl2s(cl, S)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: cl
    real(dp), dimension(1:,1:), intent(out) :: S

    integer(i4b) :: i, j, k, nmaps

    nmaps = size(S,1)

    S = 0.d0
    k = 1
    do i = 1, nmaps
       do j = i, nmaps
          S(i,j) = cl(k)
          S(j,i) = cl(k)
          k      = k+1
       end do
    end do

  end subroutine cl2s

  subroutine s2cl(S, cl)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: S
    real(dp), dimension(1:),    intent(out) :: cl

    integer(i4b) :: i, j, k, nmaps

    nmaps = size(S,1)

    k = 1
    do i = 1, nmaps
       do j = i, nmaps
          cl(k) = S(i,j)
          k      = k+1
       end do
    end do

  end subroutine s2cl

  function is_posdef(cl)
    implicit none

    real(dp),     dimension(1:), intent(in) :: cl
    logical(lgt)                            :: is_posdef

    integer(i4b) :: nmaps, ind(3), i, n
    real(dp)     :: W(3), S(3,3)

    nmaps = min(size(cl),3)
    if (nmaps == 1) then
       is_posdef = (cl(1) > 0.d0)
    else
       call cl2s(cl, S)
       n = 0
       do i = 1, nmaps
          if (S(i,i) /= 0.d0) then
             n      = n+1
             ind(n) = i
          end if
       end do
       call get_eigenvalues(S(ind(1:n),ind(1:n)),W(1:n))
       is_posdef = all(W(1:n) > 0.d0)
    end if

  end function is_posdef

  subroutine mkdir(path)
    implicit none
    character(len=*) :: path
    logical(lgt) :: exist
    integer(i4b) :: i, left
    exist = is_dir(path)
    if(exist) return
    if(mkdir_c(trim(path))) return
    do i = 1, 20
       call fsleep(1d0)
       if(is_dir(path)) return
       if(mkdir_c(trim(path))) return
    end do
    write(*,*) "Failed to create directory '" // trim(path) 
  end subroutine

  subroutine mkdirs(path, skiplast)
    implicit none
    character(len=*) :: path
    logical(lgt) :: skiplast
    integer(i4b) :: i,j,k,n,m
    n = len_trim(path)
    if(skiplast) then
       do while(path(n:n) /= '/')
          n = n-1
          if(n <= 1) return
       end do
       n = n-1
    end if

    j = 1
    do while(.true.)
       do while(path(j:j) == '/'); j = j+1; if(j > n) return; end do
       i = j
       do while(path(j:j) /= '/'); j = j+1; if(j > n) exit; end do
       call mkdir(path(1:j-1))
       if(j > n) return
    end do
  end subroutine

  subroutine mv(from, to)
    character(len=*) :: from, to
    logical(lgt) :: error
    error = .not. mv_c(trim(from), trim(to))
    if(error) write(*,*) "Failed to move '" // trim(from) // "' to '" // trim(to) // "'!"
  end subroutine

  subroutine rm(file, noerr)
    character(len=*) :: file
    logical(lgt) :: error, check
    logical(lgt), optional :: noerr
    check = .true.; if(present(noerr)) check = .not. noerr
    error = .not. rm_c(trim(file))
    if(check .and. error) write(*,*) "Failed to remove '" // trim(file) // "'!"
  end subroutine

  subroutine touch(file)
    character(len=*) :: file
    integer(i4b) :: unit
    unit = getlun()
    open(unit,file=file)
    close(unit)
  end subroutine

  ! This will be slightly slow if the compiler doesn't optimize the
  ! allocation and deallocation here. Could destroy input array if this
  ! is a problem.
  function median(array) result(res)
    implicit none
    real(dp) :: array(:), res
    real(dp), dimension(:), allocatable :: tmp

    allocate(tmp(size(array)))
    tmp = array
    call QuickSort_real(tmp)
    res = tmp(size(tmp)/2+1)
    deallocate(tmp)
  end function median

end module comm_utils
