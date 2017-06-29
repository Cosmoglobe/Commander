module comm_utils
  use alm_tools
  use rngmod
  use fitstools
  use head_fits
  use pix_tools
  use udgrade_nr
  use iso_c_binding
  use comm_system_mod
  use comm_defs
  use sort_utils
  use spline_1D_mod
  implicit none

  include "mpif.h"
  
  ! *********************************************************************
  ! *      Commander -- An MCMC code for global, exact CMB analysis     *
  ! *                                                                   *
  ! *                 Written by Hans Kristian Eriksen                  *
  ! *                                                                   *
  ! *                Copyright 2015, all rights reserved                *
  ! *                                                                   *
  ! *********************************************************************

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

  subroutine read_beam(lmax, nmaps, beam, beamfile, fwhm, pixwin)
    implicit none

    integer(i4b),                                    intent(in)            :: lmax, nmaps
    real(dp),           allocatable, dimension(:,:), intent(out)           :: beam
    real(dp),                                        intent(in), optional  :: fwhm
    character(len=*),                                intent(in), optional  :: beamfile, pixwin

    real(dp),          allocatable, dimension(:,:)   :: beam_in, pw
    character(len=80),              dimension(1:180) :: header

    allocate(beam(0:lmax,nmaps))
    beam = 1.d0

    allocate(pw(0:lmax,2), beam_in(0:lmax,4))
    if (present(pixwin)) then
       call fits2cl(pixwin, pw, lmax, 2, header)
    else
       pw = 1.d0
    end if
    
    if (present(beamfile)) then
       call fits2cl(beamfile, beam_in, lmax, 4, header)
    else if (present(fwhm)) then
       call gaussbeam(fwhm, lmax, beam_in(0:lmax,1:nmaps))
    else
       beam_in = 1.d0
    end if

    beam(:,1) = beam_in(:,1) * pw(:,1)
    if (nmaps > 1) then
       if (sum(beam_in(:,2)) < 1.d0) then
          ! Default to temperature beam+polarization pixel window
          beam(:,2) = beam(:,1)*pw(:,2)
          beam(:,3) = beam(:,1)*pw(:,2)
       else
          ! Use polarized beam+polarization pixel window
          beam(:,2) = beam_in(:,2)*pw(:,2)
          beam(:,3) = beam_in(:,3)*pw(:,2)
       end if
    end if

    deallocate(beam_in, pw)

  end subroutine read_beam

!!$  subroutine convert_RMS_to_invN(mask, rms, invN, sqrt)
!!$    implicit none
!!$
!!$    logical(lgt), dimension(0:,1:), intent(in)  :: mask
!!$    real(dp),     dimension(0:,1:), intent(in)  :: rmsmap
!!$    real(dp),     dimension(0:,1:), intent(out) :: inv_N_covar
!!$
!!$    integer(i4b) :: i, j, map_size, nmaps
!!$
!!$    map_size = size(mask(:,1))
!!$    nmaps    = size(rmsmap(0,:))
!!$
!!$    inv_N_covar = 0.d0
!!$    do j = 1, nmaps
!!$       do i = 0, map_size-1
!!$          if (mask(i,j)) then
!!$             inv_N_covar(i,j) = 1.d0 / rmsmap(i,j)**2
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end subroutine convert_RMS_to_invN


  subroutine compute_sigma_l_mat(s_i, cl_i)
    implicit none

    real(dp), dimension(1:,1:),    intent(in)   :: s_i
    real(dp), dimension(0:,1:,1:), intent(out)  :: cl_i

    integer(i4b) :: i, j, k, l, m, nmaps, lmax, nspec
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

!!$  function is_posdef(cl)
!!$    implicit none
!!$
!!$    real(dp),     dimension(1:), intent(in) :: cl
!!$    logical(lgt)                            :: is_posdef
!!$
!!$    integer(i4b) :: nmaps, ind(3), i, n
!!$    real(dp)     :: W(3), S(3,3)
!!$
!!$    nmaps = min(size(cl),3)
!!$    if (nmaps == 1) then
!!$       is_posdef = (cl(1) > 0.d0)
!!$    else
!!$       call cl2s(cl, S)
!!$       n = 0
!!$       do i = 1, nmaps
!!$          if (S(i,i) /= 0.d0) then
!!$             n      = n+1
!!$             ind(n) = i
!!$          end if
!!$       end do
!!$       call get_eigenvalues(S(ind(1:n),ind(1:n)),W(1:n))
!!$       is_posdef = all(W(1:n) > 0.d0)
!!$    end if
!!$
!!$  end function is_posdef

  subroutine mkdir(path)
    implicit none
    character(len=*) :: path
    logical(lgt) :: exist
    integer(i4b) :: i
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
    integer(i4b) :: i,j,n
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
!!$  function median(array) result(res)
!!$    implicit none
!!$    real(dp) :: array(:), res
!!$    real(dp), dimension(:), allocatable :: tmp
!!$
!!$    allocate(tmp(size(array)))
!!$    tmp = array
!!$    call QuickSort_real(tmp)
!!$    res = tmp(size(tmp)/2+1)
!!$    deallocate(tmp)
!!$  end function median

  subroutine allocate_map(comm, nside, nmaps, k, map, np, rings, pixels, filename)
    implicit none

    integer(i4b),                                  intent(in)  :: comm, nside, nmaps, k
    real(dp),         allocatable, dimension(:,:), intent(out) :: map
    integer(i4b),                                  intent(out), optional :: np
    integer(i4b),     allocatable, dimension(:),   intent(out), optional :: rings, pixels
    character(len=*),                              intent(in),  optional :: filename

    integer(i4b) :: i, npix

    npix = 12*nside**2
    allocate(map(0:npix-1,nmaps))
    if (present(filename)) then
       call read_map(filename, map)
    else
       map = 0.d0
    end if

    if (present(rings)) then
       if (allocated(rings)) deallocate(rings)
       allocate(rings(2*nside))
       do i = 1, 2*nside
          rings(i) = i
       end do
    end if

    if (present(pixels)) then
       if (allocated(pixels)) deallocate(pixels)
       allocate(pixels(1:npix))
       do i = 0, npix-1
          pixels(i+1) = i
       end do
    end if

    if (present(np)) np = npix
    
    
  end subroutine allocate_map

  
  subroutine read_map(filename, map)
    implicit none

    character(len=128),                 intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(out) :: map

    integer(i4b)   :: nside, nmaps, ordering, i, npix, nmaps_in, nside_in
    integer(i8b)   :: temp_i

    npix  = size(map(:,1))
    nmaps = size(map(0,:))
    nside = nint(sqrt(real(npix,sp)/12.))

    temp_i = getsize_fits(trim(filename), ordering=ordering, nside=nside_in, nmaps=nmaps_in)
    if ((nmaps_in < nmaps) .or. (nside_in /= nside)) then
       write(*,*) 'Incorrect nside or nmaps for file called ', trim(filename)
    end if

    call input_map(filename, map, npix, nmaps)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if

  end subroutine read_map

    ! Routine for initializing the monopole and dipoles
  subroutine initialize_mono_dipole(nside, pix, monopole, dipole)
    implicit none

    integer(i4b),                   intent(in)            :: nside
    integer(i4b), dimension(1:),    intent(in)            :: pix
    real(dp),     dimension(1:),    intent(out), optional :: monopole
    real(dp),     dimension(1:,1:), intent(out), optional :: dipole

    integer(i4b) :: i, np
    real(dp), dimension(3) :: vec

    np = size(pix)
       
    do i = 1, np
       if (present(monopole)) monopole(i) = 1.d0
       if (present(dipole)) then
          call pix2vec_ring(nside, pix(i), vec)
          dipole(i,1:3) = vec
       end if
    end do

  end subroutine initialize_mono_dipole

  function tsum(x, y)
    implicit none

    real(dp), dimension(:), intent(in) :: x, y
    real(dp)                           :: tsum

    integer(i4b) :: i

    tsum = 0.d0
    do i = 1, size(x)-1
       tsum = tsum + 0.5d0 * (y(i)+y(i+1)) * (x(i+1)-x(i))
    end do

  end function tsum

  
  subroutine report_error(message)
    implicit none

    character(len=*), intent(in) :: message

    integer(i4b) :: ierr, myid

    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    
    if (myid == 0) write(*,*) trim(message)
    call mpi_finalize(ierr)
    stop
    
  end subroutine report_error

  function mpi_dot_product(comm, x, y)
    implicit none

    integer(i4b),               intent(in) :: comm
    real(dp),     dimension(:), intent(in) :: x, y
    real(dp)                               :: mpi_dot_product

    integer(i4b) :: ierr
    real(dp)     :: prod

    prod = dot_product(x,y)
    call mpi_allreduce(MPI_IN_PLACE, prod, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    mpi_dot_product = prod

  end function mpi_dot_product



  ! Routine for reading a spectrum file
  subroutine read_spectrum(filename, spectrum)
    implicit none

    character(len=*),                               intent(in)  :: filename
    real(dp),           allocatable, dimension(:,:)             :: spectrum

    real(dp)            :: S_ref
    integer(i4b)        :: i, j, numpoint, unit, first, last
    character(len=128)  :: nu, val, string

    unit = getlun()
    open(unit, file=trim(filename))

    ! Find the number of entries
    numpoint = 0
    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*,end=1) nu, val
       numpoint = numpoint + 1
    end do
1   close(unit)

    if (numpoint == 0) then
       write(*,*) 'No valid data entries in spectrum file ', trim(filename)
       stop
    end if

    allocate(spectrum(numpoint,2))
    i = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=2) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*,end=1) nu, val
       i = i+1

       read(nu,*)  spectrum(i,1)
       read(val,*) spectrum(i,2)

       do j = 1, i-1
          if (spectrum(i,1) == spectrum(j,1)) then
             write(*,*) 'ERROR: Spectrum file ', trim(filename), ' contains double entries; j = ', j, ', spectrum = ', spectrum(i,1)
             stop
          end if
       end do
  
    end do
2   close(unit)

    ! Convert from GHz to Hz
    spectrum(:,1) = spectrum(:,1) * 1.d9

  end subroutine read_spectrum

  subroutine read_instrument_file(filename, field, label, default, val)
    implicit none
    character(len=*), intent(in)  :: filename, field, label
    real(dp),         intent(in)  :: default
    real(dp),         intent(out) :: val
    
    integer(i4b)        :: unit, i, j
    character(len=512)  :: band
    character(len=1024) :: line
    real(dp)            :: vals(2)
    logical(lgt)        :: ok

    val = default
    inquire(file=trim(filename), exist=ok)
    if (ok) then
       unit = getlun()
       open(unit, file=trim(filename), recl=1024)
       do while (.true.)
          read(unit,'(a)',end=99) line
          line = trim(adjustl(line))
          if (line(1:1) == '#') cycle
          read(line,*) band, vals
          if (trim(band) == trim(label)) then
             if (trim(field) == 'gain') then
                val = vals(1)
             else if (trim(field) == 'delta') then
                val = vals(2)
             else
                call report_error('Unsupported instrument field = '//trim(field))
             end if
             goto 99
          end if
       end do
99     close(unit)
    end if

  end subroutine read_instrument_file

  subroutine assert(condition, error_message)
    implicit none
    logical(lgt) :: condition
    character(len=*) :: error_message
    if(condition) return
    write(*,fmt="(a)") error_message
    stop
  end subroutine


  subroutine compute_radial_beam(nmaps, bl, br)
    implicit none
    integer(i4b),                                     intent(in)  :: nmaps
    real(dp),                       dimension(0:,1:), intent(in)  :: bl
    type(spline_type), allocatable, dimension(:),     intent(out) :: br

    integer(i4b)  :: i, j, k, m, n, l, lmax
    real(dp)      :: theta_max, threshold
    real(dp), allocatable, dimension(:)   :: x, pl
    real(dp), allocatable, dimension(:,:) :: y
    
    n         = 1000
    lmax      = size(bl,1)-1
    threshold = 1.d-6

    ! Find typical size
    l    = 0
    do while (bl(l,1) > 0.5d0)
       l = l+1
       if (l == lmax) call report_error('Error: Beam has not fallen off to 0.5 by lmax')
    end do
    theta_max = pi/l * 10.d0
    
    ! Compute radial beams
    allocate(x(n), y(n,nmaps), pl(0:lmax))
    do i = 1, n
       x(i) = theta_max/(n-1) * real(i-1,dp)
       call comp_normalised_Plm(lmax, 0, x(i), pl)
       do j = 1, nmaps
          y(i,j) = 0.d0
          do l = 0, lmax
             y(i,j) = y(i,j) + bl(l,j)*pl(l)/sqrt(4.d0*pi/real(2*l+1,dp))
          end do
       end do
    end do

    ! Spline significant part of beam profile
    allocate(br(nmaps))
    do j = 1, nmaps
       y(:,j) = y(:,j) / maxval(y(:,j))
       m      = 0
       do while (y(m+1,j) > threshold .and. m < n)
          m = m+1
       end do
       call spline(br(j), x(1:m), y(1:m,j))
    end do

    deallocate(x, y, pl)
    
  end subroutine compute_radial_beam

  function rand_trunc_gauss(rng_handle, mu, mu_, sigma)
    implicit none

    real(dp),         intent(in)    :: mu, mu_, sigma
    type(planck_rng), intent(inout) :: rng_handle
    real(dp)                        :: rand_trunc_gauss

    integer(i4b) :: n
    real(dp) :: eta, alpha, z, rho, u, mu0
    logical(lgt) :: ok

    mu0   = (mu_ - mu)/sigma
    alpha = 0.5d0 * (mu0 + sqrt(mu0**2 + 4.d0))
    do while (.true.)
       z   = -log(rand_uni(rng_handle)) / alpha + mu0
       rho = exp(-0.5d0 * (z-alpha)**2)
       if (rand_uni(rng_handle) < rho) then
          rand_trunc_gauss = sigma*z + mu
          exit
       end if
    end do

  end function rand_trunc_gauss
  
  
  
end module comm_utils
