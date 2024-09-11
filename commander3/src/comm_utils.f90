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
  use comm_mpi_mod
  use locate_mod
  use math_tools
  use powell_mod
  use hmc_mod
  implicit none

  !include "mpif.h"
  include 'fftw3.f'

contains

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
       !write(*,*) 'normalize beam to dipole'
       !beam_in(0:lmax,1) = beam_in(0:lmax,1) / beam_in(1,1)
    else
       beam_in = 1.d0
    end if

    beam(:,1) = beam_in(:,1) * pw(:,1)
    if (nmaps > 1) then
       if (sum(beam_in(:,2)) < 1.d0) then
          ! Default to temperature beam+polarization pixel window
          beam(:,2) = beam_in(:,1)*pw(:,2)
          beam(:,3) = beam_in(:,1)*pw(:,2)
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
  ! WARNING TO FUTURE USERS: The length of the input string is very important
  ! Only provide the length of string you think you will need, no extra space
  ! at the end. It will crash if the string is too long

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

!!$  subroutine allocate_map(comm, nside, nmaps, k, map, np, rings, pixels, filename)
!!$    implicit none
!!$
!!$    integer(i4b),                                  intent(in)  :: comm, nside, nmaps, k
!!$    real(dp),         allocatable, dimension(:,:), intent(out) :: map
!!$    integer(i4b),                                  intent(out), optional :: np
!!$    integer(i4b),     allocatable, dimension(:),   intent(out), optional :: rings, pixels
!!$    character(len=*),                              intent(in),  optional :: filename
!!$
!!$    integer(i4b) :: i, npix
!!$
!!$    npix = 12*nside**2
!!$    allocate(map(0:npix-1,nmaps))
!!$    if (present(filename)) then
!!$       call read_map(filename, map)
!!$    else
!!$       map = 0.d0
!!$    end if
!!$
!!$    if (present(rings)) then
!!$       if (allocated(rings)) deallocate(rings)
!!$       allocate(rings(2*nside))
!!$       do i = 1, 2*nside
!!$          rings(i) = i
!!$       end do
!!$    end if
!!$
!!$    if (present(pixels)) then
!!$       if (allocated(pixels)) deallocate(pixels)
!!$       allocate(pixels(1:npix))
!!$       do i = 0, npix-1
!!$          pixels(i+1) = i
!!$       end do
!!$    end if
!!$
!!$    if (present(np)) np = npix
!!$    
!!$    
!!$  end subroutine allocate_map

  
  subroutine read_map(filename, map)
    implicit none

    character(len=*),                 intent(in)  :: filename
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

    call input_map(filename, map, npix, nmaps, ignore_polcconv=.true.)
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
    if (size(x) == 1) then ! Added to handle delta bandpasses.
       tsum = y(1)
       return
    end if
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

    !write(*,*) sum(abs(x)), sum(abs(y))
    prod = dot_product(x,y)
    call mpi_allreduce(prod, mpi_dot_product, 1, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, comm, ierr)

  end function mpi_dot_product



  ! Routine for reading a spectrum file
  subroutine read_spectrum(filename, spectrum)
    implicit none

    character(len=*),                               intent(in)  :: filename
    real(dp),           allocatable, dimension(:,:)             :: spectrum

    integer(i4b)        :: i, j, numpoint, unit
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
    
    integer(i4b)        :: unit
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

    integer(i4b)  :: i, j, m, n, l, lmax, i_1pct
    real(dp)      :: theta_max, threshold, sigma, theta_1pct
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

    ! Replace beam under 1% of real-space amplitude with a Gaussian extrapolation
    do j = 1, nmaps
       y(:,j) = y(:,j) / maxval(y(:,j))

       i_1pct = 1
       do while (y(i_1pct,j) > 0.01d0 .and. i_1pct < n)
          i_1pct = i_1pct+1
       end do
       theta_1pct = x(i_1pct)
       sigma      = theta_1pct / sqrt(-2.d0*log(0.01d0))

       do i = n, i_1pct, -1
          y(i,j) = exp(-0.5*(x(i)/sigma)**2) * y(i_1pct,j) / exp(-0.5*(x(i_1pct)/sigma)**2)
       end do
    end do
    
    ! Spline significant part of beam profile
    allocate(br(nmaps))
    do j = 1, nmaps
              m      = 0
       do while (y(m+1,j) > threshold .and. m < n)
          m = m+1
       end do
       call spline(br(j), x(1:m), y(1:m,j))
    end do

!!$    open(68,file='beam.dat')
!!$    do j = 1, m
!!$       write(68,*) br(1)%x(j)*180/pi, br(1)%y(j)
!!$    end do
!!$    close(68)
!!$    stop

    deallocate(x, y, pl)
    
  end subroutine compute_radial_beam

  function rand_trunc_gauss(rng_handle, mu, mu_, sigma)
    implicit none

    real(dp),         intent(in)    :: mu, mu_, sigma
    type(planck_rng), intent(inout) :: rng_handle
    real(dp)                        :: rand_trunc_gauss

    real(dp) :: alpha, z, rho, mu0

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
  
  
  subroutine write_map2(filename, map, comptype, nu_ref, unit, ttype, spectrumfile)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map
    character(len=*),                   intent(in), optional :: comptype, unit, spectrumfile, ttype
    real(dp),                           intent(in), optional :: nu_ref

    integer(i4b)   :: npix, nlheader, nmaps, i, nside
    logical(lgt)   :: polarization

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

  end subroutine write_map2


  ! ============================================================================
  ! "WMAP_Read_NInv" reads a WMAP N-Inverse FITS file.
  !
  ! If the output array is  unassociated then this routine will allocate space
  ! for it.  This is why it is declared as a pointer.  If a destination array is
  ! supplied then be sure that it is large enough!
  !
  ! This routine requires the FITSIO library, available from
  ! http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
  !
  ! Arguments:
  !	File      - The name of the file.
  !	Status    - A status code:  0=success.  Most of these are FITS status
  !	            codes, though memory allocation errors are also passed back
  !	            unchanged.
  !	NInv      - The square N-inverse array.
  !	NElements - The number of elements on each side of NInv.  Optional.
  !	NPixels   - The number of pixels in the maps that N-inverse applies to.
  !	            Optional.
  !
  ! Written by Michael R. Greason, SSAI, 13 March 2006.
  ! ============================================================================
  Subroutine WMAP_Read_NInv (File, Status, NInv, ordering, NElements, NPixels)
    !
    Implicit None
    !
    Character (*),                 Intent(In)            :: File
    Integer (Kind=4),              Intent(Out)           :: Status
    Real (Kind=4), Dimension(:,:), Pointer               :: NInv
    Character (80),                Intent(Out), Optional :: ordering
    Integer (Kind=4),              Intent(Out), Optional :: NElements
    Integer (Kind=4),              Intent(Out), Optional :: NPixels
    !
    Character (80)                 :: comm
    Integer (Kind=4), Dimension(2) :: naxes
    Integer (Kind=4)               :: unit, lst, rwmode, tmp, npix
    Logical                        :: anyf
    ! ----------------------------------------------------------------------------
    If (Present(NPixels)  ) NPixels   = 0
    If (Present(NElements)) NElements = 0
    Status  = 0
    !
    ! Open the FITS file.  Leave it positioned at the
    ! primary header/data unit (HDU).
    !
    rwmode = 0
    Call FTGIOU (unit, Status)                    ! Get a free unit #.
    Call FTOPEN (unit, File, rwmode, tmp, Status) ! Open the file.
    If (Status .NE. 0) Return
    !
    !  How big is the array?
    !
    Call FTGISZ (unit, 2, naxes, Status)
    If ((naxes(1) .LE. 0) .OR. (naxes(1) .NE. naxes(2)) .OR. (Status .NE. 0)) Go To 99
    !
    ! Determine whether the maps are in nest ordering
    !
    Call FTGKEY (unit, 'ORDERING', ordering, comm, Status)
    !
    !  How many pixels are in the base map?  Start by looking
    !  at the NPIX keyword; if that isn't found, use LASTPIX.
    !  If neither is found, give up!
    !
    Call FTGKYJ (unit, 'NPIX', npix, comm, Status)
    If (Status .NE. 0) Then
       !
       Status = 0
       Call FTGKYJ (unit, 'LASTPIX', npix, comm, Status)
       If (Status .NE. 0) Go To 99
       npix = npix + 1
       !
    End If
    !
    !  Extract data from this first extension table.
    !
    If (.NOT. Associated(NInv)) Then
       Allocate(NInv(naxes(1), naxes(2)), Stat=Status)
       If (Status .NE. 0) Go To 99
    End If
    !
    tmp = naxes(1) * naxes(2)
    Call FTGPVE (unit, 0, 1, tmp, 0.0E0, NInv, anyf, Status)
    !
    !  Done!  Set the number of pixels, close the FITS file
    !  and return.
    !
    If (Present(NPixels)  ) NPixels   = npix
    If (Present(NElements)) NElements = naxes(1)
    !
99  Continue
    !
    lst = 0
    Call FTCLOS (unit, lst)  ! Close the file.
    Call FTFIOU (unit, lst)  ! Release the unit #.
    If ((lst .NE. 0) .AND. (Status .EQ. 0)) Status = lst
    !
    Return
    ! ----------------------------------------------------------------------------
  End Subroutine WMAP_Read_NInv


  subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
    implicit none

    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix

    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi

    sphi = sin(phi)
    cphi = cos(phi)

    sth  = sin(theta)
    cth  = cos(theta)

    spsi = sin(psi)
    cpsi = cos(psi)

    euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
    euler_matrix(1,3) =                sth * cphi
    euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
    euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(2,3) =                sth * sphi
    euler_matrix(3,1) =              - sth * cpsi
    euler_matrix(3,2) =                sth * spsi
    euler_matrix(3,3) =                cth

  end subroutine compute_euler_matrix_zyz

  subroutine linspace(from, to, array)  ! Hat tip: https://stackoverflow.com/a/57211848/5238625
    implicit none
    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
       array(1) = from
       return
    end if

    do i=1, n
       array(i) = from + range * (i - 1) / (n - 1)
    end do
  end subroutine linspace


  subroutine moving_average(input_data, output_data, window_size, weights, &
        & output_summed_weights)
     implicit none
     real(dp),  dimension(:), intent(in)                :: input_data
     real(dp),  dimension(:), intent(in), optional      :: weights
     real(dp),  dimension(:), intent(out)               :: output_data
     real(dp),  dimension(:), intent(out), optional     :: output_summed_weights
     integer(i4b)           , intent(in)                :: window_size

     integer(i4b)       :: start_ind, i, data_len, end_ind
     real(dp)           :: curr_mean

     data_len = size(input_data)
     output_data = 0.d0
     output_summed_weights = 0.d0

     do i = 1, data_len
        if (input_data(i) == 0.d0) cycle
         start_ind = max(1, int(i - window_size/2))
         end_ind = min(data_len, int(i + window_size/2))
         if (present(weights)) then
            if (sum(weights(start_ind:end_ind)) == 0) then
               curr_mean = 0.d0
            else
               curr_mean = sum(input_data(start_ind:end_ind) * weights(start_ind:end_ind)) /sum(weights(start_ind:end_ind))
            end if
            if (present(output_summed_weights)) then
               output_summed_weights(i) = sum(weights(start_ind:end_ind))
            end if
         else
            curr_mean = sum(input_data(start_ind:end_ind)) / (end_ind - start_ind + 1)
         end if
         output_data(i) = curr_mean
     end do

  end subroutine moving_average

  subroutine moving_average_variable_window(input_data, output_data, &
        & window_sizes, weights, output_summed_weights, kernel_type)
     implicit none
     real(dp),  dimension(:), intent(in)                :: input_data
     real(dp),  dimension(:), intent(in), optional      :: weights
     real(dp),  dimension(:), intent(out)               :: output_data
     real(dp),  dimension(:), intent(out), optional     :: output_summed_weights
     integer(i4b), dimension(:), intent(in)             :: window_sizes
     character(len=128)                                 :: kernel_type

     integer(i4b)       :: i, j, data_len,  max_window_size
     integer(i4b)       :: range_start, range_end, window_size, kernel_start
     integer(i4b)       :: kernel_end
     real(dp)           :: curr_mean
     real(dp), allocatable, dimension(:)        :: kernel

     data_len = size(input_data)
     output_data = 0.d0
     output_summed_weights = 0.d0

     max_window_size = maxval(window_sizes)

     do i = 1, data_len
         if (input_data(i) == 0.d0) cycle
         window_size = window_sizes(i)
         allocate(kernel(window_size + 1))
         if (trim(kernel_type) == 'boxcar') then
            kernel = 1.d0
         else if (trim(kernel_type) == 'gaussian') then
            do j = 1, window_size + 1
               kernel(j) = exp(-0.5 * ((real(j, dp) - 0.5d0 * window_size + 1) / &
                  & (window_size/4.d0)) ** 2)
            end do
         else
            write(*, *) 'Smoothing kernel invalid -- using boxcar'
            kernel = 1.d0
         end if
         range_start = max(1, int(i - window_size/2))
         range_end = min(data_len, int(i + window_size/2))
         kernel_start = window_size/2 + 1 - (i - range_start)
         kernel_end = window_size/2 + 1 + (range_end - i)
         if (present(weights)) then
            if (sum(weights(range_start:range_end)) == 0) then
               curr_mean = 0.d0
            else
               curr_mean = sum(input_data(range_start:range_end) * &
                  & weights(range_start:range_end) * &
                  & kernel(kernel_start:kernel_end)) / &
                  & sum(weights(range_start:range_end) * &
                  & kernel(kernel_start:kernel_end))

!!$               if (i == 3652) then
!!$                  do j = range_start, range_end
!!$                     write(*,*) j, padded_weights(j), padded_data(j)
!!$                  end do
!!$               end if
!!$               write(*,*) 'y', i, range_start, range_end
!!$               write(*,*) 'x', i, lbound(padded_weights),ubound(padded_weights)
!!$               write(*,*) 'r', i, padded_data(range_start:range_end)
!!$               !write(*,*) 's', i, input_data
!!$               if (size(input_data)==3) write(*,*) 'a', i, sum(padded_weights(range_start:range_end))
!!$               if (size(input_data)==3) write(*,*) 'b', i, sum(padded_data(range_start:range_end))
!!$               write(*,*) 'c', i, sum(padded_data(range_start:range_end) * &
!!$                  & padded_weights(range_start:range_end))
!!$               curr_mean = sum(padded_data(range_start:range_end) * &
!!$                  & padded_weights(range_start:range_end)) / &
!!$                  & sum(padded_weights(range_start:range_end))
            end if
            if (present(output_summed_weights)) then
               output_summed_weights(i) = sum(weights(range_start:range_end))
            end if
         else
            curr_mean = sum(input_data(range_start:range_end)*kernel(kernel_start:kernel_end)) / sum(kernel(kernel_start:kernel_end))
         end if
         output_data(i) = curr_mean
         deallocate(kernel)
     end do

  end subroutine moving_average_variable_window


  subroutine moving_variance(input_data, output_data, window_size)
     implicit none
     real(dp),  dimension(:), intent(in)        :: input_data
     real(dp),  dimension(:), intent(out)       :: output_data
     integer(i4b)           , intent(in)        :: window_size


     integer(i4b)       :: start_ind, i, data_len, end_ind
     real(dp)           :: curr_var, curr_mean

     data_len = size(input_data)

     do i = 1, data_len
         start_ind = max(1, int(i - window_size/2))
         end_ind = min(data_len, int(i + window_size/2))
         curr_mean = sum(input_data(start_ind:end_ind)) / (end_ind - start_ind + 1)
         curr_var = sum((input_data(start_ind:end_ind) - curr_mean) ** 2) / (end_ind - start_ind +1)

         output_data(i) = curr_var
     end do

  end subroutine moving_variance

  function masked_variance(data, mask)
   implicit none
   real(dp)                                         :: masked_variance

   real(sp),    dimension(:), intent(in)             :: data
   real(sp),    dimension(:), intent(in)             :: mask

   real(dp)         :: currmean, currvar
   integer(i4b)     :: i, n_unmasked

   n_unmasked = count(mask /= 0)
   if(n_unmasked == 0) then
     masked_variance = 9999999999999d0
     return
   end if
   currmean = sum(data * mask) / n_unmasked
   currvar = 0
   do i = 1, size(data)
      if (mask(i) == 0) cycle
      currvar = currvar + (data(i) - currmean) ** 2
   end do
   currvar = currvar / n_unmasked
   masked_variance = currvar

  end function masked_variance


  !*************************************************
  !    Convert integer to string
  !*************************************************
  character(len=20) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
  end function str


  function calc_corr_coeff(spec_arr, n_samp)
    implicit none
    integer(i4b),               intent(in), optional :: n_samp
    real(dp),     dimension(:), intent(in)           :: spec_arr
    real(dp)                                         :: calc_corr_coeff

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y, covarr
    real(dp), dimension(:), allocatable :: xarr, yarr
    !
    !  Function to calculate the Pearson correlation coefficient between samples in a set of a sampled spectral parameter.
    !
    !  Arguments:
    !  ------------------
    !  spec_arr: double precision array, unknown length
    !     An array containing the sampled spectral parameters
    !  n_spec: integer(i4b)
    !     The last sample number/index
    !  n_lim: integer(i4b)
    !     The maximum index to consider when computing the correlation coefficient
    !  Returns:
    !  -----------------
    !  calc_corr_coeff: double precision real number
    !     The computed correlation coefficient of the sample set.
    !

    !default if no corr coeff is found, return a value that is non-physical.
    !correlation coefficient is -1 <= coeff <= 1, by definition
    calc_corr_coeff = -2.d0 

    ns=min(size(spec_arr),n_samp)

    if (ns > 2) then
       allocate(xarr(ns),yarr(ns))
       do i = 1,ns
          xarr(i)=i*1.d0
          yarr(i)=spec_arr(i)
       end do
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       if (sig_x > 0.d0 .and. sig_y > 0.d0) then
          covarr = sum((xarr-x_mean)*(yarr-y_mean))/ns
          calc_corr_coeff = covarr/(sig_x*sig_y)
       end if
       deallocate(xarr,yarr)
    end if

  end function calc_corr_coeff

  function calc_corr_len(spec_arr,n_samp) 
    implicit none
    integer(i4b),               intent(in), optional :: n_samp
    real(dp),     dimension(:), intent(in)           :: spec_arr
    integer(i4b)                                     :: calc_corr_len

    integer(i4b) :: i, j, maxcorr, ns
    real(dp)     :: x_mean, y_mean, sig_x, sig_y
    real(dp), dimension(:), allocatable :: xarr, yarr,covarr
    !
    !  Function to calculate the sampling correlation length of a spectral parameter sampling.
    !
    !  Arguments:
    !  ------------------
    !  spec_arr: double precision array, unknown length
    !     An array containing the sampled spectralparameters
    !  n_samp: integer(i4b)
    !     the total length of spec_arr of which to compute the correlation length. must be 0 < n_samp <= len(spec_arr).
    !
    !  Returns:
    !  -----------------
    !  calc_corr_len: integer(i4b)
    !     The computed sample correlation length of the sampled spectral parameters
    !

    calc_corr_len = -1 !default if no corr length is found

    maxcorr = min(n_samp/2,size(spec_arr)/2)
    
    allocate(xarr(maxcorr),yarr(maxcorr),covarr(maxcorr))

    do i = 1,maxcorr
       ns=maxcorr-i
       xarr=spec_arr(1:ns)
       yarr=spec_arr(1+i:ns+i)
       x_mean=sum(xarr)/ns
       y_mean=sum(yarr)/ns
       sig_x = sqrt(sum((xarr-x_mean)**2)/ns)
       sig_y = sqrt(sum((yarr-y_mean)**2)/ns)
       covarr(i) = sum((xarr-x_mean)*(yarr-y_mean))/(ns*sig_x*sig_y)
       if (covarr(i) < 0.d0) then
          calc_corr_len = i
          exit
       end if
    end do

    deallocate(xarr,yarr,covarr)

  end function calc_corr_len


   subroutine leggaus(deg, x, w, x1, x2)
      ! Computes the sample points and weights for Gauss-Legendre quadrature.
      !
      ! Given lower and upper integration limits `x1`, and `x2, and `n`, the order of 
      ! quadrature, this routine returns the integration grid/abscissas `x` and 
      ! corresponding weights `w` of the Gauss-Legendre n-point quadrature formula. 
      ! Given the returned grid and weights, the integral of some function f is computed 
      ! as follows: integral = sum(f(x) * w).
      ! 
      ! NOTE: Both x1 and x2 must be present to make use of either.
      !
      ! This code is from NUMERICAL RECIPES in FORTRAN 77 the second edition page 145.
      !
      ! Parameters:
      ! -----------
      ! deg: 
      !     Gaussian quadrature order.    
      ! x: 
      !     1-D ndarray containing the sample points
      ! w: 
      !     1-D ndarray containing the weights.
      ! x1: optional
      !     Lower integration limit.         
      ! x2: optional 
      !     Upper integration limit.      

      implicit none
   
      integer(i4b), intent(in) :: deg
      real(dp), intent(out) :: x(deg), w(deg)
      real(dp), intent(in), optional :: x1, x2

      integer(i4b) :: i, j, m 
      real(dp) :: p1, p2, p3, pp, xl, xm, z, z1, EPS

      EPS = 3.d-14
      m = (deg + 1) / 2
      if (present(x1) .and. present(x2)) then
         xm = 0.5d0 * (x2 + x1)
         xl = 0.5d0 * (x2 - x1)
      else
         xm = 0.d0
         xl = 1.d0
      end if

      do i = 1, m
         z = cos(pi * (i - .25d0)/(deg + .5d0))
      1  continue
            p1 = 1.d0
            p2 = 0.d0
            do j = 1, deg
               p3 = p2
               p2 = p1
               p1 = ((2.d0 * j - 1.d0) * z * p2 - (j-1.d0) * p3) / j
            end do

            pp = deg * (z * p1 - p2) / (z * z - 1.d0)
            z1 = z
            z = z1 - p1/pp
         if (abs(z - z1) .gt. EPS) goto 1
         x(i) = xm - xl * z
         x(deg + 1 - i) = xm + xl * z
         w(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)
         w(deg + 1 - i) = w(i)
      enddo 
   end subroutine leggaus

   subroutine ecl_to_gal_rot_mat(m)
      implicit none
      real(dp), dimension(3,3) :: m

      m(1,1) =  -0.054882486d0
      m(1,2) =  -0.993821033d0
      m(1,3) =  -0.096476249d0
      m(2,1) =   0.494116468d0
      m(2,2) =  -0.110993846d0
      m(2,3) =   0.862281440d0
      m(3,1) =  -0.867661702d0
      m(3,2) =  -0.000346354d0
      m(3,3) =   0.497154957d0
   end subroutine ecl_to_gal_rot_mat

   function get_string_index(arr, str, allow_missing)
     implicit none
     character(len=*), dimension(:), intent(in) :: arr
     character(len=*),               intent(in) :: str
     logical(lgt),                   intent(in), optional :: allow_missing
     integer(i4b)                               :: get_string_index

     integer(i4b) :: i
     character(len=128) :: str1, str2
     logical(lgt) :: allow

     allow = .false.; if (present(allow_missing)) allow = allow_missing

     str1 = str
     call toupper(str1)
     do i = 1, size(arr)
        str2 = arr(i)
        call toupper(str2)
        if (trim(str1) == trim(str2)) then
           get_string_index = i
           exit
        end if
     end do
     if (i > size(arr)) then
        if (allow) then
           get_string_index = -1
        else
           write(*,*) 'get_string_index: String not found = ', trim(str)
           stop
        end if
     end if

   end function get_string_index


   ! from https://rosettacode.org/wiki/Determine_if_a_string_is_numeric#Fortran
   FUNCTION is_numeric(string)
     IMPLICIT NONE
     CHARACTER(len=*), INTENT(IN) :: string
     LOGICAL :: is_numeric
     REAL :: x
     INTEGER :: e
     READ(string,*,IOSTAT=e) x
     is_numeric = e == 0
   END FUNCTION is_numeric

end module comm_utils
