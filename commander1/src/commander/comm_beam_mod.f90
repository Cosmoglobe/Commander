module comm_beam_mod
  use healpix_types
  use comm_utils
  implicit none

  integer(i4b),                              private :: lmax, nside, npix, nmaps
  integer(i4b),                              private :: lmax_lowres
  integer(i4b),                              private :: comm_alms, myid_alms, root, numprocs, ierr
  integer(i4b),                              private :: comm_chain, myid_chain
  character(len=128),                        private :: convolution_method

  real(dp),     allocatable, dimension(:,:)  :: beam, beam_lowres


contains

  subroutine initialize_beam_mod(paramfile, comm_alms_in, comm_chain_in, map_id)
    implicit none

    integer(i4b),       intent(in) :: comm_alms_in, comm_chain_in, map_id
    character(len=128), intent(in) :: paramfile

    integer(i4b)       :: unit
    logical(lgt)       :: polarization
    character(len=2)   :: map_text
    character(len=128) :: beamfile, paramtext, pixwinfile

    real(dp) :: fwhm_lowres
    real(dp), pointer, dimension(:,:)   :: pixwin

    unit = getlun()
    comm_chain = comm_chain_in
    call mpi_comm_rank(comm_chain, myid_chain, ierr)
    comm_alms = comm_alms_in
    call mpi_comm_rank(comm_alms, myid_alms, ierr)
    call mpi_comm_size(comm_alms, numprocs, ierr)
    root = 0

    call int2string(map_id, map_text)

    call get_parameter(paramfile, 'NSIDE',                   par_int=nside)       
    call get_parameter(paramfile, 'LMAX',                    par_int=lmax)       
    call get_parameter(paramfile, 'LMAX_LOWRES',             par_int=lmax_lowres)       
    call get_parameter(paramfile, 'FWHM_LOWRES',             par_dp=fwhm_lowres)       
    call get_parameter(paramfile, 'POLARIZATION',            par_lgt=polarization)     
    call get_parameter(paramfile, 'BEAM_CONVOLUTION_METHOD', par_string=convolution_method)     
    call get_parameter(paramfile, 'PIXEL_WINDOW',            par_string=pixwinfile)
    npix = 12*nside**2

   if (polarization) then
       nmaps = 3
    else
       nmaps = 1
    end if

    paramtext = 'BEAM' // map_text
    call get_parameter(paramfile, trim(paramtext), par_string=beamfile)


    ! Read beam transfer function; needed even for effective beam for preconditioning
    allocate(beam(0:lmax, nmaps))
    call read_beam(beamfile, beam)
    
    ! Merge pixel and beam windows into one quantity
    allocate(pixwin(0:lmax,1:nmaps))
    call read_beam(pixwinfile, pixwin)
    beam = beam * pixwin(0:lmax,1:nmaps)
    deallocate(pixwin)
    

    ! Read low-resolution beam transfer function
    allocate(beam_lowres(0:lmax_lowres, nmaps))
    call gaussbeam(fwhm_lowres, lmax_lowres, beam_lowres)       
    
    ! Merge pixel and beam windows into one quantity
    call read_pixwin(nside, nmaps, pixwin)
    beam_lowres = beam_lowres * pixwin(0:lmax_lowres,1:nmaps)
    deallocate(pixwin)
    

    if (trim(convolution_method) == 'transfer_function') then


    else if (trim(convolution_method) == 'effective_beam') then


    else

       write(*,*) 'comm_beam_mod error: Unknown beam convolution method = ', trim(convolution_method)
       stop

    end if

  end subroutine initialize_beam_mod


  subroutine cleanup_beam_mod
    implicit none

    if (allocated(beam)) deallocate(beam)
    
  end subroutine cleanup_beam_mod



  subroutine multiply_with_beam(alms, deconvolve, use_transfer, transpose)
    implicit none

    real(dp),     dimension(1:,1:), intent(inout)          :: alms
    logical(lgt),                   intent(in),   optional :: deconvolve, use_transfer, transpose

    logical(lgt) :: deconv, transfer
    integer(i4b) :: j, l, m

    ! Check whether we want to use the transfer function approach or the effective beam approach
    transfer = (trim(convolution_method) == 'transfer_function')
    if (present(use_transfer)) then
       transfer = use_transfer .or. transfer
    end if

    ! Check if we want to deconvolve the beam. If so, we have to use the transfer function method
    if (present(deconvolve)) then
       deconv = deconvolve
       transfer  = transfer .or. deconvolve
    else
       deconv = .false.
    end if


    ! Do the convolution
    if (transfer) then
       j = 1
       do l = 0, lmax
          do m = -l, l
             if (deconv) then
                alms(j,:) = alms(j,:) / beam(l,:)
             else
                alms(j,:) = alms(j,:) * beam(l,:)
             end if
             j = j+1
          end do
       end do
    else 

    end if

  end subroutine multiply_with_beam



  subroutine multiply_with_beam_lowres(alms, transpose)
    implicit none

    logical(lgt),                   intent(in),   optional :: transpose
    real(dp),     dimension(0:,1:), intent(inout)          :: alms

    integer(i4b) :: j, l, m

    if (trim(convolution_method) == 'transfer_function') then

       j = 1
       do l = 0, lmax
          do m = -l, l
             alms(j,:) = alms(j,:) * beam_lowres(l,:)
             j = j+1
          end do
       end do

    else if (trim(convolution_method) == 'effective_beam') then

    end if

  end subroutine multiply_with_beam_lowres
  

  function get_b_l(l,spec)
    implicit none

    integer(i4b), intent(in) :: l, spec
    real(dp)                 :: get_b_l

    get_b_l = beam(l,spec)

  end function get_b_l


end module comm_beam_mod
