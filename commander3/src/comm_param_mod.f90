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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Commander parameter module !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module comm_param_mod
  use comm_utils
  use hashtbl
  implicit none

  ! Note: This module reads in the Commander parameter file as the first operation
  !       in the program. This is primarily intended to avoid crashes after hours
  !       of running because of user errors; catch these early, and report back
  !       problems. Then, copy parameters over to module structures if convenient
  !       at a later stage. 

  integer(i4b), parameter, private :: MAXPAR           = 10
  integer(i4b), parameter, private :: MAXAUXPAR        = 10
  integer(i4b), parameter, private :: MAXSAMPGROUP     = 100
  integer(i4b), parameter, private :: MAXZODICOMPS     = 100
  integer(i4b), parameter, private :: MAXZODIPARAMS    = 100


  type comm_params

     ! MPI info
     integer(i4b) :: myid, numprocs, root = 0
     integer(i4b) :: myid_chain, numprocs_chain, comm_chain, mychain
     integer(i4b) :: myid_shared, comm_shared, myid_inter, comm_inter
     integer(i4b), dimension(MPI_STATUS_SIZE)          :: status

     ! Global parameters
     character(len=24)  :: operation
     logical(lgt)       :: resamp_CMB
     integer(i4b)       :: first_samp_resamp, last_samp_resamp, numsamp_per_resamp
     integer(i4b)       :: verbosity, base_seed, base_seed_noise, numchain, num_smooth_scales
     integer(i4b)       :: num_gibbs_iter, thinning, num_init_chains
     character(len=512) :: chain_status, init_chain_prefix
     real(dp)           :: T_CMB
     character(len=512) :: MJysr_convention
     character(len=512) :: fft_magic_number_file
     character(len=512) :: output_comps
     logical(lgt)       :: only_pol, only_I
     logical(lgt)       :: enable_TOD_analysis
     logical(lgt)       :: enable_TOD_simulations !< start commander in simulation regime
     integer(i4b)       :: tod_freq
     integer(i4b)       :: resamp_hard_gain_prior_nth_iter
     integer(i4b)       :: output_4D_map_nth_iter, output_aux_maps
     logical(lgt)       :: include_tod_zodi, sample_zodi, incl_zodi_solar_comp
     integer(i4b)       :: zodi_solar_nside
     character(len=512) :: zodi_solar_initmap, zodi_static_bands
     character(len=512), allocatable, dimension(:)     :: init_chain_prefixes
     character(len=512)                                :: sims_output_dir !< simulations directory

     ! alm-sampler
     integer(i4b)       :: almsamp_nsamp, almsamp_nside_chisq_lowres, almsamp_prior_fwhm, almsamp_burnin
     logical(lgt)       :: almsamp_optimize, almsamp_apply_prior, almsamp_pixreg, almsamp_priorsamp_frozen

     ! Output parameters
     character(len=512) :: outdir
     integer(i4b)       :: nside_chisq, nmaps_chisq
     logical(lgt)       :: pol_chisq, output_mixmat, output_residuals, output_chisq, output_cg_eigenvals
     integer(i4b)       :: output_cg_freq
     logical(lgt)       :: output_input_model, ignore_gain_bp, output_debug_seds, output_sig_per_band
     logical(lgt)       :: sample_signal_amplitudes, sample_specind, sample_powspec

     ! Numerical parameters
     character(len=512) :: cg_conv_crit, cg_precond
     integer(i4b)       :: cg_lmax_precond, cg_maxiter, cg_num_samp_groups, cg_num_user_samp_groups, cg_miniter, cg_check_conv_freq, cg_samp_group_md
     logical(lgt)       :: cg_init_zero, set_noise_to_mean
     real(dp)           :: cg_tol
     integer(i4b)       :: num_bp_prop
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group_mask
     integer(i4b),       dimension(MAXSAMPGROUP) :: cg_samp_group_maxiter
     character(len=512), dimension(MAXSAMPGROUP) :: cg_samp_group_bands

     ! Data parameters
     integer(i4b)       :: numband
     character(len=512) :: datadir, ds_sourcemask, ds_procmask

     character(len=512), allocatable, dimension(:)   :: ds_maskfile
     character(len=512), allocatable, dimension(:)   :: ds_noisefile
     character(len=512), allocatable, dimension(:)   :: ds_instlabel
     logical(lgt),       allocatable, dimension(:)   :: ds_active
     logical(lgt),       allocatable, dimension(:)   :: ds_polarization
     integer(i4b),       allocatable, dimension(:)   :: ds_nside
     integer(i4b),       allocatable, dimension(:)   :: ds_lmax
     character(len=512), allocatable, dimension(:)   :: ds_label


  end type comm_params


contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************
  subroutine read_comm_params(cpar)
    implicit none
    type(hash_tbl_sll) :: htable
    type(comm_params), intent(inout) :: cpar

    integer(i4b)       :: paramfile_len, ierr, i, idx
    character(len=512) :: paramfile, paramfile_name
    character(len=512), allocatable, dimension(:) :: paramfile_cache

    call getarg(1, paramfile)

    write(*,*) cpar%myid, cpar%root, 'stuff is here'

    ! Read parameters into cache
    if (cpar%myid == cpar%root) then
       paramfile_len = 512
       allocate(paramfile_cache(paramfile_len))
       call read_paramfile_to_ascii(paramfile,paramfile_cache,paramfile_len)
    end if
    call mpi_bcast(paramfile_len, 1, MPI_INTEGER, cpar%root, MPI_COMM_WORLD, ierr)
    if (cpar%myid /= cpar%root) then 
       allocate(paramfile_cache(paramfile_len))
    end if
    do i=1,paramfile_len
       call mpi_bcast(paramfile_cache(i), 512, MPI_CHAR, cpar%root, MPI_COMM_WORLD, ierr)
    end do
    call init_hash_tbl_sll(htable,tbl_len=10*paramfile_len)
    call put_ascii_into_hashtable(paramfile_cache,htable)

    ! read the data directory and store it so it can be appended to other paths
    call get_parameter_hashtable(htable, 'DATA_DIRECTORY', par_string=cpar%datadir)


    ! Read parameters from the hash table
    call read_global_params_hash(htable,cpar)
    call read_data_params_hash(htable,cpar)

    !output parameter file to output directory
    if (cpar%myid == cpar%root) then
       idx = index(paramfile, '/', back=.true.) 
       paramfile_name = trim(cpar%outdir)//'/'//paramfile(idx+1:len(paramfile))
    end if
    !deallocate ascii cache
    deallocate(paramfile_cache)

    !Deallocate hash table
    call free_hash_tbl_sll(htable)
  end subroutine read_comm_params

  subroutine initialize_mpi_struct(cpar, handle)
    implicit none
    type(comm_params), intent(inout) :: cpar
    type(planck_rng),  intent(out)   :: handle

    integer(i4b) :: i, j, m, n, ierr, mpistat(MPI_STATUS_SIZE)
    integer(i4b), allocatable, dimension(:,:) :: ind

    ierr = 0

    cpar%numchain = min(cpar%numchain, cpar%numprocs)

    allocate(ind(0:cpar%numprocs-1,2))
    n = 0
    do i = 1, cpar%numchain
       m = cpar%numprocs / cpar%numchain
       if ((cpar%numprocs-(cpar%numprocs/cpar%numchain)*cpar%numchain) >= i) m = m+1
       ind(n:n+m-1,1) = i
       do j = 0, m-1
          ind(n+j,2) = j
       end do
       n = n+m
    end do

    cpar%mychain    = ind(cpar%myid,1)
    cpar%myid_chain = ind(cpar%myid,2)
    cpar%init_chain_prefix = cpar%init_chain_prefixes(mod(cpar%mychain-1,cpar%num_init_chains)+1)

    write(*,*) 'pre-comm_split', cpar%mychain, cpar%myid_chain, cpar%comm_chain, ierr
    call mpi_comm_split(MPI_COMM_WORLD, cpar%mychain, cpar%myid_chain, cpar%comm_chain,  ierr) 
    write(*,*) 'post-comm_split', cpar%mychain, cpar%myid_chain, cpar%comm_chain, ierr
    call mpi_comm_size(cpar%comm_chain, cpar%numprocs_chain, ierr)

    !Communicators for shared memory access
    call mpi_comm_split_type(cpar%comm_chain, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, cpar%comm_shared, ierr) 
    call mpi_comm_rank(cpar%comm_shared, cpar%myid_shared, ierr)
    call mpi_comm_split(cpar%comm_chain, cpar%myid_shared, 0, cpar%comm_inter, ierr)
    call mpi_comm_rank(cpar%comm_inter, cpar%myid_inter, ierr)

    deallocate(ind)

    ! Initialize random number generator
    if (cpar%myid == cpar%root) then
       call rand_init(handle, cpar%base_seed)
       do i = 1, cpar%numprocs-1
          j = nint(rand_uni(handle)*1000000.d0)
          call mpi_send(j, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
       end do
    else 
       call mpi_recv(j, 1, MPI_INTEGER, cpar%root, 98, MPI_COMM_WORLD, mpistat, ierr)
       call rand_init(handle, j)
    end if

  end subroutine initialize_mpi_struct

  ! ********************************************************
  !              Specialized routines; one per module
  ! ********************************************************

  subroutine read_global_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: i
    character(len=2) :: itext


    call get_parameter_hashtable(htbl, 'NUM_INIT_CHAINS',          par_int=cpar%num_init_chains)
    allocate(cpar%init_chain_prefixes(cpar%num_init_chains))
    do i = 1, cpar%num_init_chains
       call int2string(i,itext)
       call get_parameter_hashtable(htbl, 'INIT_CHAIN'//itext,     par_string=cpar%init_chain_prefixes(i))
    end do

  end subroutine read_global_params_hash


  subroutine read_data_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: i, j, n,len_itext
    character(len=3) :: itext
    character(len=2) :: jtext

    call int2string(1, itext)
    len_itext=len(trim(itext))
    call get_parameter_hashtable(htbl, 'NUMBAND',             par_int=cpar%numband)

    n = cpar%numband

    allocate(cpar%ds_polarization(n), cpar%ds_nside(n), cpar%ds_lmax(n))
    allocate(cpar%ds_active(n), cpar%ds_label(n), cpar%ds_instlabel(n))
    allocate(cpar%ds_noisefile(n), cpar%ds_maskfile(n))
    do i = 1, n
       call int2string(i, itext)
       call get_parameter_hashtable(htbl, 'INCLUDE_BAND'//itext, len_itext=len_itext, par_lgt=cpar%ds_active(i))
       if (.not. cpar%ds_active(i)) cycle
       call get_parameter_hashtable(htbl, 'BAND_NSIDE'//itext, len_itext=len_itext, par_int=cpar%ds_nside(i))
       call get_parameter_hashtable(htbl, 'BAND_LMAX'//itext, len_itext=len_itext, par_int=cpar%ds_lmax(i))
       call get_parameter_hashtable(htbl, 'BAND_NOISEFILE'//itext, len_itext=len_itext, par_string=cpar%ds_noisefile(i), path=.true.)
    end do


  end subroutine read_data_params_hash


  ! ********************************************************
  !                     Utility routines
  ! ********************************************************


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

    integer(i4b) :: n

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
          write(*,*) "parse_parameter: Reached unreachable point! ", present(par_string)
       end if
       found = .true.
    else
       found = .false.
    end if
  end subroutine parse_parameter

  !gets parameters from input arguments in Commander call
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
    do i = 1, command_argument_count() !iargc()
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
       write(*,*) "get_parameter_arg: Fatal error: Cannot find " // trim(parname) // " in argument list!"
       if(present(desc)) write(*,*) trim(desc)
       stop
    end if
  end subroutine get_parameter_arg

  ! Loops through the parameter files and children, counting lines.
  ! No error reporting.

  subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int
  end subroutine str2int

  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = string(ext(1):ext(2))
  end function get_token

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
  end subroutine get_tokens

  subroutine get_detectors(filename, detectors, num_dets)
    !
    ! Reads detector names from a text file and saves them in a character array.
    !
    ! Arguments:
    ! ----------
    ! filename:  character string
    !            Filename of the file where detector names are stored.
    ! num_dets:  integer (optional)
    !            Number of detectors
    !
    ! Return:
    ! -------
    ! detectors: character array
    !            Initially empty array is filled with detector names. 
    ! 
    implicit none
    character(len=*), intent(in)           :: filename
    character(len=*), intent(inout)        :: detectors(:)
    integer(i4b),     intent(in), optional :: num_dets

    character(len=500)           :: detector_list_file
    integer(i4b)                 :: unit,io_error,counter, ndet, i
    character(len=30)             :: line

  end subroutine get_detectors

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
  end function has_token

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
  end function num_tokens

  integer(i4b) function count_detectors(filename)
    ! 
    ! Takes in the filename and directory of a detector list and returns the number of 
    ! detectors in that list. Each detector has to be written on a separate line, as 
    ! the function simply counts the lines of the file that don't start in '#'.
    !
    ! Arguments:
    ! ----------
    ! filename:    character string
    !              Filename of the detector list             
    !
    ! Returns:
    ! --------
    ! count_detectors: integer
    !                  Number of lines in the file that are not commented out using '#'.
    !
    implicit none
    character(len=*) :: filename

    character(len=500)           :: detector_list_file
    integer(i4b)                 :: unit,io_error,counter
    logical                      :: counting
    character(len=8)             :: line
    count_detectors=0

  end function count_detectors

  subroutine tokenize(string, sep, ext, group, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional   :: group
    character(len=256)  :: op, cl
    integer(i4b), save           :: level(256), nl
    integer(i4b), intent(inout)  :: ext(2)
    logical(lgt), optional       :: allow_empty

    integer(i4b) :: i, j, o, c, ng
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
    j     = 1
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
  end subroutine tokenize


  subroutine read_paramfile_to_ascii(paramfile,paramfile_cache, paramfile_len)
    implicit none
    character(len=512),                            intent(in)  :: paramfile
    character(len=512), allocatable, dimension(:), intent(inout) :: paramfile_cache
    integer(i4b),intent(inout) :: paramfile_len

    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), line_nr,i, stat, pos
    character(len=512)         :: key, value, filenames(maxdepth), line
    character(len=1024)        :: default_path
    character(len=3)           :: band_num

    character(len=512), allocatable, dimension(:) :: new_cache

    ! read file to ascii array


    call get_environment_variable("COMMANDER_PARAMS_DEFAULT", default_path, status=stat)

    band_num = 'XXX'
    line_nr = 0
    depth = 1
    units(depth) = getlun()
    filenames(depth) = paramfile
    open(units(depth),file=trim(paramfile),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),fmt="(a)",end=1) line
             pos = index(line, ' ')
             value = adjustl(line(pos:len(line)))
             pos = index(value, ' ')
             value = trim(value(1:pos))

             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else if(key == '@DEFAULT') then
             if(stat /= 0) then
                write(*,*) "Paramater file uses @DEFAULT command but the environment variable COMMANDER_PARAMS_DEFAULT returns ", stat
                stop
             end if
             ! Recurse to the default new file
             read(units(depth),fmt="(a)",end=1) line
             pos = index(line, ' ')
             value = adjustl(line(pos:len(line)))
             pos = index(value, ' ')
             value = trim(value(1:pos))

             depth = depth+1
             units(depth) = getlun()
             filenames(depth) = trim(default_path)//'/'//trim(value)
             open(units(depth),file=filenames(depth), status="old", err=2)
          else if(key == '@START') then
             read(units(depth),*,end=1) key, value
             if(band_num /= 'XXX') then
                write(*,*) "Error starting band number ", trim(value), ", band ", band_num, " has not ended in file ", trim(filenames(depth))
                stop
             end if
             band_num = value(1:len(value))
          else if(key == '@END') then
             read(units(depth),*,end=1) key, value
             if(value(1:len(value)) /= band_num) then
                write(*,*) "Error ending band ", trim(value), ", current band is ", band_num, " in file ", trim(filenames(depth))
                stop
             end if
             band_num = 'XXX'

          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          !if we get here we have read a new line from the parameter file(s)
          line_nr = line_nr + 1
          if(line_nr > paramfile_len) then !we need to resize the cache array
             allocate(new_cache(2*paramfile_len))

             new_cache(1:paramfile_len) = paramfile_cache(1:paramfile_len)
             deallocate(paramfile_cache)
             call move_alloc(new_cache, paramfile_cache)

             paramfile_len = paramfile_len*2
          end if
          if(band_num /= 'XXX') then !active @START directive
             !replace the string &&& with the band number given by @START
             pos = index(line, '&&&') !if this is a band
             if(pos > 0)  line(pos:pos+2)=band_num
             pos = index(line, '&&') !this could be a component
             if(pos > 0) line(pos:pos+1)=band_num
          else
             pos = index(line, '&&') !check for the special chars outside START
             if(pos > 0) then
                write(*,*) "Warning: parameter line ", line, " found outside of a START-END block"
             end if
          end if
          write(paramfile_cache(line_nr),fmt="(a)") line

       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do
    !resize the cache to be the exact length
    allocate(new_cache(line_nr))

    new_cache(1:line_nr) = paramfile_cache(1:line_nr)
    deallocate(paramfile_cache)

    call move_alloc(new_cache, paramfile_cache)
    paramfile_len = line_nr
    return

    ! ===== Error handling section ======

    ! Case 1: Include file error
2   write(*,*) "Error: Cannot open include file '" // trim(filenames(depth)) // "'"
    write(*,*) " in file " // trim(filenames(depth-1))
    do i = depth-2, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    stop

    ! Case 2: Directive error
3   write(*,*) "Error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    stop

    ! Case 3: Top level parameter file unreadable
4   write(*,*) "Error: Cannot open parameter file '" // trim(paramfile) // "'"
    stop

    end subroutine read_paramfile_to_ascii


    ! read parameter from input argument or hash table
    subroutine get_parameter_hashtable(htbl, parname, len_itext, par_int, par_char, &
         & par_string, par_sp, par_dp, par_lgt, par_present, desc, path)
      implicit none
      type(hash_tbl_sll), intent(in) :: htbl 
      character(len=*),   intent(in) :: parname
      integer(i4b),     optional :: len_itext
      integer(i4b),     optional :: par_int
      character(len=*), optional :: par_char
      character(len=*), optional :: par_string
      real(sp),         optional :: par_sp
      real(dp),         optional :: par_dp
      logical(lgt),     optional :: par_lgt
      logical(lgt),     optional :: par_present
      character(len=*), optional :: desc
      logical(lgt),     optional :: path

      logical(lgt)               :: found

      found = .false.
      call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
      if(found) then
         if(present(par_present)) par_present = .true.
      else
         call get_parameter_from_hash(htbl, parname, len_itext, par_int, &
              & par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc, path)
      end if
    end subroutine get_parameter_hashtable

    ! getting parameter value from hash table
    subroutine get_parameter_from_hash(htbl, parname, len_itext, par_int, par_char, &
         & par_string, par_sp, par_dp, par_lgt, par_present, desc, path)
      implicit none
      type(hash_tbl_sll), intent(in) :: htbl
      character(len=*),   intent(in) :: parname
      integer(i4b),     optional :: len_itext
      integer(i4b),     optional :: par_int
      character(len=*), optional :: par_char
      character(len=*), optional :: par_string
      real(sp),         optional :: par_sp
      real(dp),         optional :: par_dp
      logical(lgt),     optional :: par_lgt
      logical(lgt),     optional :: par_present
      character(len=*), optional :: desc
      logical(lgt),     optional :: path
      character(len=256)         :: key
      character(len=:), ALLOCATABLE   :: itext,jtext
      CHARACTER(len=:), ALLOCATABLE   :: val,val2,val3
      character(len=512)              :: val4
      integer(i4b)                    :: i,j
      logical(lgt)                    :: loc_path
    
      if(.not. present(path)) then 
        loc_path = .false.
      else
        loc_path = path
      end if
      key=trim(parname)
      call tolower(key)
      call get_hash_tbl_sll(htbl,trim(key),val)
      if (.not. allocated(val)) then
         goto 1
         if (.not. present(len_itext)) goto 1
         allocate(character(len=len_itext) :: itext,jtext)
         itext=key(len(trim(key))-(len_itext-1):len(trim(key)))
         call get_hash_tbl_sll(htbl,'band_default_params'//trim(itext),val2)
         if (allocated(val2)) then
            read(val2,*) j
            if (j /= 0) then
               call int2string(j, jtext)
               call get_hash_tbl_sll(htbl,'band_default_params'//trim(jtext),val3)
               if (allocated(val3)) then
                  read(val3,*) i
                  if (i /= 0) goto 2
               end if
               call get_hash_tbl_sll(htbl,key(1:len(trim(key))-len_itext)//trim(jtext),val)
               if (.not. allocated(val)) goto 3
            else
               goto 1
            end if
         else
            goto 1
         end if
         deallocate(itext,jtext)
      end if

      if (val == '#') then
        write(*,*) trim(parname), ' has invalid value #, double check your parameter file'
        stop
      end if

      if (present(par_int)) then
         read(val,*) par_int
      elseif (present(par_char)) then
         read(val,*) par_char
      elseif (present(par_string)) then
         !append data directory if required
         if(len(val) > 0) then
           if(loc_path .and. trim(val) /= 'fullsky' .and. trim(val) /= 'none' .and. trim(val) /= 'native' .and. trim(val) /= 'default') then
             if(val(1:1) /= '/') then
              call get_parameter_hashtable(htbl, "DATA_DIRECTORY", par_string=val4, path=.false.)
              val = trim(val4) // '/' // trim(val)
             end if
           end if
         end if
         !read(val,*) par_string
         par_string = val
      elseif (present(par_sp)) then
         read(val,*) par_sp
      elseif (present(par_dp)) then
         read(val,*) par_dp
      elseif (present(par_lgt)) then
         if (trim(val) == '.true.' .or. trim(val) == '.false.') then
            read(val,*) par_lgt
         else
            write(*,*) "Error: parameter "//trim(parname)//" should be .true. or .false."
            stop
         end if
      else
         write(*,*) "get_parameter: Reached unreachable point! ", val, present(par_string)
      end if

      deallocate(val)
      return

1     write(*,*) "Error: Could not find parameter '" // trim(parname) // "'"
      write(*,*) ""
      stop


2     write(*,*) "Error: Recursive default parameters, bands " // &
           & trim(jtext) // " and " //trim(itext)
      write(*,*) ""
      stop

3     write(*,*) "Error: Could not find parameter '" // trim(parname)//&
         & "' from default '"//key(1:len(trim(key))-len_itext)//trim(jtext)//"'"
      write(*,*) ""
      stop
  end subroutine get_parameter_from_hash

  ! filling the hash table with elements from the parameter file (ascii array) 
  subroutine put_ascii_into_hashtable(asciitbl,htbl)
    implicit none
    character(len=512), allocatable, dimension(:), intent(in) :: asciitbl
    type(hash_tbl_sll), intent(inout) :: htbl
    character(len=512) :: key, val
    character(len=256) :: toks(2)
    integer            :: i, n
    do i = 1,size(asciitbl)
       call get_tokens(trim(asciitbl(i)), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
       if(n < 2) then ! only need the lines where one has 'key'='value'
          cycle
       end if
       key = get_token(toks(1), " ", 1, group="''" // '""')
       val = get_token(toks(2), " ", 1, group="''" // '""')
       call tolower(key)  ! we don't differentiate btw. upper and lower case
       if (key=="") cycle !we don't need blank lines
       call put_hash_tbl_sll(htbl,trim(key),trim(val)) 
    end do
    return
    
    write(*,*) "Error: Cannot read ascii line:", i, "line = '" // trim(asciitbl(i)) // "'"
    stop
    
  end subroutine put_ascii_into_hashtable
  
  subroutine get_chainfile_and_samp(string, chainfile, initsamp)
    implicit none
    character(len=*),   intent(in)  :: string
    character(len=512), intent(out) :: chainfile
    integer(i4b),       intent(out) :: initsamp
    
    integer(i4b) :: i, num, e
    character(len=512), dimension(2) :: toks


    call get_tokens(string, ":", toks, num)    
    chainfile = toks(1)
    read(toks(2),*, iostat=e) initsamp
    if (e .ne. 0) then
      write(*,*) 'Issue with chain file formatting, got ', initsamp, trim(toks(2))
    end if

    if (index(chainfile, '.h5') == 0) then
        write(*,*) "poorly formatted naming of chain file", trim(string)
        write(*,*) "Should be filename:sample, e.g., data/chain_c0001.h5:12"
        stop
    end if
    
  end subroutine get_chainfile_and_samp
  
  subroutine define_cg_samp_groups(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar
    
    integer(i4b) :: i, j, k, n
    character(len=16), dimension(1000) :: comp_label
    
    
  end subroutine define_cg_samp_groups
  

  
end module comm_param_mod
