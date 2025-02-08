
module comm_param_mod
  use comm_utils
  use hashtbl
  implicit none

  ! Note: This module reads in the Commander parameter file as the first operation
  !       in the program. This is primarily intended to avoid crashes after hours
  !       of running because of user errors; catch these early, and report back
  !       problems. Then, copy parameters over to module structures if convenient
  !       at a later stage. 


  type comm_params

     ! MPI info
     integer(i4b) :: myid, numprocs, root = 0
     integer(i4b) :: myid_chain, comm_chain, mychain

     ! Data parameters

     integer(i4b),       allocatable, dimension(:)   :: ds_nside
     integer(i4b),       allocatable, dimension(:)   :: ds_lmax


  end type comm_params


contains

  ! ********************************************************
  !                     Driver routines
  ! ********************************************************
  subroutine read_comm_params(cpar)
    implicit none
    type(hash_tbl_sll) :: htable
    type(comm_params), intent(inout) :: cpar

    integer(i4b)       :: paramfile_len, ierr, i
    character(len=512) :: paramfile
    character(len=512), allocatable, dimension(:) :: paramfile_cache

    call getarg(1, paramfile)

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



    ! Read parameters from the hash table
    call read_data_params_hash(htable,cpar)

    deallocate(paramfile_cache)

    !Deallocate hash table
    call free_hash_tbl_sll(htable)
  end subroutine read_comm_params

  subroutine initialize_mpi_struct(cpar)
    implicit none
    type(comm_params), intent(inout) :: cpar

    integer(i4b) :: i, j, m, n, ierr
    integer(i4b), allocatable, dimension(:,:) :: ind

    ierr = 0

    allocate(ind(0:cpar%numprocs-1,2))
    n = 0
    do i = 1, 1
       m = cpar%numprocs / 1
       if ((cpar%numprocs-(cpar%numprocs/1)*1) >= i) m = m+1
       ind(n:n+m-1,1) = i
       do j = 0, m-1
          ind(n+j,2) = j
       end do
       n = n+m
    end do

    cpar%mychain    = ind(cpar%myid,1)
    cpar%myid_chain = ind(cpar%myid,2)

    call mpi_comm_split(MPI_COMM_WORLD, cpar%mychain, cpar%myid_chain, cpar%comm_chain,  ierr) 

    deallocate(ind)


  end subroutine initialize_mpi_struct

  ! ********************************************************
  !              Specialized routines; one per module
  ! ********************************************************

  subroutine read_data_params_hash(htbl, cpar)
    implicit none

    type(hash_tbl_sll), intent(in) :: htbl
    type(comm_params),  intent(inout) :: cpar

    integer(i4b)     :: len_itext
    character(len=3) :: itext

    call int2string(1, itext)
    len_itext=len(trim(itext))

    allocate(cpar%ds_nside(1), cpar%ds_lmax(1))

    call get_parameter_hashtable(htbl, 'BAND_NSIDE'//itext, len_itext=len_itext, par_int=cpar%ds_nside(1))
    call get_parameter_hashtable(htbl, 'BAND_LMAX'//itext, len_itext=len_itext, par_int=cpar%ds_lmax(1))


  end subroutine read_data_params_hash


  ! ********************************************************
  !                     Utility routines
  ! ********************************************************


  subroutine parse_parameter(line, parname, found, par_int, par_string)
    implicit none
    character(len=*)           :: line, parname
    character(len=256)         :: toks(2), key, value, par
    logical(lgt)               :: found
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_string

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
       elseif (present(par_string)) then
          read(value,*) par_string
       else
          write(*,*) "parse_parameter: Reached unreachable point! ", present(par_string)
       end if
       found = .true.
    else
       found = .false.
    end if
  end subroutine parse_parameter

  !gets parameters from input arguments in Commander call
  subroutine get_parameter_arg(parname, par_int,  &
       & par_string, par_present, desc)
    implicit none
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_string
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    character(len=512) :: line
    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, command_argument_count() !iargc()
       call getarg(i, line)
       if(line(1:2) /= "--") cycle
       call parse_parameter(line(3:), parname, found, par_int, par_string)
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
    integer(i4b)               :: depth, units(maxdepth), line_nr, stat
    character(len=512)         :: key, filenames(maxdepth), line
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
    open(units(depth),file=trim(paramfile),status="old")
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

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
       write(paramfile_cache(line_nr),fmt="(a)") line

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

    end subroutine read_paramfile_to_ascii


    ! read parameter from input argument or hash table
    subroutine get_parameter_hashtable(htbl, parname, len_itext, par_int,  &
         & par_string, par_present, desc, path)
      implicit none
      type(hash_tbl_sll), intent(in) :: htbl 
      character(len=*),   intent(in) :: parname
      integer(i4b),     optional :: len_itext
      integer(i4b),     optional :: par_int
      character(len=*), optional :: par_string
      logical(lgt),     optional :: par_present
      character(len=*), optional :: desc
      logical(lgt),     optional :: path

      logical(lgt)               :: found

      found = .false.
      call get_parameter_arg(parname, par_int, par_string, found, desc)
      if(found) then
         if(present(par_present)) par_present = .true.
      else
         call get_parameter_from_hash(htbl, parname, len_itext, par_int, &
              & par_string, path)
      end if
    end subroutine get_parameter_hashtable

    ! getting parameter value from hash table
    subroutine get_parameter_from_hash(htbl, parname, len_itext, par_int,  &
         & par_string, path)
      implicit none
      type(hash_tbl_sll), intent(in) :: htbl
      character(len=*),   intent(in) :: parname
      integer(i4b),     optional :: len_itext
      integer(i4b),     optional :: par_int
      character(len=*), optional :: par_string
      logical(lgt),     optional :: path
      character(len=256)         :: key
      character(len=:), ALLOCATABLE   :: itext,jtext
      CHARACTER(len=:), ALLOCATABLE   :: val,val2,val3
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
  
end module comm_param_mod
