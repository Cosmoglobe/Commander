module comm_output_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_chisq_mod
  use comm_hdf_mod
  implicit none

contains

  subroutine output_FITS_sample(cpar, iter, output_hdf)
    implicit none
    
    type(comm_params), intent(in) :: cpar
    integer(i4b),      intent(in) :: iter
    logical(lgt),      intent(in) :: output_hdf

    integer(i4b)                 :: i, hdferr, ierr
    real(dp)                     :: chisq, t1, t2, t3, t4
    logical(lgt), save           :: first_call=.true.
    logical(lgt)                 :: exist, init
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix, chainfile, hdfpath
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: map
    class(comm_comp),    pointer :: c
    type(hdf_file) :: file
    TYPE(h5o_info_t) :: object_info
    
    call int2string(cpar%mychain, ctext)
    call int2string(iter,         itext)
    postfix = 'c'//ctext//'_k'//itext
    chainfile = trim(adjustl(cpar%outdir)) // '/' // trim(adjustl(cpar%chain_prefix)) // &
         & '_c' // trim(adjustl(ctext)) // '.h5'

    ! Delete existing chain file if necessary; create new file if necessary; open file
    if (first_call .and. cpar%myid_chain == 0 .and. output_hdf) then
       if (trim(cpar%init_chain_prefix) /= trim(cpar%chain_prefix)) then
          inquire(file=trim(chainfile), exist=exist)
          if (exist) call rm(trim(chainfile))
       end if
       inquire(file=trim(chainfile), exist=exist)
       if (.not. exist) then
          call open_hdf_file(chainfile, file, 'w')
          call close_hdf_file(file)
       end if
       first_call = .false.
    end if
    if (cpar%myid_chain == 0 .and. output_hdf) then
       call open_hdf_file(chainfile, file, 'b')
       ! Delete group if it already exists
       call h5eset_auto_f(0, hdferr)
       call h5oget_info_by_name_f(file%filehandle, trim(adjustl(itext)), object_info, hdferr)
       if (hdferr == 0) call h5gunlink_f(file%filehandle, trim(adjustl(itext)), hdferr)
       call create_hdf_group(file, trim(adjustl(itext)))
       call create_hdf_group(file, trim(adjustl(itext))//'/md')             
    end if

    ! Output instrumental parameters
    call output_inst_params(cpar, file, itext, &
         & trim(cpar%outdir)//'/instpar_'//trim(postfix)//'.dat', output_hdf)

    ! Output component results
    c => compList
    call wall_time(t1)
    do while (associated(c))
       call c%dumpFITS(iter, file, output_hdf, postfix, cpar%outdir)
       select type (c)
       class is (comm_diffuse_comp)
          if (output_hdf) then
             hdfpath = trim(adjustl(itext))//'/'//trim(adjustl(c%label))
             call c%Cl%writeFITS(cpar%mychain, iter, hdffile=file, hdfpath=hdfpath)
          else
             call c%Cl%writeFITS(cpar%mychain, iter)
          end if
       end select
       c => c%next()
    end do
    !call wall_time(t2)
    !if (cpar%myid == 0) write(*,*) 'components = ', t2-t1

    ! Output channel-specific residual maps
    if (cpar%output_residuals) then
       !call wall_time(t1)
       do i = 1, numband
          call wall_time(t3)
          map => compute_residual(i)
          call wall_time(t4)
          !if (cpar%myid == 0) write(*,*) 'compute = ', t4-t3
          call wall_time(t3)
          call map%writeFITS(trim(cpar%outdir)//'/res_'//trim(data(i)%label)//'_'// &
               & trim(postfix)//'.fits')
          call wall_time(t4)
          !if (cpar%myid == 0) write(*,*) 'write = ', t4-t3
          call map%dealloc()
       end do
       !call wall_time(t2)
       !if (cpar%myid == 0) write(*,*) 'residuals = ', t2-t1
    end if

    ! Compute and output chi-square
    if (cpar%output_chisq) then
       info => comm_mapinfo(cpar%comm_chain, cpar%nside_chisq, 0, cpar%nmaps_chisq, cpar%pol_chisq)
       map  => comm_map(info)
       call compute_chisq(chisq_map=map, chisq_fullsky=chisq)
       call map%writeFITS(trim(cpar%outdir)//'/chisq_'// trim(postfix) //'.fits')
       call map%dealloc()
       if (cpar%myid == cpar%root) write(*,fmt='(a,i4,a,e16.3)') &
            & 'Chain = ', cpar%mychain, ' -- chisq = ', chisq
    end if

    if (cpar%myid_chain == 0 .and. output_hdf) call close_hdf_file(file)    
  end subroutine output_FITS_sample

  subroutine output_inst_params(cpar, chainfile, iter, filename, output_hdf)
    implicit none
    type(comm_params), intent(in) :: cpar
    type(hdf_file),    intent(in) :: chainfile
    character(len=*),  intent(in) :: iter
    character(len=*),  intent(in) :: filename
    logical(lgt),      intent(in) :: output_hdf
    
    integer(i4b) :: unit, i

    if (cpar%myid /= 0) return

    if (output_hdf) then
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/gain')
       call create_hdf_group(chainfile, trim(adjustl(iter))//'/bandpass')
    end if
    
    unit = getlun()
    open(unit, file=trim(filename), recl=1024)
    write(unit,*) '# Band          Gain         delta_bp'
    do i = 1, numband
       write(unit,'(a15,f10.4,f10.2)') adjustl(trim(data(i)%label)), data(i)%gain, data(i)%bp%delta
       if (output_hdf) then
          call write_hdf(chainfile, trim(adjustl(iter))//'/gain/'//trim(adjustl(data(i)%label)), &
               & data(i)%gain)
          call write_hdf(chainfile, trim(adjustl(iter))//'/bandpass/'//trim(adjustl(data(i)%label)), &
               & data(i)%bp%delta)
       end if
    end do
    close(unit)

  end subroutine output_inst_params

end module comm_output_mod
