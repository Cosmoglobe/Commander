module comm_output_mod
  use comm_param_mod
  use comm_comp_mod
  use comm_chisq_mod
  implicit none

contains

  subroutine output_FITS_sample(cpar, iter)
    implicit none
    
    type(comm_params), intent(in) :: cpar
    integer(i4b),      intent(in) :: iter

    integer(i4b)                 :: i
    real(dp)                     :: chisq, t1, t2, t3, t4
    character(len=4)             :: ctext
    character(len=6)             :: itext
    character(len=512)           :: postfix
    class(comm_mapinfo), pointer :: info
    class(comm_map),     pointer :: map
    class(comm_comp),    pointer :: c

    call int2string(cpar%mychain, ctext)
    call int2string(iter,         itext)
    postfix = 'c'//ctext//'_k'//itext

    ! Output instrumental parameters
    !call wall_time(t1)
    call output_inst_params(cpar, trim(cpar%outdir)//'/instpar_'//trim(postfix)//'.dat')    
    !call wall_time(t2)
    !if (cpar%myid == 0) write(*,*) 'instrument = ', t2-t1

    ! Output component results
    c => compList
    call wall_time(t1)
    do while (associated(c))
       call c%dumpFITS(postfix, cpar%outdir)
       select type (c)
       class is (comm_diffuse_comp)
          call c%Cl%writeFITS(cpar%mychain, iter)
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
          if (cpar%myid == 0) write(*,*) 'compute = ', t4-t3
          call wall_time(t3)
          call map%writeFITS(trim(cpar%outdir)//'/res_'//trim(data(i)%label)//'_'// &
               & trim(postfix)//'.fits')
          call wall_time(t4)
          if (cpar%myid == 0) write(*,*) 'write = ', t4-t3
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
    
  end subroutine output_FITS_sample

  subroutine output_inst_params(cpar, filename)
    implicit none
    type(comm_params), intent(in) :: cpar
    character(len=*), intent(in) :: filename
    
    integer(i4b) :: unit, i

    if (cpar%myid /= 0) return

    unit = getlun()
    open(unit, file=trim(filename), recl=1024)
    write(unit,*) '# Band          Gain         delta_bp'
    do i = 1, numband
       write(unit,'(a15,f10.4,f10.2)') adjustl(trim(data(i)%label)), data(i)%gain, data(i)%bp%delta
    end do
    close(unit)

  end subroutine output_inst_params

end module comm_output_mod
