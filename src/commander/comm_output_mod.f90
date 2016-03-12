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
    real(dp)                     :: chisq
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
    call output_inst_params(cpar, trim(cpar%outdir)//'/instpar_'//trim(postfix)//'.dat')    

    ! Output component results
    c => compList
    do while (associated(c))
       call c%dumpFITS(postfix, cpar%outdir)
       select type (c)
       class is (comm_diffuse_comp)
          call c%Cl%writeFITS(cpar%mychain, iter)
       end select
       c => c%next()
    end do

    ! Output channel-specific residual maps
    do i = 1, numband
       map => compute_residual(i)
       call map%writeFITS(trim(cpar%outdir)//'/res_'//trim(data(i)%label)//'_'// &
            & trim(postfix)//'.fits')
       call map%dealloc()
    end do

    ! Compute and output chi-square
    info => comm_mapinfo(cpar%comm_chain, cpar%nside_chisq, 0, cpar%nmaps_chisq, cpar%pol_chisq)
    map  => comm_map(info)
    call compute_chisq(chisq_map=map, chisq_fullsky=chisq)
    call map%writeFITS(trim(cpar%outdir)//'/chisq_'// trim(postfix) //'.fits')
    call map%dealloc()
    if (cpar%myid == cpar%root) write(*,fmt='(a,i4,a,e16.3)') &
         & 'Chain = ', cpar%mychain, ' -- chisq = ', chisq


    
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
