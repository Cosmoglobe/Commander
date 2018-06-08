program cambtest
  use comm_camb_mod
  implicit none

  integer(i4b) :: l
  type(camb), pointer :: theta
  real(dp),   allocatable, dimension(:,:) :: cls

  theta => camb(3000, .false.)
  call theta%setParam(h=0.72d0, tau=0.05d0)
  call theta%getCls(cls)
  
  open(58,file='cl.dat',recl=1024)
  write(58,*) '#    l          D_l_TT          D_l_TE          D_l_EE          D_l_BB'
  do l = 2, size(cls,1)-1
     write(58,fmt='(i8,4f16.5)') l, cls(l,CAMB_TT), cls(l,CAMB_TE), cls(l,CAMB_EE), cls(l,CAMB_BB)
  end do
  close(58)

  deallocate(cls)

end program cambtest
