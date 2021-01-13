!! Fortran 90 interface routines for FFTW3
module planck_fft
!use planck_config
  use healpix_types
implicit none
private

public real_fft, complex_fft, make_plan, destroy_plan, planck_fft_plan

logical, parameter, public :: fft_forward=.false., fft_backward=.true.

!! Data type containing information about the requested FFT (length, direction
!! etc.). Do not manipulate its contents directly; always use make_plan() and
!! destroy_plan() instead.
type planck_fft_plan
  logical :: direction=fft_forward
  integer(i4b) :: length=-1
  integer(i8b) :: fftw_plan=0
end type

interface complex_fft
  module procedure s_c_complex_fft, s_r_complex_fft, &
                   d_c_complex_fft, d_r_complex_fft
end interface

interface real_fft
  module procedure s_real_fft, d_real_fft
end interface

include "fftw3.f"

contains

subroutine sanity_check (plan, len)
  type(planck_fft_plan), intent(in) :: plan
  integer, intent(in) :: len

  if (len/=plan%length) &
    call exit_with_status(1,"planck_fft: invalid plan for this transform")
end subroutine

!! Create a plan for a FFT with the given length and direction.
!! Do not forget to call destroy_plan() on your plan after use, because
!! this will lead to memory leaks.
subroutine make_plan (plan,length,direction)
  type(planck_fft_plan), intent(out) :: plan
  integer, intent(in) :: length
  logical, intent(in) :: direction
  complex(dp) arr_in(length), arr_out(length)

  integer dir

  plan%length=length
  plan%direction=direction

  dir = merge (FFTW_BACKWARD, FFTW_FORWARD, direction)
!$OMP CRITICAL
  call dfftw_plan_dft_1d (plan%fftw_plan, length, arr_in, arr_out, &
                          dir, FFTW_ESTIMATE+FFTW_UNALIGNED+FFTW_DESTROY_INPUT)
!$OMP END CRITICAL
end subroutine

!! Free the resources associated with a plan via make_plan().
!! Do not use your plan in an FFT after destroying it!
!! It is safe to call destroy_plan() on an uninitialized plan.
subroutine destroy_plan (plan)
  type(planck_fft_plan), intent(inout) :: plan

!$OMP CRITICAL
  if (plan%length>0) call dfftw_destroy_plan (plan%fftw_plan)
!$OMP END CRITICAL
  plan%fftw_plan=0
  plan%direction=fft_forward
  plan%length=-1
end subroutine

subroutine d_c_complex_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  complex(dp), intent(inout) :: data(:)
  real(dp) data2(2*size(data)), data3(2*size(data))

  call sanity_check (plan, size(data))
  data2 = transfer (data, data2)
  call dfftw_execute_dft (plan%fftw_plan, data2, data3)
  data = transfer (data3, data)
end subroutine

subroutine d_r_complex_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)
  real(dp) data2(2*size(data))

  call sanity_check (plan, size(data)/2)
  call dfftw_execute_dft (plan%fftw_plan, data, data2)
  data = transfer (data2, data)
end subroutine

subroutine s_c_complex_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  complex(sp), intent(inout) :: data(:)
  real(dp) data2(2*size(data)), data3(2*size(data))
  integer i

  call sanity_check (plan, size(data))
  do i=1,plan%length
    data2(2*i-1) = real (data(i))
    data2(2*i  ) = aimag(data(i))
  end do
  call dfftw_execute_dft (plan%fftw_plan, data2, data3)
  do i=1,plan%length
    data(i) = cmplx(data3(2*i-1),data3(2*i))
  end do
end subroutine

subroutine s_r_complex_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  real(sp), intent(inout) :: data(:)
  real(dp) data2(size(data)), data3(size(data))

  call sanity_check (plan, size(data)/2)
  data2=data
  call dfftw_execute_dft (plan%fftw_plan, data2, data3)
  data=data3
end subroutine

subroutine d_real_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  real(dp), intent(inout) :: data(:)

  real(dp) data2(2*size(data)), data3(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft_forward) then
    data2 (2:2*n:2) = 0
    data2 (1:2*n-1:2) = data
    call dfftw_execute_dft (plan%fftw_plan, data2, data3)
    data(1) = data3(1)
    data(2:n) = data3(3:n+1)
  else
    data2(1) = data(1)
    data2(2) = 0
    data2(3:n+1) = data(2:n)
    data2(n+2:) = 0
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call dfftw_execute_dft (plan%fftw_plan, data2, data3)
    data = data3(1:2*n-1:2)
  endif
end subroutine

subroutine s_real_fft (plan, data)
  type(planck_fft_plan), intent(in) :: plan
  real(sp), intent(inout) :: data(:)

  real(dp) data2(2*size(data)), data3(2*size(data))
  integer n, i

  call sanity_check (plan, size(data))
  n = size(data)

  if (plan%direction .eqv. fft_forward) then
    data2 (2:2*n:2) = 0
    data2 (1:2*n-1:2) = data
    call dfftw_execute_dft (plan%fftw_plan, data2, data3)
    data(1) = data3(1)
    data(2:n) = data3(3:n+1)
  else
    data2(1) = data(1)
    data2(2) = 0
    data2(3:n+1) = data(2:n)
    data2(n+2:) = 0
    do i=1, n/2
      data2(2*n-2*i+1)=  data2(2*i+1)
      data2(2*n-2*i+2)= -data2(2*i+2)
    end do
    call dfftw_execute_dft (plan%fftw_plan, data2, data3)
    data = data3(1:2*n-1:2)
  endif
end subroutine

end module planck_fft
