module qutraj_general

  implicit none

  !
  ! Global constants and general purpose subroutines
  !

  !
  ! Constants
  !

  integer, parameter :: sp = kind(1.0e0) ! single precision
  integer, parameter :: dp = kind(1.0d0) ! double precision
  integer, parameter :: wp = dp ! working precision
  integer, parameter :: wpc = wp ! working precision complex numbers
  integer, parameter :: blas_error_param    = -23

  ! small and large number
  real, parameter :: epsi=5*epsilon(1.0)
  real, parameter :: huge1=0.2*huge(1.0)

  ! imaginary unit
  complex(wp), parameter :: ii = (0._wp,1._wp)

  contains

  subroutine error(errormsg,ierror)
    character(len=*), intent(in), optional :: errormsg
    integer, intent(in), optional :: ierror
    if (present(errormsg)) then
      write(*,*) 'error: ',errormsg
    endif
    if (present(ierror)) then
      write(*,*) 'error flag=',ierror
    endif
  end subroutine

  subroutine fatal_error(errormsg,ierror)
    character(len=*), intent(in), optional :: errormsg
    integer, intent(in), optional :: ierror
    if (present(errormsg)) then
      write(*,*) 'fatal error: ',errormsg
    endif
    if (present(ierror)) then
      write(*,*) 'error flag=',ierror
    endif
    write(*,*) 'halting'
    stop 1
  end subroutine

end module
