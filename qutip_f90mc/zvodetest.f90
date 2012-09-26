! linpack: zgefa zgesl zgbfa zgbsl 
! blas: zcopy
program zvodetest

  implicit none

  ! Defines the RHS
  external rhs
  external dummy_jac
  ! No. of ODES
  integer,parameter :: neq=1
  ! work array zwork should have length 15*neq for non-stiff
  integer, parameter :: lzw = 15
  double complex, dimension(lzw) :: zwork
  ! work array rwork should have length 20+neq for non-siff
  integer, parameter :: lrw = 21
  double precision, dimension(lrw) :: rwork
  ! work array iwrok should have length 30 for non-stiff
  integer, parameter :: liw = 30
  integer, dimension(liw) :: iwork
  ! method flag mf should be 10 for non-stiff
  integer, parameter :: mf = 10
  ! arbitrary real/complex and int array for user def input to rhs
  double complex :: rpar(1)
  !complex :: rpar(1)
  integer :: ipar(1)

  double complex :: y(neq)
  double precision :: t, tout, atol(1), rtol(1)
  integer :: itask, istate, iopt, itol

  ! Initial values
  y(1) = (1.,0.)
  ! Initial value of indep. variable
  t = 0.
  ! Solution wanted at
  tout = 1.
  ! atol is a scalar
  itol = 1
  ! absolute tolerance
  atol(1) = 1e-5
  ! relative tolerance
  rtol(1) = 1e-5
  ! normal output
  itask = 1
  ! first call to zvode
  istate = 1
  ! no optional input
  iopt = 0


  call zvode(rhs,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt, &
    zwork,lzw,rwork,lrw,iwork,liw,dummy_jac,mf,rpar,ipar)

  write(*,*) t, y(1)
  if (istate.lt.0) then
    write(*,*) "error: istate=",istate
    stop
  endif

end program

subroutine rhs (neq, t, y, ydot, rpar, ipar)
  double complex y(neq), ydot(neq),rpar
  double precision t
  integer ipar
  !ydot(1) = sin(t)*(1.,0.)
  ydot(1) = y(1)
end subroutine

subroutine dummy_jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
  double complex y(neq), pd(nrpd,neq), rpar
  double precision t
  return
end

