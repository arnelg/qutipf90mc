program zvodetest

  ! Defines the RHS
  external F
  ! No. of ODES
  integer,parameter :: neq=1
  ! Initial values
  double complex :: y(neq)
  ! Initial value of indep. variable
  double precision :: t = 0


end program

subroutine F (NEQ, T, Y, YDOT, RPAR, IPAR)
  DOUBLE COMPLEX Y(NEQ), YDOT(NEQ)
  DOUBLE PRECISION T
end subroutine

