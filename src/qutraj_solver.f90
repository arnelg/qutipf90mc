module qutraj_solver

  use qutraj_general
  use qutraj_hilbert

  implicit none

  ! Defines the RHS, to be sent to zvode
  external rhs
  external dummy_jac

  !
  ! Types
  !

  type odeoptions
    ! No. of ODES
    integer :: neq=1
    ! work array zwork should have length 15*neq for non-stiff
    integer :: lzw = 0
    double complex, allocatable :: zwork(:)
    ! work array rwork should have length 20+neq for non-siff
    integer :: lrw = 0
    double precision, allocatable :: rwork(:)
    ! work array iwrok should have length 30 for non-stiff
    integer :: liw = 0
    integer, allocatable :: iwork(:)
    ! method flag mf should be 10 for non-stiff
    integer :: mf = 10
    ! arbitrary real/complex and int array for user def input to rhs
    double complex :: rpar(1)
    integer :: ipar(1)
    ! abs. tolerance, rel. tolerance 
    double precision, allocatable :: atol(:), rtol(:)
    ! itask=1 for normal output
    ! iopt=number of optional inputs, itol=1 for atol scalar, 2 otherwise
    integer :: itask, iopt, itol
  end type

  !
  ! (Public) Data defining the problem
  !

  real(sp), allocatable :: tlist(:)
  type(state) :: psi0,psi
  type(operat) :: hamilt
  type(odeoptions) :: ode

  !
  ! Interfaces
  !

  interface new
    module procedure sp_array_init
    module procedure sp_array_init2
  end interface

  interface finalize
    module procedure sp_array_finalize
    module procedure odeoptions_finalize
  end interface

  !
  ! Subs and funcs
  !

  contains

  !
  ! Initializers and finalizers
  !

  subroutine sp_array_init(this,n)
    real(sp), allocatable, intent(inout) :: this(:)
    integer, intent(in) :: n
    integer :: istat
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("sp_array_init: could not allocate.",istat)
    endif
  end subroutine
  subroutine sp_array_init2(this,val)
    real(sp), allocatable, intent(inout) :: this(:)
    real(sp), intent(in), dimension(:) :: val
    call sp_array_init(this,size(val))
    this = val
  end subroutine

  subroutine sp_array_finalize(this)
    real(sp), allocatable, intent(inout) :: this
    integer :: istat
    deallocate(this,stat=istat)
    if (istat.ne.0) then
      call error("sp_array_finalize: could not deallocate.",istat)
    endif
  end subroutine

  subroutine odeoptions_finalize(this)
    type(odeoptions), intent(inout) :: this
    integer :: istat
    deallocate(this%zwork,this%rwork,this%iwork,this%atol,this%rtol,stat=istat)
    if (istat.ne.0) then
      call error("odeoptions_finalize: could not deallocate.",istat)
    endif
  end subroutine

  !
  ! Evolution subs
  !

  subroutine nojump(y,t,tout,istate)
    double complex, intent(inout) :: y(:)
    double precision, intent(inout) :: t
    double precision, intent(in) :: tout
    integer, intent(inout) :: istate
    integer :: istat

    call zvode(rhs,ode%neq,y,t,tout,ode%itol,ode%rtol,ode%atol,&
      ode%itask,istate,ode%iopt,ode%zwork,ode%lzw,ode%rwork,ode%lrw,&
      ode%iwork,ode%liw,dummy_jac,ode%mf,ode%rpar,ode%ipar)
  end subroutine

end module

!
! RHS for zvode
!

subroutine rhs (neq, t, y, ydot, rpar, ipar)
  use qutraj_hilbert
  use qutraj_solver
  double complex y(neq), ydot(neq),rpar
  double precision t
  integer ipar,neq
  !type(state) :: dpsi
  !ydot(1) = y(1)
  psi%x = y
  psi = (-ii)*(hamilt*psi)
  !write(*,*) psi%x
  ydot = psi%x
end subroutine

subroutine dummy_jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
  double complex y(neq), pd(nrpd,neq), rpar
  double precision t
  integer neq,ml,mu,nrpd,ipar
  return
end



