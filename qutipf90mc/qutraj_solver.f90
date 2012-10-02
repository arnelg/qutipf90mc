!
! Deprecated
!

module qutraj_solver

  use qutraj_precision
  use qutraj_general
  use qutraj_hilbert

  implicit none

  ! Defines the RHS, to be sent to zvode
  !external rhs
  !external dummy_jac

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
    ! work array iwork should have length 30 for non-stiff
    integer :: liw = 0
    integer, allocatable :: iwork(:)
    ! method flag mf should be 10 for non-stiff
    integer :: mf = 10
    ! arbitrary real/complex and int array for user def input to rhs
    double complex :: rpar(1)
    integer :: ipar(1)
    ! abs. tolerance, rel. tolerance 
    double precision, allocatable :: atol(:), rtol(:)
    ! iopt=number of optional inputs, itol=1 for atol scalar, 2 otherwise
    integer :: iopt, itol
  end type

  !
  ! (Public) Data defining the problem
  !

  type(operat) :: hamilt
  type(operat), allocatable :: c_ops(:), e_ops(:)
  type(odeoptions) :: ode


  !
  ! Interfaces
  !

  interface new
    module procedure int_array_init
    module procedure int_array_init2
    module procedure r_array_init
    module procedure sp_array_init2
    module procedure dp_array_init2
  end interface

  interface finalize
    module procedure r_array_finalize
    module procedure odeoptions_finalize
  end interface

  !
  ! Subs and funcs
  !

  contains

  !
  ! Initializers and finalizers
  !

  subroutine int_array_init(this,n)
    integer, allocatable, intent(inout) :: this(:)
    integer, intent(in) :: n
    integer :: istat
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("int_array_init: could not allocate.",istat)
    endif
  end subroutine
  subroutine int_array_init2(this,val)
    integer, allocatable, intent(inout) :: this(:)
    integer, intent(in), dimension(:) :: val
    call int_array_init(this,size(val))
    this = val
  end subroutine

  subroutine r_array_init(this,n)
    real(wp), allocatable, intent(inout) :: this(:)
    integer, intent(in) :: n
    integer :: istat
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("sp_array_init: could not allocate.",istat)
    endif
  end subroutine
  subroutine sp_array_init2(this,val)
    real(wp), allocatable, intent(inout) :: this(:)
    real(sp), intent(in), dimension(:) :: val
    call r_array_init(this,size(val))
    this = val
  end subroutine
  subroutine dp_array_init2(this,val)
    real(wp), allocatable, intent(inout) :: this(:)
    real(dp), intent(in), dimension(:) :: val
    call r_array_init(this,size(val))
    this = val
  end subroutine

  subroutine r_array_finalize(this)
    real(wp), allocatable, intent(inout) :: this(:)
    integer :: istat=0
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    if (istat.ne.0) then
      call error("sp_array_finalize: could not deallocate.",istat)
    endif
  end subroutine

  subroutine odeoptions_finalize(this)
    type(odeoptions), intent(inout) :: this
    integer :: istat
    if (allocated(this%zwork)) then
      deallocate(this%zwork,stat=istat)
      if (istat.ne.0) then
        call error("odeoptions_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%rwork)) then
      deallocate(this%rwork,stat=istat)
      if (istat.ne.0) then
        call error("odeoptions_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%iwork)) then
      deallocate(this%iwork,stat=istat)
      if (istat.ne.0) then
        call error("odeoptions_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%atol)) then
      deallocate(this%atol,stat=istat)
      if (istat.ne.0) then
        call error("odeoptions_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%rtol)) then
      deallocate(this%rtol,stat=istat)
      if (istat.ne.0) then
        call error("odeoptions_finalize: could not deallocate.",istat)
      endif
    endif
  end subroutine

  !
  ! Evolution subs
  !

  subroutine nojump(y,t,tout,itask,istate,ode)
    ! evolve with effective hamiltonian
    type(odeoptions), intent(in) :: ode
    double complex, intent(inout) :: y(:)
    double precision, intent(inout) :: t
    double precision, intent(in) :: tout
    integer, intent(in) :: itask
    integer, intent(inout) :: istate
    integer :: istat

    call zvode(rhs,ode%neq,y,t,tout,ode%itol,ode%rtol,ode%atol,&
      itask,istate,ode%iopt,ode%zwork,ode%lzw,ode%rwork,ode%lrw,&
      ode%iwork,ode%liw,dummy_jac,ode%mf,ode%rpar,ode%ipar)
  end subroutine

  !
  ! RHS for zvode
  !

  subroutine rhs (neq, t, y, ydot, rpar, ipar)
    ! evolve with effective hamiltonian
    complex(wp) :: y(neq), ydot(neq),rpar
    real(wp) :: t
    integer :: ipar,neq
    !integer :: pntr,i
    ydot = -ii*(hamilt*y)
    !ydot = 0.
    !do i = 1, hamilt%M
    !  pntr = hamilt%pb(i)
    !  do while(pntr.lt.hamilt%pe(i))
    !    ydot(i) = ydot(i) &
    !        + hamilt%A(pntr) * y(hamilt%IA1(pntr))
    !     pntr = pntr + 1
    !  end do
    !end do
    !ydot = -(0.,1.)*ydot
  end subroutine

  subroutine dummy_jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
    complex(wp) :: y(neq), pd(nrpd,neq), rpar
    real(wp) :: t
    integer :: neq,ml,mu,nrpd,ipar
    return
  end


end module


