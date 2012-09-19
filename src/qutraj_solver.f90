module qutraj_solver

  use qutraj_general
  use qutraj_hilbert

  implicit none

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

  type(state) :: psi0
  type(operat) :: hamilt
  type(odeoptions) :: ode

  !
  ! Interfaces
  !

  interface finalize
    module procedure odeoptions_finalize
  end interface

  !
  ! Subs and funcs
  !

  contains

  !
  ! Initializers and finalizers
  !

  subroutine odeoptions_finalize(this)
    type(odeoptions), intent(inout) :: this
    integer :: istat
    deallocate(this%zwork,this%rwork,this%iwork,this%atol,this%rtol,stat=istat)
    if (istat.ne.0) then
      call error("odeoptions_finalize: could not deallocate.",istat)
    endif
  end subroutine

end module
