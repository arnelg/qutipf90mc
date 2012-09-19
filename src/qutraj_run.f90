module qutraj_run

  use qutraj_solver

  implicit none

  ! Defines the RHS, to be sent to zvode
  external rhs
  external dummy_jac

  contains

  ! Initialize problem

  subroutine init_psi0(val,n)
    complex, intent(in), dimension(n) :: val
    integer, intent(in) :: n
    call new(psi0,n,val)
    call new(work,n)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,nnz,nrows,ncols)
    integer, intent(in) :: nnz,nrows,ncols
    complex, intent(in), dimension(nnz) :: val
    integer, intent(in), dimension(nnz) :: col,ptr
    call new(hamilt,nnz,val,col,ptr,nrows,ncols)
  end subroutine

  subroutine init_odedata(neq,atol,rtol,mf,lzw,lrw,liw,ml,mu,natol,nrtol)
    integer, intent(in) :: neq
    integer, intent(in), optional :: lzw,lrw,liw,mf
    integer, intent(in) :: natol,nrtol
    double precision, intent(in) :: atol(natol), rtol(nrtol)
    !double precision, optional :: atol,rtol
    integer, intent(in), optional :: ml,mu
    integer :: istat

    ode%neq = neq
    !allocate(ode%y(neq),stat=istat)
    if (lzw.ne.0) then
      ode%lzw = lzw
    endif
    if (lrw.ne.0) then
      ode%lrw = lrw
    endif
    if (liw.ne.0) then
      ode%liw = liw
    endif
    if (lrw.eq.0) then
      ode%lrw = 20+neq
    endif

    if (mf==0 .or. mf==10) then
      ! assuming non-stiff by default
      ode%mf=10
      if (lzw.eq.0) then
        ode%lzw = 15*neq
      endif
      if (liw.eq.0) then
        ode%liw = 30
      endif
    elseif (mf==21.or.mf==22) then
      ode%mf = mf
      if (lzw.eq.0) then
        ode%lzw = 8*neq+2*neq**2
      endif
      if (liw.eq.0) then
        ode%liw = 30+neq
      endif
    elseif (mf==24.or.mf==25) then
      ode%mf = mf
      if (lzw.eq.0) then
        ! mf=24,25 requires ml and mu
        ode%lzw = 10*neq + (3*ml + 2*mu)*neq
      endif
      if (liw.eq.0) then
        ode%liw = 30+neq
      endif
    endif

    allocate(ode%zwork(ode%lzw),stat=istat)
    allocate(ode%rwork(ode%lrw),stat=istat)
    allocate(ode%iwork(ode%liw),stat=istat)

    allocate(ode%atol(size(atol)),stat=istat)
    ode%atol = atol
    allocate(ode%rtol(size(rtol)),stat=istat)
    ode%rtol = rtol
    if (size(ode%atol)==1) then
      ode%itol=1
    else
      ode%itol=2
    endif

    ode%itask = 1
    ode%iopt = 0

    if (istat.ne.0) then
      call fatal_error("init_odedata: could not allocate.",istat)
    endif
  end subroutine

  ! Deallocate everything

  subroutine finalize_all
    !deallocate(ode%zwork,ode%rwork,ode%iwork,ode%atol,ode%rtol)
    call finalize(psi0)
    call finalize(work)
    call finalize(hamilt)
    call finalize(ode)
  end subroutine

  ! Evolution

  subroutine evolve
    type(state) :: psi
    double precision :: t, tout
    double complex, allocatable :: y(:)
    integer :: istate
    integer :: istat

    allocate(y(ode%neq),stat=istat)
    if (istat.ne.0) then
      call error("evolve: could not allocate.",istat)
    endif

    ! Initial values
    y = psi0%x
    ! Initial value of indep. variable
    t = 0.
    ! Solution wanted at
    tout = 1.
    ! first call to zvode
    istate = 1

    !write(*,*) ode%neq,ode%itol,ode%rtol,ode%atol
    !write(*,*) ode%itask,istate,ode%iopt
    !write(*,*) ode%lzw,ode%lrw,ode%liw
    !write(*,*) ode%mf

    call zvode(rhs,ode%neq,y,t,tout,ode%itol,ode%rtol,ode%atol,&
      ode%itask,istate,ode%iopt,ode%zwork,ode%lzw,ode%rwork,ode%lrw,&
      ode%iwork,ode%liw,dummy_jac,ode%mf,ode%rpar,ode%ipar)

    write(*,*) t, y(1)
    if (istate.lt.0) then
      write(*,*) "error: istate=",istate
      stop
    endif
  end subroutine

end module

subroutine rhs (neq, t, y, ydot, rpar, ipar)
  double complex y(neq), ydot(neq),rpar
  double precision t
  integer ipar
  type(state) :: dpsi
  !ydot(1) = sin(t)*(1.,0.)
  !ydot(1) = y(1)
  dpsi = -ii*(hamilt*psi)
  ydot = dpsi%x
end subroutine

subroutine dummy_jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
  double complex y(neq), pd(nrpd,neq), rpar
  double precision t
  return
end

