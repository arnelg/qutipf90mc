module qutraj_run

  use qutraj_solver

  implicit none

  ! Solution
  complex, allocatable :: sol(:,:)

  contains

  ! Initialize problem

  subroutine init_tlist(val,n)
    real(sp), intent(in), dimension(n) :: val
    integer, intent(in) :: n
    call new(tlist,val)
  end subroutine

  subroutine init_psi0(val,n)
    complex, intent(in), dimension(n) :: val
    integer, intent(in) :: n
    call new(psi0,val)
    call new(psi,n)
    call new(work,n)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,nnz,nrows,ncols)
    integer, intent(in) :: nnz,nrows,ncols
    complex, intent(in), dimension(nnz) :: val
    integer, intent(in), dimension(nnz) :: col,ptr
    call new(hamilt,nnz,val,col+1,ptr+1,nrows,ncols)
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
    integer :: istat
    !deallocate(ode%zwork,ode%rwork,ode%iwork,ode%atol,ode%rtol)
    deallocate(tlist,sol,stat=istat)
    if (istat.ne.0) then
      call error("finalize_all: could not deallocate.",istat)
    endif
    !call finalize(tlist)
    call finalize(psi0)
    call finalize(psi)
    call finalize(work)
    call finalize(hamilt)
    call finalize(ode)
  end subroutine

  ! Evolution

  subroutine evolve
    double precision :: t, tout
    double complex, allocatable :: y(:)
    integer :: istate
    integer :: istat,i

    ! Allocate solution array
    allocate(sol(size(tlist),ode%neq),stat=istat)
    ! Allocate solution
    allocate(y(ode%neq),stat=istat)
    if (istat.ne.0) then
      call fatal_error("evolve: could not allocate.",istat)
    endif

    !write(*,*) ode%neq,ode%itol,ode%rtol,ode%atol
    !write(*,*) ode%itask,istate,ode%iopt
    !write(*,*) ode%lzw,ode%lrw,ode%liw
    !write(*,*) ode%mf

    ! first call to zvode
    istate = 1
    ! Initial values
    y = psi0

    ! Initial value of indep. variable
    t = tlist(i)
    do i=1,size(tlist)
      ! Solution wanted at
      tout = tlist(i)

      call nojump(y,t,tout,istate)
      !write(*,*) t, y
      sol(i,:) = y
      if (istate.lt.0) then
        write(*,*) "error: istate=",istate
        stop
      endif
    enddo
  end subroutine

end module
