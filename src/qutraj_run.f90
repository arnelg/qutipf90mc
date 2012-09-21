module qutraj_run

  use qutraj_precision
  use qutraj_general
  use qutraj_solver

  implicit none

  !
  ! (Public) Data defining the problem
  !
  ! Some data is hidden in qutraj_solver instead.
  !

  double precision, allocatable :: tlist(:)
  complex(wp), allocatable :: psi0(:)
  integer :: ntraj=1
  integer :: norm_steps = 5
  real(wp) :: norm_tol = 0.001

  ! Solution
  complex(wp), allocatable :: sol(:,:)

  contains

  ! Initialize problem

  subroutine init_tlist(val,n)
    real(sp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(tlist,val)
  end subroutine

  subroutine init_psi0(val,n)
    complex(sp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(psi0,val)
    call new(work,n)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,nnz,nrows,ncols)
    integer, intent(in) :: nnz,nrows,ncols
    complex(sp), intent(in), dimension(nnz) :: val
    integer, intent(in), dimension(nnz) :: col,ptr
    call new(hamilt,nnz,val,col+1,ptr+1,nrows,ncols)
  end subroutine

  subroutine init_c_ops(i,n,val,col,ptr,nnz,nrows,ncols)
    integer, intent(in) :: i,n
    integer, intent(in) :: nnz,nrows,ncols
    complex(sp), intent(in), dimension(nnz) :: val
    integer, intent(in), dimension(nnz) :: col,ptr
    if (.not.allocated(c_ops)) then
      call new(c_ops,n)
    endif
    call new(c_ops(i),nnz,val,col+1,ptr+1,nrows,ncols)
    write(*,*) c_ops(i)%a
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
    call finalize(work)
    call finalize(hamilt)
    call finalize(c_ops)
    call finalize(ode)
  end subroutine

  ! Evolution

  subroutine evolve
    double precision :: t, tout
    double complex, allocatable :: y(:)
    integer :: istate,itask
    integer :: istat,i
    real(wp) :: nu
    ! ITASK  = An index specifying the task to be performed.
    !          Input only.  ITASK has the following values and meanings.
    !          1  means normal computation of output values of y(t) at
    !             t = TOUT (by overshooting and interpolating).
    !          2  means take one step only and return.
    !          3  means stop at the first internal mesh point at or
    !             beyond t = TOUT and return.
    !          4  means normal computation of output values of y(t) at
    !             t = TOUT but without overshooting t = TCRIT.
    !             TCRIT must be input as RWORK(1).  TCRIT may be equal to
    !             or beyond TOUT, but not behind it in the direction of
    !             integration.  This option is useful if the problem
    !             has a singularity at or beyond t = TCRIT.
    !          5  means take one step, without passing TCRIT, and return.
    !             TCRIT must be input as RWORK(1).
    !
    !          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
    !          (within roundoff), it will return T = TCRIT (exactly) to
    !          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
    !          in which case answers at T = TOUT are returned first).

    ! Allocate solution array
    allocate(sol(size(tlist),ode%neq),stat=istat)
    ! Allocate solution
    allocate(y(ode%neq),stat=istat)
    if (istat.ne.0) then
      call fatal_error("evolve: could not allocate.",istat)
    endif

    ! integrate one step at the time, w/o overshooting
    itask = 5
    ! first call to zvode
    istate = 1
    ! Initial values
    y = psi0

    !write(*,*) ode%neq,ode%itol,ode%rtol,ode%atol
    !write(*,*) itask,istate,ode%iopt
    !write(*,*) ode%lzw,ode%lrw,ode%liw
    !write(*,*) ode%mf

    ! Initalize rng
    call init_genrand(1)
    nu = grnd()

    ! Initial value of indep. variable
    t = tlist(i)
    do i=1,size(tlist)
      ! Solution wanted at
      if (i==1) then
        tout = t
      else
        tout = tlist(i)
        ode%rwork(1) = tout
      endif

      do while(t<tout)
        call nojump(y,t,tout,itask,istate)
        write(*,*) t, y
      enddo
      sol(i,:) = y
      if (istate.lt.0) then
        write(*,*) "error: istate=",istate
        stop
      endif
    enddo
  end subroutine

  ! Misc
  subroutine test_real_precision
    real(wp) :: b,a
    integer :: i
    write(*,*) "wp=",wp
    b = 1.0
    a = 1.0
    i = 1
    do while (b.ne.b+a)
      a = a*0.1
      if (b==b+a) then
        write(*,*) "number of decimals precision: ",i-1
      endif
      i = i+1
    enddo
  end subroutine

end module
