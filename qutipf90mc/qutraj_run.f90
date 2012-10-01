!
! TODO:
!
! - return array of density matrices when mc_avg=.true.
!
module qutraj_run
  use qutraj_precision
  use qutraj_general
  use qutraj_solver
  use qutraj_hilbert
  use mt19937

  implicit none

  !
  ! (Public) Data defining the problem
  !
  ! Some data is hidden in qutraj_solver instead.
  !

  real(wp), allocatable :: tlist(:)
  complex(wp), allocatable :: psi0(:)
  integer :: ntraj=1
  integer :: norm_steps = 5
  real(wp) :: norm_tol = 0.001
  integer :: n_c_ops = 0
  integer :: n_e_ops = 0
  logical :: mc_avg = .true.

  ! Ode options, 0 means use default values
  integer :: order=0,nsteps=0
  double precision :: first_step=0,min_step=0,max_step=0

  ! Solution
  ! format:
  ! all states: sol(1,ntraj,size(tlist),1,y(t))
  ! avg. density mat: sol(1,1,size(tlist),rho(t))
  ! all expect: sol(n_e_ops,ntraj,size(tlist),1,1)
  ! avg. expect: sol(n_e_ops,1,size(tlist),1,1)
  complex(wp), allocatable :: sol(:,:,:,:,:)

  contains

  !
  ! Initialize problem
  !

  subroutine init_tlist(val,n)
    real(sp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(tlist,val)
  end subroutine

  subroutine init_psi0(val,n)
    complex(sp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(psi0,val)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,m,k,nnz,nptr)
    integer, intent(in) :: nnz,nptr,m,k
    complex(sp), intent(in)  :: val(nnz)
    integer, intent(in) :: col(nnz),ptr(nptr)
    call new(hamilt,nnz,nptr,val,col+1,ptr+1,m,k)
  end subroutine

  subroutine init_c_ops(i,n,val,col,ptr,m,k,nnz,nptr,first)
    integer, intent(in) :: i,n
    integer, intent(in) :: nnz,nptr,m,k
    complex(sp), intent(in) :: val(nnz)
    integer, intent(in) :: col(nnz),ptr(nptr)
    logical, optional :: first
    if (.not.present(first)) then
      first = .false.
    endif
    if (first) then
      call new(c_ops,n)
    endif
    if (.not.allocated(c_ops)) then
      call error('init_c_ops: c_ops not allocated. call with first=True first.')
    endif
    n_c_ops = n
    call new(c_ops(i+1),nnz,nptr,val,col+1,ptr+1,m,k)
  end subroutine

  subroutine init_e_ops(i,n,val,col,ptr,m,k,nnz,nptr,first)
    integer, intent(in) :: i,n
    integer, intent(in) :: nnz,nptr,m,k
    complex(sp), intent(in) :: val(nnz)
    integer, intent(in) :: col(nnz),ptr(nptr)
    logical, optional :: first
    if (.not.present(first)) then
      first = .false.
    endif
    if (first) then
      call new(e_ops,n)
    endif
    if (.not.allocated(e_ops)) then
      call error('init_e_ops: e_ops not allocated. call with first=True first.')
    endif
    n_e_ops = n
    call new(e_ops(i+1),nnz,nptr,val,col+1,ptr+1,m,k)
  end subroutine

  subroutine init_odedata(neq,atol,rtol,mf,&
      lzw,lrw,liw,ml,mu,natol,nrtol)
    integer, intent(in) :: neq
    integer, intent(in), optional :: lzw,lrw,liw,mf
    integer, intent(in) :: natol,nrtol
    !real(sp), intent(in) :: atol(natol), rtol(nrtol)
    double precision, optional :: atol(1),rtol(1)
    integer, intent(in), optional :: ml,mu
    integer :: istat

    ode%neq = neq
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

    call new(ode%zwork,ode%lzw)
    call new(ode%rwork,ode%lrw)
    call new(ode%iwork,ode%liw)
    call new(ode%atol,atol)
    call new(ode%rtol,rtol)
    if (size(ode%atol)==1) then
      ode%itol=1
    else
      ode%itol=2
    endif

    ode%iopt = 0
  end subroutine

  !
  ! Evolution
  !

  subroutine evolve(states,instanceno)
    ! Save states or expectation values?
    logical, intent(in) :: states
    integer, intent(in) :: instanceno
    double precision :: t, tout, t_prev, t_final, t_guess
    double complex, allocatable :: y(:),y_prev(:),y_tmp(:),rho(:,:)
    integer :: istate,itask
    integer :: istat=0,i,j,k,traj,progress
    integer :: l,m,n,cnt
    real(wp) :: nu,mu,norm2_psi,norm2_prev,norm2_guess,sump
    real(wp), allocatable :: p(:)
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
    if (allocated(sol)) then
      deallocate(sol,stat=istat)
      if (istat.ne.0) then
        call error("evolve: could not deallocate.",istat)
      endif
    endif
    if (states) then
      if (mc_avg) then
        allocate(sol(1,1,size(tlist),ode%neq,ode%neq),stat=istat)
        if (istat.ne.0) call fatal_error("evolve: could not allocate.")
        allocate(rho(ode%neq,ode%neq),stat=istat)
        if (istat.ne.0) call fatal_error("evolve: could not allocate.")
        rho = (0.,0.)
      else
        allocate(sol(1,ntraj,size(tlist),1,ode%neq),stat=istat)
        if (istat.ne.0) call fatal_error("evolve: could not allocate.")
      endif
    else 
      if (mc_avg) then
        allocate(sol(n_e_ops,1,size(tlist),1,1),stat=istat)
        if (istat.ne.0) call fatal_error("evolve: could not allocate.")
      else
        allocate(sol(n_e_ops,ntraj,size(tlist),1,1),stat=istat)
        if (istat.ne.0) call fatal_error("evolve: could not allocate.")
      endif
    endif
    sol = (0.,0.)
    ! Allocate solution
    call new(y,ode%neq)
    call new(y_prev,ode%neq)
    call new(y_tmp,ode%neq)
    ! Allocate tmp array for jump probabilities
    call new(p,n_c_ops)

    ! integrate one step at the time, w/o overshooting
    itask = 5

    ! set optinal arguments
    ode%rwork = 0.0
    ode%iwork = 0
    ode%rwork(5) = first_step
    ode%rwork(6) = max_step
    ode%rwork(7) = min_step
    ode%iwork(5) = order
    ode%iwork(6) = nsteps
    ode%iopt = 1

    ! Loop over trajectories
    progress = 1
    ! Initalize rng
    call init_genrand(instanceno)
    do traj=1,ntraj
      ! two random numbers
      mu = grnd()
      nu = grnd()

      ! first call to zvode
      istate = 1
      ! Initial values
      y = psi0
      ! Initial value of indep. variable
      t = tlist(1)
      do i=1,size(tlist)
        ! Solution wanted at
        if (i==1) then
          tout = t
        else
          tout = tlist(i)
        endif
        ode%rwork(1) = tout
        norm2_psi = abs(braket(y,y))
        do while(t<tout)
          t_prev = t
          y_prev = y
          norm2_prev = norm2_psi
          call nojump(y,t,tout,itask,istate,ode)
          if (istate.lt.0) then
            write(*,*) "zvode error: istate=",istate
            !stop
          endif
          ! prob of nojump
          norm2_psi = abs(braket(y,y))
          if (norm2_psi.le.mu) then
            ! jump happened
            ! find collapse time to specified tolerance
            t_final = t
            cnt=1
            do k=1,norm_steps
              !t_guess=t_prev+(mu-norm2_prev)&
              !  /(norm2_psi-norm2_prev)*(t_final-t_prev)
              t_guess=t_prev+log(norm2_prev/mu)&
                /log(norm2_prev/norm2_psi)*(t_final-t_prev)
              if (t_guess<t_prev .or. t_guess>t_final) then
                t_guess = t_prev+0.5*(t_final-t_prev)
              endif
              y = y_prev
              t = t_prev
              call nojump(y,t,t_guess,1,istate,ode)
              if (istate.lt.0) then
                write(*,*) "zvode failed after adjusting step size. istate=",istate
                !stop
              endif
              norm2_guess = abs(braket(y,y))
              if (abs(mu-norm2_guess) < norm_tol*mu) then
                  exit
              elseif (norm2_guess < mu) then
                  ! t_guess is still > t_jump
                  t_final=t_guess
                  norm2_psi=norm2_guess
              else
                  ! t_guess < t_jump
                  t_prev=t_guess
                  y_prev=y
                  norm2_prev=norm2_guess
              endif
              cnt = cnt+1
            enddo
            if (cnt > norm_steps) then
              call error("Norm tolerance not reached. Increase accuracy of ODE solver or norm_steps.")
            endif
            ! determine which jump
            do j=1,n_c_ops
              y_tmp = c_ops(j)*y
              p(j) = abs(braket(y_tmp,y_tmp))
            enddo
            p = p/sum(p)
            sump = 0
            do j=1,n_c_ops
              if ((sump <= nu) .and. (nu < sump+p(j))) then
                y = c_ops(j)*y
              endif
              sump = sump+p(j)
            enddo
            ! new random numbers
            mu = grnd()
            nu = grnd()
            ! normalize y
            call normalize(y)
            ! reset, first call to zvode
            istate = 1
          endif
        enddo
        y_tmp = y
        call normalize(y_tmp)
        if (states) then
          if (mc_avg) then
            ! construct density matrix
            call densitymatrix(y_tmp,rho)
            sol(1,1,i,:,:) = sol(1,1,i,:,:) + rho
          else
            sol(1,traj,i,1,:) = y_tmp
          endif
        else
          if (mc_avg) then
            do l=1,n_e_ops
              sol(l,1,i,1,1) = sol(l,1,i,1,1)+braket(y_tmp,e_ops(l)*y_tmp)
            enddo
          else
            do l=1,n_e_ops
              sol(l,traj,i,1,1) = braket(y_tmp,e_ops(l)*y_tmp)
            enddo
          endif
        endif
        ! End time loop
      enddo
      ! Indicate progress
      if (instanceno == 1 .and. traj.ge.progress*ntraj/10.0) then
        write(*,*) "progress of process 1: ", progress*10, "%"
        progress=progress+1
      endif
      ! End loop over trajectories
    enddo
    ! normalize
    if (mc_avg) then
      sol = sol/ntraj
    endif
    ! Deallocate
    call finalize(y)
    call finalize(y_prev)
    call finalize(y_tmp)
    call finalize(p)
  end subroutine

  !
  ! Misc
  !

  ! Deallocate stuff

  subroutine finalize_work
    integer :: istat=0
    call finalize(psi0)
    call finalize(hamilt)
    call finalize(c_ops)
    call finalize(e_ops)
    call finalize(ode)
  end subroutine

  subroutine finalize_sol
    integer :: istat=0
    call finalize(tlist)
    if (allocated(sol)) then
      deallocate(sol,stat=istat)
    endif
    if (istat.ne.0) then
      call error("finalize_all: could not deallocate.",istat)
    endif
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
        write(*,*) "number of decimals working precision: ",i-1
      endif
      i = i+1
    enddo
  end subroutine

end module
