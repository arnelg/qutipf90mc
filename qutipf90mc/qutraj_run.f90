!
! TODO:
!
!
!

module qutraj_run
  !
  ! This is the main module on the fortran side of things
  !

  use qutraj_precision
  use qutraj_general
  use qutraj_hilbert
  use mt19937
  use linked_list

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
    ! task and state of solver
    integer :: itask, istate
  end type

  !
  ! Data defining the problem
  !
  ! Invisible to python: hamilt, c_ops, e_opts, ode
  ! (because f2py can't handle derived types)
  !

  type(operat) :: hamilt
  type(operat), allocatable :: c_ops(:), e_ops(:)
  type(odeoptions) :: ode

  real(wp), allocatable :: tlist(:)
  complex(wp), allocatable :: psi0(:)

  integer :: ntraj=1
  integer :: norm_steps = 5
  real(wp) :: norm_tol = 0.001
  integer :: n_c_ops = 0
  integer :: n_e_ops = 0
  logical :: mc_avg = .true.

  ! Optional ode options, 0 means use default values
  integer :: order=0,nsteps=0
  double precision :: first_step=0,min_step=0,max_step=0

  ! Solution
  ! format:
  ! all states: sol(1,trajectory,time,y(:))
  ! avg. kets: sol(1,1,time,y(:))
  ! all expect: sol(e_ops(i),trajectory,time,expecation value)
  ! avg. expect: sol(e_ops(i),1,time,expectation value)
  ! if returning averaged dense density matrices:
  ! sol(1,time,rho_i,rho_j)
  complex(wp), allocatable :: sol(:,:,:,:)
  ! if returning averaged density matrices in sparse CSR format,
  ! use the following solution array and get_rho_sparse instead.
  type(operat), allocatable :: sol_rho(:)

  ! return density martices or averaged kets?
  logical :: return_kets = .true.
  ! return dense or sparse density matrices?
  logical :: rho_return_sparse = .true.

  ! temporary storage for csr matrix, available for python
  ! this is needed because you can't send assumed
  ! shape arrays to python
  complex(wp), allocatable :: csr_val(:)
  integer, allocatable :: csr_col(:), csr_ptr(:)
  integer :: csr_nrows,csr_ncols

  ! Collapse times and integer denoting which operator did it
  real(wp), allocatable :: col_times(:)
  integer, allocatable :: col_which(:)
  ! data temporarily stored in linked lists...
  type(linkedlist_real) :: ll_col_times
  type(linkedlist_int) :: ll_col_which
  type(llnode_real), pointer :: realnode
  type(llnode_int), pointer :: intnode

  ! Integer denoting the type of unravelling
  ! 1 for jump unravelling
  ! diffusive unravellings to be implemented
  integer :: unravel_type = 1

  !
  ! Interfaces
  !

  interface finalize
    module procedure odeoptions_finalize
  end interface

  contains

  !
  ! Initialize problem
  !

  subroutine init_tlist(val,n)
    use qutraj_precision
    real(wp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(tlist,val)
  end subroutine

  subroutine init_psi0(val,n)
    use qutraj_precision
    complex(wp), intent(in) :: val(n)
    integer, intent(in) :: n
    call new(psi0,val)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,m,k,nnz,nptr)
    use qutraj_precision
    integer, intent(in) :: nnz,nptr,m,k
    complex(wp), intent(in)  :: val(nnz)
    integer, intent(in) :: col(nnz),ptr(nptr)
    call new(hamilt,val,col,ptr,m,k)
  end subroutine

  subroutine init_c_ops(i,n,val,col,ptr,m,k,first,nnz,nptr)
    use qutraj_precision
    integer, intent(in) :: i,n
    integer, intent(in) :: nnz,nptr,m,k
    complex(wp), intent(in) :: val(nnz)
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
    call new(c_ops(i),val,col,ptr,m,k)
  end subroutine

  subroutine init_e_ops(i,n,val,col,ptr,m,k,first,nnz,nptr)
    use qutraj_precision
    integer, intent(in) :: i,n
    integer, intent(in) :: nnz,nptr,m,k
    complex(wp), intent(in) :: val(nnz)
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
    call new(e_ops(i),val,col,ptr,m,k)
  end subroutine

  subroutine init_odedata(neq,atol,rtol,mf,&
      lzw,lrw,liw,ml,mu,natol,nrtol)
    use qutraj_precision
    integer, intent(in) :: neq
    integer, intent(in), optional :: lzw,lrw,liw,mf
    integer, intent(in) :: natol,nrtol
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

  subroutine get_rho_sparse(i)
    integer, intent(in) :: i
    call new(csr_val,sol_rho(i)%a)
    call new(csr_col,sol_rho(i)%ia1)
    call new(csr_ptr,sol_rho(i)%pb)
    csr_nrows = sol_rho(i)%m
    csr_ncols = sol_rho(i)%k
  end subroutine

  !
  ! Evolution
  !

  subroutine evolve(states,instanceno,rngseed)
    ! Save states or expectation values?
    logical, intent(in) :: states
    ! What process # am I?
    integer, intent(in) :: instanceno,rngseed
    double precision :: t, tout
    double complex, allocatable :: y(:),y_tmp(:),rho(:,:)
    type(operat) :: rho_sparse
    integer :: istat=0,istat2=0,traj,progress
    integer :: i,j,l,m,n
    real(wp) :: mu,nu
    real(wp), allocatable :: p(:)

    ! Allocate solution array
    if (allocated(sol)) then
      deallocate(sol,stat=istat)
      if (istat.ne.0) then
        call error("evolve: could not deallocate.",istat)
      endif
    endif
    if (states) then
      if (mc_avg) then
        if (return_kets) then
          allocate(sol(1,1,size(tlist),ode%neq),stat=istat)
          sol = (0.,0.)
        else
          if (rho_return_sparse) then
            call new(sol_rho,size(tlist))
            call new(rho_sparse,1,1)
          else
            allocate(sol(1,size(tlist),ode%neq,ode%neq),stat=istat)
            allocate(rho(ode%neq,ode%neq),stat=istat2)
            sol = (0.,0.)
            rho = (0.,0.)
          endif
        endif
      else
        allocate(sol(1,ntraj,size(tlist),ode%neq),stat=istat)
        sol = (0.,0.)
      endif
    else 
      if (mc_avg) then
        allocate(sol(n_e_ops,1,size(tlist),1),stat=istat)
        sol = (0.,0.)
      else
        allocate(sol(n_e_ops,ntraj,size(tlist),1),stat=istat)
        sol = (0.,0.)
      endif
    endif
    if (istat.ne.0) call fatal_error("evolve: could not allocate solution.",&
      istat)
    if (istat2.ne.0) call fatal_error("evolve: could not allocate rho.",&
      istat2)

    ! Allocate work arrays
    call new(y,ode%neq)
    call new(y_tmp,ode%neq)
    ! Allocate tmp array for jump probabilities
    call new(p,n_c_ops)
    ! Initalize rng
    call init_genrand(rngseed)

    ! Loop over trajectories
    progress = 1
    do traj=1,ntraj
      ! two random numbers
      mu = grnd()
      nu = grnd()
      ! First call to zvode
      ode%istate = 1
      ! Initial values
      y = psi0
      ! Initial value of indep. variable
      t = tlist(1)
      do i=1,size(tlist)
        ! Solution wanted at
        if (i==1) then
          ! found this to be necessary due to round off error
          tout = t
        else
          tout = tlist(i)
        endif
        select case(unravel_type)
        case(1)
          call evolve_jump(t,tout,y,y_tmp,p,mu,nu)
        case default
          call fatal_error('Unknown unravel type.')
        end select
        y_tmp = y
        call normalize(y_tmp)
        ! Compute solution
        if (states) then
          if (mc_avg) then
            if (return_kets) then
              sol(1,1,i,:) = sol(1,1,i,:) + y_tmp
            else
              ! construct density matrix
              if (rho_return_sparse) then
                call densitymatrix_sparse(y_tmp,rho_sparse)
                if (traj==1) then
                  sol_rho(i) = rho_sparse
                else
                  sol_rho(i) = sol_rho(i) + rho_sparse
                endif
              else
                call densitymatrix_dense(y_tmp,rho)
                sol(1,i,:,:) = sol(1,i,:,:) + rho
              endif
            endif
          else
            sol(1,traj,i,:) = y_tmp
          endif
        else
          if (mc_avg) then
            do l=1,n_e_ops
              sol(l,1,i,1) = sol(l,1,i,1)+braket(y_tmp,e_ops(l)*y_tmp)
            enddo
          else
            do l=1,n_e_ops
              sol(l,traj,i,1) = braket(y_tmp,e_ops(l)*y_tmp)
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
    ! Normalize
    if (mc_avg) then
      if (states .and. .not.return_kets .and. rho_return_sparse) then
        do j=1,size(sol_rho)
          sol_rho(j) = (1._wp/ntraj)*sol_rho(j)
        enddo
      else
        sol = (1._wp/ntraj)*sol
      endif
    endif
    ! Turn linked lists into arrays
    call ll_to_array(ll_col_times,col_times)
    call ll_to_array(ll_col_which,col_which)
    ! Deallocate
    call finalize(y)
    call finalize(y_tmp)
    call finalize(p)
  end subroutine

  !
  ! Misc
  !

  ! Deallocate stuff

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
    call finalize(sol_rho)
    if (allocated(sol)) then
      deallocate(sol,stat=istat)
    endif
    if (istat.ne.0) then
      call error("finalize_sol: could not deallocate.",istat)
    endif
  end subroutine

  ! Misc

  subroutine test_real_precision
    use qutraj_precision
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


  !
  ! Stuff not visible to python
  !

  !
  ! Evolution subs
  !

  subroutine evolve_jump(t,tout,y,y_tmp,p,mu,nu)
    !
    ! Evolve quantum trajectory y(t) to y(tout) using ``jump'' method
    !
    ! Input: t, tout, y
    ! Work arrays: y_tmp, p
    ! mu, nu: two random numbers
    !
    double complex, intent(inout) :: y(:),y_tmp(:)
    double precision, intent(inout) :: t, tout
    real(wp), intent(inout) :: p(:)
    real(wp), intent(inout) :: mu,nu
    double precision :: t_prev, t_final, t_guess
    integer :: j,k
    integer :: cnt
    real(wp) :: norm2_psi,norm2_prev,norm2_guess,sump
    logical, save :: first = .true.
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
    !          indicate this (unless ITASK = 4 and TOUT comes before 
    !          TCRIT, in which case answers at T = TOUT are returned 
    !          first).
    if (first) then
      ! integrate one step at the time, w/o overshooting
      ode%itask = 5
      ! set optinal arguments
      ! see zvode.f
      ode%rwork = 0.0
      ode%iwork = 0
      ode%rwork(5) = first_step
      ode%rwork(6) = max_step
      ode%rwork(7) = min_step
      ode%iwork(5) = order
      ode%iwork(6) = nsteps
      ode%iopt = 1
      ! first call to zvode
      ode%istate = 1
      first = .false.
    endif

    ode%rwork(1) = tout
    norm2_psi = abs(braket(y,y))
    !write(*,*) y
    do while(t<tout)
      t_prev = t
      y_tmp = y
      norm2_prev = norm2_psi
      call nojump(y,t,tout,ode%itask,ode)
      if (ode%istate.lt.0) then
        write(*,*) "zvode error: istate=",ode%istate
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
          y = y_tmp
          t = t_prev
          call nojump(y,t,t_guess,1,ode)
          if (ode%istate.lt.0) then
            write(*,*) "zvode failed after adjusting step size. istate=",ode%istate
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
              y_tmp=y
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
            ! Append collapse time and operator # to linked lists
            call new(realnode,t)
            call new(intnode,j)
            call append(ll_col_times,realnode)
            call append(ll_col_which,intnode)
          endif
          sump = sump+p(j)
        enddo
        ! new random numbers
        mu = grnd()
        nu = grnd()
        ! normalize y
        call normalize(y)
        ! reset, first call to zvode
        ode%istate = 1
      endif
    enddo
  end subroutine

  subroutine nojump(y,t,tout,itask,ode)
    ! evolve with effective hamiltonian
    type(odeoptions), intent(in) :: ode
    double complex, intent(inout) :: y(:)
    double precision, intent(inout) :: t
    double precision, intent(in) :: tout
    integer, intent(in) :: itask
    integer :: istat

    call zvode(rhs,ode%neq,y,t,tout,ode%itol,ode%rtol,ode%atol,&
      itask,ode%istate,ode%iopt,ode%zwork,ode%lzw,ode%rwork,ode%lrw,&
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
    ydot = -ii*(hamilt*y)
  end subroutine

  subroutine dummy_jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
    ! dummy jacobian for zvode
    complex(wp) :: y(neq), pd(nrpd,neq), rpar
    real(wp) :: t
    integer :: neq,ml,mu,nrpd,ipar
    return
  end

end module
