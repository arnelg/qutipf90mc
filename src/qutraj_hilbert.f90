module qutraj_hilbert

  use qutraj_general

  implicit none

  !
  ! Types
  !

  !type state
  !  ! States are dense one dimensional vectors
  !  integer :: n
  !  complex(wp), allocatable :: x(:)
  !end type

  type operat
    ! Operators are represented as spare matrices
    ! stored in compressed row format (CSR)

    ! m = number of cols, k = number of rows
    integer :: m,k
    ! number of values
    integer :: nnz
    ! compression format is CSR
    character*5 :: fida = 'CSR'
    ! base: Fortran or C base
    integer :: base = 0 ! comp with python
    ! diag: 'U' for un-stored diag entries, assumed to be one
    character*11 :: diag = 'N'
    ! typem: 'S' for symmetric, 'H' for Hermitian
    character*11 :: typem = 'G'
    !both/lower/upper half of matrix specified
    character*11 :: part = 'B'
    ! values
    complex(wp), allocatable :: a(:)
    integer, allocatable :: ia1(:),pb(:),pe(:)
    !complex(wp), pointer :: a(:)
    !integer, pointer :: ia1(:),pb(:),pe(:)
    ! notice: pe(i) = pb(i+1)-1
  end type

  ! work variables
  !type(state) :: work
  complex(wp), allocatable, target :: work(:)


  !
  ! Interfaces
  !

  interface new
    module procedure state_init
    module procedure state_init2
    module procedure operat_init
    module procedure operat_init2
  end interface

  interface finalize
    module procedure state_finalize
    module procedure operat_finalize
  end interface

  !interface operator(=)
  !  module procedure state_state_eq
  !end interface

  !interface operator(+)
  !  module procedure state_state_add
  !end interface

  !interface operator(-)
  !  module procedure state_state_sub
  !end interface

  interface operator(*)
    !module procedure real_state_mult
    !module procedure complex_state_mult
    module procedure operat_state_mult
  end interface

  !
  ! Subs and funcs
  !

  contains

  !
  ! Initializers & finalizers
  !

  subroutine state_init(this,n)
    !type(state), intent(out) :: this
    complex(wp), allocatable :: this(:)
    integer, intent(in) :: n
    integer :: istat
    !this%n = n
    !allocate(this%x(n),stat=istat)
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("state_init: could not allocate.",istat)
    endif
  end subroutine
  subroutine state_init2(this,val)
    !type(state), intent(out) :: this
    complex(wp), allocatable :: this(:)
    complex, intent(in), dimension(:) :: val
    call state_init(this,size(val))
    this = val
  end subroutine

  subroutine state_finalize(this)
    !type(state), intent(inout) :: this
    complex(wp), allocatable :: this(:)
    integer :: istat
    deallocate(this,stat=istat)
    if (istat.ne.0) then
      call error("state_finalize: could not deallocate.",istat)
    endif
  end subroutine

  subroutine operat_init(this,nnz)
    ! todo: add special support for Hermitian matrix
    type(operat), intent(out) :: this
    integer, intent(in) :: nnz
    integer :: istat
    this%nnz = nnz
    allocate(this%a(nnz),stat=istat)
    allocate(this%ia1(nnz),stat=istat)
    allocate(this%pb(nnz),stat=istat)
    allocate(this%pe(nnz),stat=istat)
    if (istat.ne.0) then
      call fatal_error("operat_init: could not allocate.",istat)
    endif
    ! Set default parameters
    this%fida = 'CSR'
    this%base = 1 ! fortran base
    this%diag = 'N'
    this%typem = 'G'
    this%part = 'B'
  end subroutine

  subroutine operat_init2(this,nnz,val,col,ptr,nrows,ncols)
    type(operat), intent(out) :: this
    complex, intent(in) :: val(:)
    integer, intent(in) :: col(:),ptr(:)
    integer, intent(in) :: nnz,nrows,ncols
    integer :: i
    call operat_init(this,nnz)
    if (nrows.ne.ncols) then
      call fatal_error("operat_init2: nrows should equal ncols &
        for operator type.")
    endif
    this%m = ncols
    this%k = nrows
    this%a = val
    this%ia1 = col
    this%pb = ptr
    do i=1,nnz-1
      this%pe(i) = this%pb(i+1)
    enddo
    this%pe(nnz) = nnz+1
  end subroutine

  subroutine operat_finalize(this)
    type(operat), intent(inout) :: this
    integer :: istat
    deallocate(this%a,this%ia1,this%pb,this%pe,stat=istat)
    !nullify(this%a,this%ia1,this%pb,this%pe)
    if (istat.ne.0) then
      call error("operat_finalize: could not deallocate.",istat)
    endif
  end subroutine

  !
  ! State arithmetic
  !

  !type(state) function state_state_eq(fi,psi)
  !  type(state), intent(in) :: fi,psi
  !  work%x = psi%x
  !  state_state_eq = work
  !end function

  !type(state) function state_state_add(fi,psi)
  !  ! |fi> + |psi>
  !  type(state), intent(in) :: fi
  !  type(state), intent(in) :: psi
  !  work%x = fi%x+psi%x
  !  state_state_add = work
  !end function

  !type(state) function state_state_sub(fi,psi)
  !  ! |fi> + |psi>
  !  type(state), intent(in) :: fi
  !  type(state), intent(in) :: psi
  !  work%x = fi%x-psi%x
  !  state_state_sub = work
  !end function

  !type(state) function real_state_mult(c,psi)
  !  type(state), intent(in) :: psi
  !  real(wp), intent(in) :: c
  !  work%x = c*psi%x
  !  real_state_mult = work
  !end function

  !type(state) function complex_state_mult(c,psi)
  !  type(state), intent(in) :: psi
  !  complex(wp), intent(in) :: c
  !  work%x = c*psi%x
  !  complex_state_mult = work
  !end function

  !
  ! Matrix vector multiplicatoin
  !

  function operat_state_mult(oper,psi)
    complex(wp), pointer :: operat_state_mult(:)
    type(operat), intent(in) :: oper
    !type(state), intent(in) :: psi
    complex(wp), allocatable, intent(in) :: psi(:)
    !complex(wp), allocatable :: work
    !type(state) :: work
    integer :: ierr

    !if (psi%n.ne.work%n) then
    if (size(psi).ne.size(work)) then
      write(*,*) "operate_state_mult: state has wrong size:",size(psi)
      write(*,*) "should be:",size(work)
      write(*,*) "have you properly initialized 'work' state?"
      call fatal_error
      return
    endif
    !call sparse_mv_mult(oper,psi%x,work%x,ierr)
    call sparse_mv_mult(oper,psi,work,ierr)
    if (ierr.ne.0) then
      call error("operate_state_mult: error",ierr)
    endif
    operat_state_mult => work
  end function

  subroutine sparse_mv_mult(mat,x,y,ierr)
    ! y = Ax
    ! Adapted from sparse blas
    type(operat) :: mat
    complex(KIND=wp) , dimension(:), intent(in) :: x
    complex(KIND=wp) , dimension(:), intent(out) :: y
    integer, intent(out) :: ierr
    integer :: m,n,base,ofs,i,pntr
    character :: diag,type,part
    ierr = -1
    m = size(y)
    n = size(x)
    if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
       ierr = blas_error_param
       return
    end if
    !call get_infoa(mat%INFOA,'b',base,ierr)
    base = mat%base
    ofs = 1 - base
    !call get_descra(mat%DESCRA,'d',diag,ierr)
    diag = mat%diag
    !call get_descra(mat%DESCRA,'t',type,ierr)
    type = mat%typem
    !call get_descra(mat%DESCRA,'a',part,ierr)
    part = mat%part
    y = (0.0d0, 0.0d0) 
    if (diag.eq.'U') then !process unstored diagonal
       if (m.eq.n) then
          y = x
       else
          ierr = blas_error_param
          return
       end if
    end if
    if ((type.eq.'S').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
       if (part.eq.'U') then
          do i = 1, mat%M
             pntr = mat%pb(i)
             do while(pntr.lt.mat%pe(i))
                if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr + ofs) + ofs) =  &
            y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                end if
                pntr = pntr + 1
             end do
          end do
       else
          do i = 1, mat%M
             pntr = mat%pb(i)
             do while(pntr.lt.mat%pe(i))
                if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr + ofs) + ofs) = &
            y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                end if
                pntr = pntr + 1
             end do
          end do
       end if
       ierr = 0
    else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
       if (part.eq.'U') then
          do i = 1, mat%M
             pntr = mat%pb(i)
             do while(pntr.lt.mat%pe(i))
                if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                         + conjg (mat%A(pntr + ofs)) * x(i) 
                end if
                pntr = pntr + 1
             end do
          end do
       else
          do i = 1, mat%M
             pntr = mat%pb(i)
             do while(pntr.lt.mat%pe(i))
                if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                   y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                 y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                         + conjg (mat%A(pntr + ofs)) * x(i) 
                end if
                pntr = pntr + 1
             end do
          end do
       end if
       ierr = 0
    else
       do i = 1, mat%M
          pntr = mat%pb(i)
          do while(pntr.lt.mat%pe(i))
            y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
             pntr = pntr + 1
          end do
       end do
       ierr = 0
    end if
  end subroutine

end module
