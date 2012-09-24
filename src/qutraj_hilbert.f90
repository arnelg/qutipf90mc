module qutraj_hilbert

  use qutraj_precision
  use qutraj_general

  implicit none

  !
  ! Types
  !

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

  !
  ! work variables
  !

  complex(wp), allocatable, target :: work(:)


  !
  ! Interfaces
  !

  interface new
    module procedure state_init
    module procedure state_init2
    module procedure operat_init
    module procedure operat_init2
    module procedure operat_list_init
  end interface

  interface finalize
    module procedure state_finalize
    module procedure operat_finalize
    module procedure operat_list_finalize
  end interface

  interface operator(*)
    !module procedure state_state_mult
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
    integer :: istat=0
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("state_init: could not allocate.",istat)
    endif
  end subroutine
  subroutine state_init2(this,val)
    !type(state), intent(out) :: this
    complex(wp), allocatable :: this(:)
    complex(sp), intent(in) :: val(:)
    call state_init(this,size(val))
    this = val
  end subroutine

  subroutine state_finalize(this)
    !type(state), intent(inout) :: this
    complex(wp), allocatable :: this(:)
    integer :: istat=0
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    if (istat.ne.0) then
      call error("state_finalize: could not deallocate.",istat)
    endif
  end subroutine

  subroutine operat_init(this,nnz,nptr)
    ! todo: add special support for Hermitian matrix
    type(operat), intent(out) :: this
    integer, intent(in) :: nnz,nptr
    integer :: istat=0
    this%nnz = nnz
    if (allocated(this%a)) then
      deallocate(this%a,stat=istat)
    endif
    if (allocated(this%ia1)) then
      deallocate(this%ia1,stat=istat)
    endif
    if (allocated(this%pb)) then
      deallocate(this%pb,stat=istat)
    endif
    if (allocated(this%pe)) then
      deallocate(this%pe,stat=istat)
    endif
    allocate(this%a(nnz),stat=istat)
    if (istat.ne.0) then
      call fatal_error("operat_init: could not allocate.",istat)
    endif
    allocate(this%ia1(nnz),stat=istat)
    if (istat.ne.0) then
      call fatal_error("operat_init: could not allocate.",istat)
    endif
    allocate(this%pb(nptr),stat=istat)
    if (istat.ne.0) then
      call fatal_error("operat_init: could not allocate.",istat)
    endif
    allocate(this%pe(nptr),stat=istat)
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
  subroutine operat_init2(this,nnz,nptr,val,col,ptr,m,k)
    integer, intent(in) :: nnz,nptr,m,k
    type(operat), intent(out) :: this
    complex(sp), intent(in) :: val(nnz)
    integer, intent(in) :: col(nnz),ptr(nptr)
    integer :: i
    call operat_init(this,nnz,nptr)
    if (m.ne.k) then
      call fatal_error("operat_init2: # rows should equal # cols for operator type.")
    endif
    this%m = m
    this%k = k
    this%a = val
    this%ia1 = col
    this%pb = ptr
    do i=1,nnz-1
      this%pe(i) = this%pb(i+1)
    enddo
    this%pe(nnz) = nnz+1
  end subroutine

  subroutine operat_list_init(this,n)
    type(operat), intent(inout), allocatable :: this(:)
    integer, intent(in) :: n
    integer :: istat
    if (allocated(this)) then
      deallocate(this,stat=istat)
    endif
    allocate(this(n),stat=istat)
    if (istat.ne.0) then
      call fatal_error("operat_list_init: could not allocate.",istat)
    endif
  end subroutine

  subroutine operat_finalize(this)
    type(operat), intent(inout) :: this
    integer :: istat=0
    if (allocated(this%a)) then
      deallocate(this%a,stat=istat)
      if (istat.ne.0) then
        call error("operat_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%ia1)) then
      deallocate(this%ia1,stat=istat)
      if (istat.ne.0) then
        call error("operat_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%pb)) then
      deallocate(this%pb,stat=istat)
      if (istat.ne.0) then
        call error("operat_finalize: could not deallocate.",istat)
      endif
    endif
    if (allocated(this%pe)) then
      deallocate(this%pe,stat=istat)
      if (istat.ne.0) then
        call error("operat_finalize: could not deallocate.",istat)
      endif
    endif
  end subroutine

  subroutine operat_list_finalize(this)
    type(operat), intent(inout), allocatable :: this(:)
    integer :: istat=0,i
    if (allocated(this)) then
      do i=1,size(this)
        call finalize(this(i))
      enddo
      deallocate(this,stat=istat)
    endif
    if (istat.ne.0) then
      call error("operat_list_finalize: could not deallocate.",istat)
    endif
  end subroutine

  !
  ! State/operator arithmetic
  !

  function braket(fi,psi)
    ! return <fi|psi>
    complex(wp) :: braket
    complex(wp), intent(in) :: fi(:),psi(:)
    braket = sum(conjg(fi)*psi)
  end function

  subroutine normalize(psi)
    complex(wp), intent(inout) :: psi(:)
    real(wp) :: tmp
    tmp = sqrt(real(braket(psi,psi)))
    ! Check for division by zero
    if (abs(tmp) < epsi) then
      psi = 0.
    else
      psi = psi/tmp
    end if
  end subroutine


  function operat_state_mult(oper,psi)
    complex(wp), pointer :: operat_state_mult(:)
    type(operat), intent(in) :: oper
    complex(wp), intent(in) :: psi(:)
    !complex(wp), allocatable :: work
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
