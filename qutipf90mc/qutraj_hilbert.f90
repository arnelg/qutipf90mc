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

    ! m = number of rows, k = number of cols
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

  !complex(wp), allocatable, target :: work(:)


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
    do i=1,nptr-1
      this%pe(i) = this%pb(i+1)
    enddo
    this%pe(nptr) = nnz+1
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
    !braket = 1.
  end function

  subroutine normalize(psi)
    complex(wp), intent(inout) :: psi(:)
    real(wp) :: tmp
    tmp = sqrt(abs(braket(psi,psi)))
    ! Check for division by zero
    if (abs(tmp) < epsi) then
      psi = 0.
    else
      psi = psi/tmp
    end if
  end subroutine

  function operat_state_mult(oper,psi)
    complex(wp), intent(in) :: psi(:)
    type(operat), intent(in) :: oper
    complex(wp):: operat_state_mult(size(psi))
    complex(wp), allocatable :: tmp(:)
    integer :: ierr
    !integer, allocatable :: ost(:)
    call new(tmp,size(psi))
    call sparse_mv_mult(oper,psi,tmp,ierr)
    if (ierr.ne.0) then
      call error("operate_state_mult: error",ierr)
    endif
    !allocate(ost(size(oper%pb)+1))
    !ost = oper%pb
    !ost(size(ost)) = oper%nnz+1
    !call amux(oper%k,psi,tmp,oper%a,oper%ia1,ost)
    operat_state_mult = tmp
    call finalize(tmp)
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
    !if (diag.eq.'U') then !process unstored diagonal
    !   if (m.eq.n) then
    !      y = x
    !   else
    !      ierr = blas_error_param
    !      return
    !   end if
    !end if
    !if ((type.eq.'S').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
    !   if (part.eq.'U') then
    !      do i = 1, mat%M
    !         pntr = mat%pb(i)
    !         do while(pntr.lt.mat%pe(i))
    !            if(i.eq.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !            else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !               y(mat%IA1(pntr + ofs) + ofs) =  &
    !        y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
    !            end if
    !            pntr = pntr + 1
    !         end do
    !      end do
    !   else
    !      do i = 1, mat%M
    !         pntr = mat%pb(i)
    !         do while(pntr.lt.mat%pe(i))
    !            if(i.eq.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !            else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !               y(mat%IA1(pntr + ofs) + ofs) = &
    !        y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
    !            end if
    !            pntr = pntr + 1
    !         end do
    !      end do
    !   end if
    !   ierr = 0
    !else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
    !   if (part.eq.'U') then
    !      do i = 1, mat%M
    !         pntr = mat%pb(i)
    !         do while(pntr.lt.mat%pe(i))
    !            if(i.eq.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !            else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !              y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
    !                     + conjg (mat%A(pntr + ofs)) * x(i) 
    !            end if
    !            pntr = pntr + 1
    !         end do
    !      end do
    !   else
    !      do i = 1, mat%M
    !         pntr = mat%pb(i)
    !         do while(pntr.lt.mat%pe(i))
    !            if(i.eq.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !            else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
    !               y(i) = y(i) &
    !            + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
    !             y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
    !                     + conjg (mat%A(pntr + ofs)) * x(i) 
    !            end if
    !            pntr = pntr + 1
    !         end do
    !      end do
    !   end if
    !   ierr = 0
    !else
       do i = 1, mat%M
          pntr = mat%pb(i)
          do while(pntr.lt.mat%pe(i))
            y(i) = y(i) &
                + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
             pntr = pntr + 1
          end do
       end do
       ierr = 0
    !end if
  end subroutine


  subroutine amux ( n, x, y, a, ja, ia )
  ! Aadapted from sparsekit

  !*****************************************************************************80
  !
  !! AMUX multiplies a CSR matrix A times a vector.
  !
  !  Discussion:
  !
  !    This routine multiplies a matrix by a vector using the dot product form.
  !    Matrix A is stored in compressed sparse row storage.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer N, the row dimension of the matrix.
  !
  !    Input, real X(*), and array of length equal to the column dimension 
  !    of A.
  !
  !    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  !    Output, real Y(N), the product A * X.
  !
    implicit none

    integer n

    complex ( kind = wp ) a(*)
    integer i
    integer ia(*)
    integer ja(*)
    integer k
    complex ( kind = wp ) t
    complex ( kind = wp ) x(*)
    complex ( kind = wp ) y(n)

    do i = 1, n
  !
  !  Compute the inner product of row I with vector X.
  !
      t = (0.0,0.0)
      do k = ia(i), ia(i+1)-1
        t = t + a(k) * x(ja(k))
      end do

      y(i) = t

    end do

    return
  end

subroutine amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )
  ! Aadapted from sparsekit

!*****************************************************************************80
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer JOB, job indicator.  When JOB = 0, only the structure
!    is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow
  integer nzmax

  real ( kind = wp ) a(*)
  real ( kind = wp ) b(*)
  real ( kind = wp ) c(nzmax)
  integer ia(nrow+1)
  integer ib(ncol+1)
  integer ic(ncol+1)
  integer ierr
  integer ii
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(nzmax)
  integer jcol
  integer jj
  integer job
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  real ( kind = wp ) scal
  logical values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              if ( values ) then
                c(jpos) = c(jpos) + scal * b(kb)
              end if
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
end
end module
