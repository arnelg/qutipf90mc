module f90mcsolve

  use solver

  implicit none

  contains

  subroutine hello(x,y)
    real, intent(in) :: x
    real, intent(out) :: y
    write(*,*) 'hello from fortran'
    y = x**2
    ! try to use Bessel function from scipy
    !b = (0.0,1.0)
    !write(*,*) jn(0,b)
  end subroutine

  subroutine init_psi0(val,n)
    complex, intent(in), dimension(n) :: val
    integer, intent(in) :: n
    write(*,*) val
    call new(psi0,n,val)
    write(*,*) psi0%x
    call new(work,n)
  end subroutine

  subroutine init_hamiltonian(val,col,ptr,nnz,nrows,ncols)
    integer, intent(in) :: nnz,nrows,ncols
    complex, intent(in), dimension(nnz) :: val
    integer, intent(in), dimension(nnz) :: col,ptr
    write(*,*) val
    call new(hamilt,nnz,val,col,ptr,nrows,ncols)
    write(*,*) hamilt%a
    write(*,*) hamilt%ia1
    write(*,*) hamilt%pb
    write(*,*) hamilt%pe
  end subroutine

  subroutine evolve
    type(state) :: psi
    psi = hamilt*psi0
    write(*,*) psi%x
  end subroutine




!subroutine amux (n, x, y, a,ja,ia) 
!  real*8, intent(in) ::  x(:), a(:)
!  real*8, intent(out) ::  y(:)
!  integer, intent(in) :: n, ja(:), ia(:)
!!-----------------------------------------------------------------------
!!         A times a vector
!!----------------------------------------------------------------------- 
!! multiplies a matrix by a vector using the dot product form
!! Matrix A is stored in compressed sparse row storage.
!!
!! on entry:
!!----------
!! n     = row dimension of A
!! x     = real array of length equal to the column dimension of
!!         the A matrix.
!! a, ja,
!!    ia = input matrix in compressed sparse row format.
!!
!! on return:
!!-----------
!! y     = real array of length n, containing the product y=Ax
!!
!!-----------------------------------------------------------------------
!! local variables
!!
!  real*8 t
!  integer i, k
!!-----------------------------------------------------------------------
!  do i = 1,n
!!
!!     compute the inner product of row i with vector x
!!
!     t = 0.0d0
!     do k=ia(i), ia(i+1)-1 
!        t = t + a(k)*x(ja(k))
!     enddo
!!
!!     store result in y(i) 
!!
!     y(i) = t
!  enddo
!!c---------end-of-amux---------------------------------------------------
!!c-----------------------------------------------------------------------
!end subroutine

end module
