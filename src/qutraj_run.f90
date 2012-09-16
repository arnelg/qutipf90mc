module qutraj_run

  use qutraj_solver

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

end module
