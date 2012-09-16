from numpy import *
from scipy import *
from qutip import *
import f90mcsolve as mc

psi0 = basis(2,0)
psi0d = psi0.data.toarray()
a = sigmax()

mc.f90mcsolve.init_psi0(psi0d,size(psi0d))
mc.f90mcsolve.init_hamiltonian(a.data.data,a.data.indices,a.data.indptr[0:2],a.data.nnz,a.data.shape[0],a.data.shape[1])

mc.f90mcsolve.evolve()
