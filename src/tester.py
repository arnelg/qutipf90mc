from numpy import *
from scipy import *
from qutip import *
import qutraj_run as qt

psi0 = basis(2,0)
psi0d = psi0.data.toarray()
a = sigmax()

qt.qutraj_run.init_psi0(psi0d,size(psi0d))
qt.qutraj_run.init_hamiltonian(a.data.data,a.data.indices,a.data.indptr[0:2],a.data.nnz,a.data.shape[0],a.data.shape[1])

qt.qutraj_run.evolve()
