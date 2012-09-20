from numpy import *
from scipy import *
from qutip import *
import qutraj_run as qt

def init_psi0(psi0):
    psi0d = psi0.data.toarray()
    qt.qutraj_run.init_psi0(psi0d,size(psi0d))

def init_hamilt(H):
    qt.qutraj_run.init_hamiltonian(H.data.data,H.data.indices,H.data.indptr[0:size(H.data.indptr)-1],H.data.nnz,H.data.shape[0],H.data.shape[1])

neq = 2
psi0 = basis(neq,0)
H = sigmaz()

# Times
T = 1.0
dt = 0.1
nstep = T/dt
tlist = linspace(0,T,nstep)

qt.qutraj_run.init_tlist(tlist)

init_psi0(psi0)
init_hamilt(H)
qt.qutraj_run.init_odedata(neq,1e-5,1e-5,mf=10)

qt.qutraj_run.evolve()

qt.qutraj_run.finalize_all()

#

sol = mcsolve(H,psi0,tlist,[],[],ntraj=1)

