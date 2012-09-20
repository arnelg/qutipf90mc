from numpy import *
from scipy import *
from qutip import *
from pylab import *
import qutraj_run as qt

def init_psi0(psi0):
    psi0d = psi0.data.toarray()
    qt.qutraj_run.init_psi0(psi0d,size(psi0d))

def init_hamilt(H):
    qt.qutraj_run.init_hamiltonian(H.data.data,H.data.indices,H.data.indptr[0:size(H.data.indptr)-1],H.data.nnz,H.data.shape[0],H.data.shape[1])

def get_states():
    states=array([Qobj()]*nstep)
    for i in range(len(tlist)):
        states[i] = Qobj(matrix(qt.qutraj_run.sol[i]).transpose())
    return states

def finalize():
    qt.qutraj_run.finalize_all()

neq = 2
psi0 = basis(neq,0)
H = sigmay()

# Times
T = 10.0
dt = 0.1
nstep = int(T/dt)
tlist = linspace(0,T,nstep)

qt.qutraj_run.init_tlist(tlist)

init_psi0(psi0)
init_hamilt(H)
atol = odeconfig.options.atol
rtol = odeconfig.options.rtol
qt.qutraj_run.init_odedata(neq,atol,rtol,mf=10)

qt.qutraj_run.evolve()


#qt.qutraj_run.finalize_all()

#

sol = mcsolve(H,psi0,tlist,[],[],ntraj=1)
#sol = mcsolve(H,psi0,tlist,[sigmax()],[],ntraj=2)

states = get_states()
states2 = sol.states

figure()
plot(sol.times,real(expect(sigmaz(),states)))
plot(sol.times,real(expect(sigmaz(),states2)))

finalize()
