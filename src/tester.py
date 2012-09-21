from numpy import *
from scipy import *
from qutip import *
from pylab import *
import qutraj_run as qt

# Working precision
wpr = dtype(float32)
wpc = dtype(complex64)

def init_tlist(tlist):
    qt.qutraj_run.init_tlist(array(tlist,dtype=wpr))
    #qt.qutraj_run.init_tlist(tlist)

def init_psi0(psi0):
    psi0d = array(psi0.data.toarray(),dtype=wpc).transpose()
    qt.qutraj_run.init_psi0(psi0d,size(psi0d))

def init_hamilt(H):
    datad = array(H.data.data,dtype=wpc)
    qt.qutraj_run.init_hamiltonian(datad,H.data.indices,H.data.indptr[0:size(H.data.indptr)-1],H.data.nnz,H.data.shape[0],H.data.shape[1])

def init_c_ops(c_ops):
    n = len(c_ops)
    for i in range(n):
        datad = array(c_ops[i].data.data,dtype=wpc)
        qt.qutraj_run.init_c_ops(i+1,n,datad,c_ops[i].data.indices,c_ops[i].data.indptr[0:size(c_ops[i].data.indptr)-1],c_ops[i].data.nnz,c_ops[i].data.shape[0],c_ops[i].data.shape[1])

def get_states():
    states=array([Qobj()]*nstep)
    for i in range(len(tlist)):
        states[i] = Qobj(matrix(qt.qutraj_run.sol[i]).transpose())
    return states

def finalize():
    qt.qutraj_run.finalize_all()

neq = 2
psi0 = basis(neq,0)
H = sigmax()
gamma = 0.1
c_ops = [gamma*sigmax()]

# Times
T = 10.0
dt = 0.1
nstep = int(T/dt)
tlist = linspace(0,T,nstep)

init_tlist(tlist)
init_psi0(psi0)
init_hamilt(H)
#init_c_ops(c_ops)

opts = Odeoptions()
atol = opts.atol
rtol = opts.rtol
qt.qutraj_run.init_odedata(neq,atol,rtol,mf=10)

qt.qutraj_run.evolve()
#sol = mcsolve(H,psi0,tlist,c_ops,[],ntraj=1)
sol = mcsolve(H,psi0,tlist,[],[],ntraj=1)

states = get_states()
states2 = sol.states

figure()
plot(sol.times,real(expect(sigmaz(),states)))
plot(sol.times,real(expect(sigmaz(),states2)))

finalize()
