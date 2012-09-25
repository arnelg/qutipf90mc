import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import qutraj_run as qt
import time

# Working precision
wpr = dtype(float64)
wpc = dtype(complex64)

def init_tlist(tlist):
    qt.qutraj_run.init_tlist(array(tlist,dtype=wpr))
    #qt.qutraj_run.init_tlist(tlist)

def init_psi0(psi0):
    psi0d = array(psi0.data.toarray(),dtype=wpc).transpose()
    qt.qutraj_run.init_psi0(psi0d,size(psi0d))

def init_hamilt(H,c_ops):
    # construct effective non-Hermitian Hamiltonian
    H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
        for i in range(len(c_ops))])
    #H_eff = H
    datad = array(H_eff.data.data,dtype=wpc)
    qt.qutraj_run.init_hamiltonian(datad,H_eff.data.indices,
            H_eff.data.indptr[0:size(H_eff.data.indptr)-1],
            H_eff.data.shape[0],H_eff.data.shape[1],
            H_eff.data.nnz,size(H_eff.data.indptr)-1)

def init_c_ops(c_ops):
    n = len(c_ops)
    first = True
    for i in range(n):
        datad = array(c_ops[i].data.data,dtype=wpc)
        qt.qutraj_run.init_c_ops(i+1,n,datad,c_ops[i].data.indices,
                c_ops[i].data.indptr[0:size(c_ops[i].data.indptr)-1],
                c_ops[i].data.shape[0],c_ops[i].data.shape[1],
                c_ops[i].data.nnz,size(c_ops[i].data.indptr)-1,first)
        first = False

def init_e_ops(e_ops):
    n = len(e_ops)
    first = True
    for i in range(n):
        datad = array(e_ops[i].data.data,dtype=wpc)
        qt.qutraj_run.init_e_ops(i+1,n,datad,e_ops[i].data.indices,
                e_ops[i].data.indptr[0:size(e_ops[i].data.indptr)-1],
                e_ops[i].data.shape[0],e_ops[i].data.shape[1],
                e_ops[i].data.nnz,size(e_ops[i].data.indptr)-1,first)
        first = False


def get_states(nstep,ntraj):
    states=array([array([Qobj()]*nstep)]*ntraj)
    for traj in range(ntraj):
        for i in range(nstep):
            states[traj][i] = Qobj(matrix(qt.qutraj_run.sol[0,traj,i]).transpose())
    return states

def get_expect(nstep,n_e_ops):
    expect=[array([0.+0.j]*nstep)]*n_e_ops
    for j in range(n_e_ops):
        expect[j] = qt.qutraj_run.sol[j,0,:,0]
    return expect


def finalize():
    qt.qutraj_run.finalize_work()
    qt.qutraj_run.finalize_sol()

def mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=500,options=Odeoptions()):
    init_tlist(tlist)
    init_psi0(psi0)
    init_hamilt(H,c_ops)
    if (c_ops == []):
        # force one trajectory if no collapse operators
        ntraj=1
    else:
        init_c_ops(c_ops)
    if (e_ops == []):
        states=True
    else:
        states=False
        init_e_ops(e_ops)
    qt.qutraj_run.ntraj = ntraj
    qt.qutraj_run.mc_avg = options.mc_avg
    qt.qutraj_run.init_odedata(psi0.shape[0],
            options.atol,options.rtol,options.max_step,mf=10)
    qt.qutraj_run.norm_steps=1
    #qt.qutraj_run.norm_tol=0.01
    qt.qutraj_run.evolve(states)
    sol = Odedata()
    sol.ntraj=ntraj
    sol.num_collapse=size(c_ops)
    sol.num_expect=size(e_ops)
    sol.solver='Fortran 90 Monte Carlo Solver'
    sol.times = tlist
    if (states):
        sol.states = get_states(size(tlist),ntraj)
    else:
        sol.expect = get_expect(size(tlist),size(e_ops))
    return sol




gamma = 1
neq = 2
psi0 = basis(neq,0)
#psi0 = Qobj([[1],[1]])
#a = destroy(neq)
#ad = a.dag()
#H = ad*a
#c_ops = [gamma*a]
#e_ops = [ad*a]
H = sigmax()
c_ops = [sqrt(gamma)*sigmax()]
#e_ops = [sigmam()*sigmap(),sigmap()*sigmam()]
#e_ops = [sigmam()*sigmap()]

# Times
T = 10.0
dt = 0.1
nstep = int(T/dt)
tlist = np.linspace(0,T,nstep)

ntraj=10

#init_tlist(tlist)
#init_psi0(psi0)
#init_hamilt(H,c_ops)
#init_c_ops(c_ops)
#init_e_ops(e_ops)

# set options
opts = Odeoptions()
opts.num_cpus=1
opts.gui=False
#opts.max_step=1000
#atol = opts.atol
#rtol = opts.rtol
#
#qt.qutraj_run.ntraj = ntraj
#
#qt.qutraj_run.init_odedata(neq,atol,rtol,mf=10)
#
#qt.qutraj_run.evolve()

#start_time = time.time()
#sol_f90 = mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
##sol_f90 = mcsolve_f90(H,psi0,tlist,[],e_ops,ntraj=ntraj,options=opts)
##sol_f90 = mcsolve_f90(H,psi0,tlist,c_ops,[],ntraj=ntraj,options=opts)
#print "solutiton took", time.time()-start_time, "s"

#start_time = time.time()
##sol= mesolve(H,psi0,tlist,[],e_ops)
##sol= mesolve(H,psi0,tlist,c_ops,[])
#sol_me = mesolve(H,psi0,tlist,c_ops,e_ops)
#print "solutiton took", time.time()-start_time, "s"
#start_time = time.time()
#sol_mc = mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
#sol = mcsolve(H,psi0,tlist,c_ops,[],ntraj=ntraj)
sol_mc = mcsolve(H,psi0,tlist,c_ops,[],ntraj=ntraj,options=opts)
#print "solutiton took", time.time()-start_time, "s"

#H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
#        for i in range(len(c_ops))])
##sol_me = mesolve(H_eff,psi0,tlist,[],e_ops)
#
#states = sol_f90.states
### avg over trajectories
##avg = zeros(size(tlist))
##for i in range(ntraj):
##    avg = avg+real(expect(e_ops[0],states[i]))
##avg = avg/ntraj
#
#sol_bare = qt.qutraj_run.sol
#
#plt.figure()
#for i in range(len(e_ops)):
#    #plt.plot(tlist,avg)
#    #plt.plot(tlist,sol_f90.expect[i])
#    plt.plot(tlist,sol_mc.expect[i],'--')
#    #plt.plot(tlist,sol_me.expect[i],'--')
#
#plt.plot(tlist,sol_f90.expect[0]+sol_f90.expect[1])
#plt.plot(tlist,sol.expect[0]+sol.expect[1],'--')


#finalize()
