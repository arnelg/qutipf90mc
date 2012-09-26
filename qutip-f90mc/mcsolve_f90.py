import numpy as np
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
