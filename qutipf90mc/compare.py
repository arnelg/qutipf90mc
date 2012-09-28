import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import qutipf90mc as mcf90
import time

def run(neq,ntraj,f90only=False):
    gamma = 1
    # sparse initial state
    #psi0 = basis(neq,neq-1)
    # dense initial state
    psi0 = Qobj(ones((neq,1))).unit()
    a = destroy(neq)
    ad = a.dag()
    H = ad*a
    #c_ops = [gamma*a]
    c_ops = [qeye(neq)]
    e_ops = [ad*a]

    # Times
    T = 10.0
    dt = 0.1
    nstep = int(T/dt)
    tlist = np.linspace(0,T,nstep)

    # set options
    opts = Odeoptions()
    opts.num_cpus = 1

    start_time = time.time()
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve_f90 solutiton took", time.time()-start_time, "s"

    if (not f90only):
        start_time = time.time()
        sol_mc = mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
        print "mcsolve solutiton took", time.time()-start_time, "s"

if __name__ == '__main__':
    run(10,100)
