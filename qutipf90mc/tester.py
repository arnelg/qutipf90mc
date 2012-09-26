import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import qutipf90mc as mcf90
import time

def test():
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
    e_ops = [sigmam()*sigmap(),sigmap()*sigmam()]
    #e_ops = [sigmam()*sigmap()]

    # Times
    T = 10.0
    dt = 0.1
    nstep = int(T/dt)
    tlist = np.linspace(0,T,nstep)

    ntraj=10

    # set options
    opts = Odeoptions()
    #opts.num_cpus=1
    #opts.gui=False
    #opts.max_step=1000
    #opts.atol =
    #opts.rtol =

    start_time = time.time()
    #sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,[],ntraj=ntraj,options=opts)
    print "solutiton took", time.time()-start_time, "s"

    start_time = time.time()
    sol_me = mesolve(H,psi0,tlist,c_ops,e_ops,options=opts)
    print "solutiton took", time.time()-start_time, "s"


    #plt.figure()
    #for i in range(len(e_ops)):
    #    plt.plot(tlist,sol_f90.expect[i])
    #    #plt.plot(tlist,sol_mc.expect[i],'--')
    #    plt.plot(tlist,sol_me.expect[i],'--')

    return sol_f90

