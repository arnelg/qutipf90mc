import numpy as np
import matplotlib.pyplot as plt
from qutip import *
import qutipf90mc as mcf90
import time

def run(neq=10,ntraj=100,solver='both',ncpus=1):
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
    opts.num_cpus = ncpus
    opts.gui = False

    mcf90_time = 0.
    mc_time = 0.
    if (solver=='mcf90' or solver=='both'):
        start_time = time.time()
        sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)

        mcf90_time = time.time()-start_time
        print "mcsolve_f90 solutiton took", mcf90_time, "s"

    if (solver=='mc' or solver=='both'):
        start_time = time.time()
        sol_mc = mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
        mc_time = time.time()-start_time
        print "mcsolve solutiton took", mc_time, "s"
    return mcf90_time, mc_time

def compare_system_size_scaling(dims=[10,20,30,40,50,60,70,80,90,100]):
    #dims = [10,100,200,300,400,500,600,700,800,900,1000]
    #dims = [10,100,200,300,400,500,1000,2000,4000,6000,8000]
    ntraj = 100
    t1 = [0.]*size(dims)
    t2 = [0.]*size(dims)
    for i,d in enumerate(dims):
        print 'dimension=',dims[i]
        times = run(d,ntraj)
        t1[i] = times[0]
        t2[i] = times[1]
    out_data = np.vstack((dims,t1,t2))
    file_data_store('compare_system_size_scaling.dat', out_data.T)
    plt.figure()
    plt.plot(dims,t1,label='fortran')
    plt.plot(dims,t2,label='python')
    plt.legend()
    plt.ylabel('Computation Time (sec)')
    plt.xlabel('Hilbert Space Dimension')
    plt.figure()
    plt.plot(dims,[t2[i]/t1[i] for i in range(size(t1))])
    plt.ylabel('(Python Comp. Time)/(Fortran Comp. Time)')
    plt.xlabel('Hilbert Space Dimension')

def plot_system_size_scaling():
    data =file_data_read('compare_system_size_scaling.dat')
    dims = data[:,0]
    t1 = data[:,1]
    t2 = data[:,2]
    plt.figure()
    plt.plot(dims,t1,label='fortran')
    plt.plot(dims,t2,label='python')
    plt.legend()
    plt.ylabel('Computation Time (sec)')
    plt.xlabel('Hilbert Space Dimension')
    plt.figure()
    plt.plot(dims,[t2[i]/t1[i] for i in range(size(t1))])
    plt.ylabel('(Python Comp. Time)/(Fortran Comp. Time)')
    plt.xlabel('Hilbert Space Dimension')

if __name__ == '__main__':
    run(10,100)
