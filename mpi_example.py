"""
This examplifies how to use MPI with qutipf90mc. It requires MPI and the 
python package mpi4py. You can install e.g. openmpi on ubuntu by

    sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev

Then mpi4py can be installed in the usual way you install python packages. 
For example with pip

    pip install mpi4py

To run this example on 4 CPUS:

    mpiexec -n 4 python mpi_example.py
"""

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import qutipf90mc as mcf90

def run():
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    print "Process number", rank, "of", size, "total."

    neq = 2
    gamma = 1.0
    psi0 = qt.basis(neq,neq-1)
    H = qt.sigmax()
    c_ops = [np.sqrt(gamma)*qt.sigmax()]
    e_ops = [qt.sigmam()*qt.sigmap()]

    tlist = np.linspace(0,10,100)

    ntraj=1000

    # One CPU per MPI process
    opts = qt.Odeoptions()
    opts.num_cpus = 1

    # Solve
    sols = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    # Gather data
    sols = comm.gather(sols,root=0)
    if (rank==0):
        sol = sols[0]
        sol.expect = np.array(sols[0].expect)
        plt.figure()
        plt.plot(tlist,sols[0].expect[0],'r',label='proc '+str(0))
        for i in range(1,size):
            plt.plot(tlist,sols[i].expect[0],'r',label='proc '+str(i))
            sol.expect += np.array(sols[i].expect)
        sol.expect = sol.expect/size
        plt.plot(tlist,sol.expect[0],'b',label='average')
        plt.legend()
        plt.show()

if (__name__ == '__main__'):
    run()

