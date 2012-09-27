import numpy as np
from qutip import *
import time

# Working precision
wpr = dtype(float64)
wpc = dtype(complex64)

def mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=500,
        options=Odeoptions()):
    from multiprocessing import cpu_count
    mc = _MC_class()
    mc.H = H
    mc.psi0 = psi0
    mc.tlist = tlist
    mc.c_ops = c_ops
    mc.e_ops = e_ops
    mc.ntraj = ntraj
    mc.options = options
    if (options.num_cpus == 0):
        mc.ncpus = cpu_count()
    else:
        mc.ncpus = options.num_cpus
    mc.ncpus=4
    mc.sols = [Odedata()]*mc.ncpus
    mc.run()
    return mc.sol

class _MC_class():
    def __init__(self):
        self.H = Qobj()
        self.psi0 = Qobj()
        self.tlist = []
        self.c_ops = [Qobj()]
        self.e_ops = [Qobj()]
        self.ntraj = 0
        self.options = Odeoptions()
        self.ncpus = 0
        self.sols = [Odedata()]
        self.sol = [Odedata()]
        self.states = True

    def parallel(self):
        from multiprocessing import Process, Queue
        processes = []
        queue = Queue()
        ntrajs = []
        for i in range(self.ncpus):
            ntrajs.append(min(int(ceil(float(self.ntraj)/self.ncpus)),
                self.ntraj-sum(ntrajs)))
        print "running in parallell on ", self.ncpus, " cpus."
        print "number of trajectories for each process:"
        print ntrajs
        for i in range(self.ncpus):
            nt = ntrajs[i]
            p = Process(target=self.solve_serial,
                    args=((queue,ntrajs[i],i),))
            p.start()
            processes.append(p)
        for i,proc in enumerate(processes):
            try:
                proc.join()
                self.sols[i] = queue.get()
            except KeyboardInterrupt:
                print("Cancel thread on keyboard interrupt")
                proc.terminate()
                proc.join()
        return

    def run(self):
        if (self.c_ops == []):
            # force one trajectory if no collapse operators
            self.ntraj=1
        if (self.e_ops == []):
            self.states=True
        else:
            self.states=False
        # run in paralell
        self.parallel()
        # gather data
        self.sol = Odedata()
        self.sol.ntraj=self.ntraj
        self.sol.num_collapse=size(self.c_ops)
        self.sol.num_expect=size(self.e_ops)
        self.sol.solver='Fortran 90 Monte Carlo Solver'
        self.sol.times = self.tlist
        if (self.states):
            pass
            #self.sol.states = 
        else:
            self.sol.expect = self.sols[0].expect
            for i in range(size(self.e_ops)):
                for j in range(1,self.ncpus):
                    self.sol.expect[i] += self.sols[j].expect[i]
                self.sol.expect[i] = self.sol.expect[i]/self.ncpus

    def solve_serial(self,args):
        # run ntraj trajectories in parallell via fortran
        import qutraj_run as qt
        # get args
        queue,ntraj,rngseed = args
        def init_tlist(tlist):
            qt.qutraj_run.init_tlist(array(tlist,dtype=wpr))
        def init_psi0(psi0):
            psi0d = array(psi0.data.toarray(),dtype=wpc).transpose()
            qt.qutraj_run.init_psi0(psi0d,size(psi0d))
        def init_hamilt(H,c_ops):
            # construct effective non-Hermitian Hamiltonian
            H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
                for i in range(len(c_ops))])
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

        # Initialize stuff
        init_tlist(self.tlist)
        init_psi0(self.psi0)
        init_hamilt(self.H,self.c_ops)
        if (self.c_ops != []):
            init_c_ops(self.c_ops)
        if (self.e_ops != []):
            init_e_ops(self.e_ops)
        # set options
        qt.qutraj_run.ntraj = ntraj
        qt.qutraj_run.mc_avg = self.options.mc_avg
        qt.qutraj_run.init_odedata(self.psi0.shape[0],
                self.options.atol,self.options.rtol,self.options.max_step,
                mf=10)
        #qt.qutraj_run.norm_steps=1
        #qt.qutraj_run.norm_tol=0.01

        #run
        qt.qutraj_run.evolve(self.states,rngseed)

        # construct Odedata instance
        sol = Odedata()
        sol.ntraj=ntraj
        sol.num_collapse=size(self.c_ops)
        sol.num_expect=size(self.e_ops)
        sol.solver='Fortran 90 Monte Carlo Solver'
        sol.times = self.tlist
        if (self.states):
            sol.states = get_states(size(self.tlist),ntraj)
        else:
            sol.expect = get_expect(size(self.tlist),size(self.e_ops))
        queue.put(sol)
        return
