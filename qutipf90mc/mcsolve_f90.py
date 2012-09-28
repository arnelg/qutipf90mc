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
    if (options.method == 'adams'):
        mc.mf = 10
    else:
        print 'support for stiff "bdf"-method not implemented, using "adams".'
        mc.mf = 10
    if (options.num_cpus == 0):
        mc.ncpus = cpu_count()
    else:
        mc.ncpus = options.num_cpus
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
        self.mf = 10

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
        # run ntraj trajectories in for one process via fortran
        import qutraj_run as qt
        # get args
        queue,ntraj,instanceno = args
        # initalizers
        def init_tlist(tlist):
            qt.qutraj_run.init_tlist(_realarray_to_fortran(tlist),
                    size(tlist))
        def init_psi0(psi0):
            psi0d = _qobj_to_fortranfull(psi0)
            qt.qutraj_run.init_psi0(psi0d,size(psi0d))
        def init_hamilt(H,c_ops):
            # construct effective non-Hermitian Hamiltonian
            H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
                for i in range(len(c_ops))])
            Of = _qobj_to_fortrancsr(H_eff)
            qt.qutraj_run.init_hamiltonian(Of[0],Of[1],
                    Of[2],Of[3],Of[4],Of[5],Of[6])
        def init_c_ops(c_ops):
            n = len(c_ops)
            first = True
            for i in range(n):
                Of = _qobj_to_fortrancsr(c_ops[i])
                qt.qutraj_run.init_c_ops(i,n,Of[0],Of[1],
                        Of[2],Of[3],Of[4],Of[5],Of[6],first)
                first = False
        def init_e_ops(e_ops):
            n = len(e_ops)
            first = True
            for i in range(n):
                Of = _qobj_to_fortrancsr(e_ops[i])
                qt.qutraj_run.init_e_ops(i,n,Of[0],Of[1],
                        Of[2],Of[3],Of[4],Of[5],Of[6],first)
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
                self.options.atol,self.options.rtol,mf=self.mf)
        # set optional arguments
        qt.qutraj_run.order = self.options.order
        qt.qutraj_run.nsteps = self.options.nsteps
        qt.qutraj_run.first_step = self.options.first_step
        qt.qutraj_run.min_step = self.options.min_step
        qt.qutraj_run.max_step = self.options.max_step
        qt.qutraj_run.norm_steps=self.options.norm_steps
        qt.qutraj_run.norm_tol=self.options.norm_tol

        #run
        qt.qutraj_run.evolve(self.states,instanceno)

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
        # put to queue
        queue.put(sol)
        return

#
# Misc. converison functions
#

def _realarray_to_fortran(a):
    datad = np.asfortranarray(np.array(a,dtype=wpr))
    return datad

def _qobj_to_fortranfull(A):
    datad = np.asfortranarray(np.array(A.data.toarray(),dtype=wpc))
    return datad

def _qobj_to_fortrancsr(A):
    datad = np.asfortranarray(np.array(A.data.data,dtype=wpc))
    indices = np.asfortranarray(np.array(A.data.indices))
    indptr = np.asfortranarray(np.array(A.data.indptr[0:size(A.data.indptr)-1]))
    m = A.data.shape[0]
    k = A.data.shape[1]
    nnz = A.data.nnz
    nptr = size(A.data.indptr)-1
    return datad,indices,indptr,m,k,nnz,nptr

