import numpy as np
from qutip import *

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
    mc.nprocs = mc.ncpus
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
        self.nprocs = 0
        self.sols = [Odedata()]
        #self.sols = []
        self.sol = [Odedata()]
        self.states = True
        self.mf = 10
        # If returning states, return averaged kets or density matrices?
        self.return_kets = True
        # If returning density matrices, return as sparse or dense?
        self.rho_return_sparse = True

    def parallel(self):
        from multiprocessing import Process, Queue, JoinableQueue
        ntrajs = []
        for i in range(self.ncpus):
            ntrajs.append(min(int(floor(float(self.ntraj)/self.ncpus)),
                self.ntraj-sum(ntrajs)))
        cnt = sum(ntrajs)
        while cnt<self.ntraj:
            for i in range(self.ncpus):
                ntrajs[i] += 1
                cnt+=1
                if (cnt>=self.ntraj):
                    break
        ntrajs = np.array(ntrajs)
        self.nprocs = len(ntrajs[np.where(ntrajs>0)])
        sols = []
        processes = []
        resq = JoinableQueue()
        print "Number of cpus:", self.ncpus
        print "Trying to start", self.nprocs, "process(es)."
        print "Number of trajectories for each process:"
        print ntrajs
        for i in range(self.nprocs):
            nt = ntrajs[i]
            p = Process(target=self.evolve_serial,
                    args=((resq,ntrajs[i],i),))
            p.start()
            processes.append(p)
        resq.join()
        cnt = 0
        while True:
            try:
                sols.append(resq.get())
                resq.task_done()
                cnt += 1
                if (cnt >= self.nprocs): break
            except KeyboardInterrupt:
                break
            except:
                pass
        resq.join()
        for proc in processes:
            try:
                proc.join()
            except KeyboardInterrupt:
                print("Cancel thread on keyboard interrupt")
                proc.terminate()
                proc.join()
        resq.close()
        return sols

    def run(self):
        if (self.c_ops == []):
            # force one trajectory if no collapse operators
            self.ntraj=1
        if (self.e_ops == []):
            self.states=True
        else:
            self.states=False
        # run in paralell
        sols = self.parallel()
        # construct Odedata object
        self.sol = Odedata()
        self.sol.ntraj=self.ntraj
        self.sol.num_collapse=size(self.c_ops)
        self.sol.num_expect=size(self.e_ops)
        self.sol.solver='Fortran 90 Monte Carlo Solver'
        self.sol.times = self.tlist
        # gather data
        if (self.states):
            if (self.options.mc_avg):
                # collect states, averaged over trajectories
                self.sol.states = sols[0].states
                for j in range(1,self.nprocs):
                    self.sol.states += sols[j].states
                self.sol.states = self.sol.states/self.nprocs
                # convert to list to be consistent with qutip mcsolve
                self.sol.states = list(self.sol.states)
            else:
                # collect states, all trajectories
                self.sol.states = np.concatenate([sols[j].states 
                    for j in range(self.nprocs)])
        else:
            if (self.options.mc_avg):
                # collect expectation values, averaged
                self.sol.expect = sols[0].expect
                for i in range(size(self.e_ops)):
                    for j in range(1,self.nprocs):
                        self.sol.expect[i] += sols[j].expect[i]
                    self.sol.expect[i] = self.sol.expect[i]/self.nprocs
            else:
                # collect expectation values, all trajectories
                #self.sol.expect = sols[0].expect
                self.sol.expect = np.concatenate([sols[j].expect 
                    for j in range(self.nprocs)])
                # convert to list to be consistent with qutip mcsolve
                self.sol.expect = list(self.sol.expect)

    def evolve_serial(self,args):
        # run ntraj trajectories for one process via fortran
        import qutraj_run as qt
        # get args
        queue,ntraj,instanceno = args
        # Initalizers
        def init_tlist(tlist):
            Of = _realarray_to_fortran(tlist)
            qt.qutraj_run.init_tlist(Of,
                    size(Of))
        def init_psi0(psi0):
            Of = _qobj_to_fortranfull(psi0)
            qt.qutraj_run.init_psi0(Of,size(Of))
        def init_hamilt(H,c_ops):
            # construct effective non-Hermitian Hamiltonian
            H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
                for i in range(len(c_ops))])
            Of = _qobj_to_fortrancsr(H_eff)
            qt.qutraj_run.init_hamiltonian(Of[0],Of[1],
                    Of[2],Of[3],Of[4])
        def init_c_ops(c_ops):
            n = len(c_ops)
            first = True
            for i in range(n):
                Of = _qobj_to_fortrancsr(c_ops[i])
                qt.qutraj_run.init_c_ops(i+1,n,Of[0],Of[1],
                        Of[2],Of[3],Of[4],first)
                first = False
        def init_e_ops(e_ops):
            n = len(e_ops)
            first = True
            for i in range(n):
                Of = _qobj_to_fortrancsr(e_ops[i])
                qt.qutraj_run.init_e_ops(i+1,n,Of[0],Of[1],
                        Of[2],Of[3],Of[4],first)
                first = False
        def finalize():
            # not in use...
            qt.qutraj_run.finalize_work()
            qt.qutraj_run.finalize_sol()
        # get data
        def get_states(nstep,ntraj):
            from scipy.sparse import csr_matrix
            if (self.options.mc_avg):
                states=np.array([Qobj()]*nstep)
                if (self.return_kets):
                    for i in range(nstep):
                        states[i] = Qobj(matrix(
                            qt.qutraj_run.sol[0,0,i,:]).transpose(),
                            dims=self.psi0.dims)
                elif (self.rho_return_sparse):
                    for i in range(nstep):
                        qt.qutraj_run.get_rho_sparse(i+1)
                        val = qt.qutraj_run.csr_val
                        col = qt.qutraj_run.csr_col-1
                        ptr = qt.qutraj_run.csr_ptr-1
                        m = qt.qutraj_run.csr_nrows
                        k = qt.qutraj_run.csr_ncols
                        states[i] = Qobj(csr_matrix((val,col,ptr),
                            (m,k)).toarray(),
                            dims=self.psi0.dims)
                else:
                    for i in range(nstep):
                        states[i] = Qobj(qt.qutraj_run.sol[0,i,:,:],
                            dims=self.psi0.dims)
            else:
                states=np.array([np.array([Qobj()]*nstep)]*ntraj)
                for traj in range(ntraj):
                    for i in range(nstep):
                        states[traj][i] = Qobj(matrix(
                            qt.qutraj_run.sol[0,traj,i,:]).transpose(),
                            dims=self.psi0.dims)
            return states
        def get_expect(nstep,n_e_ops):
            if (self.options.mc_avg):
                expect=[np.array([0.+0.j]*nstep)]*n_e_ops
                for j in range(n_e_ops):
                    expect[j] = qt.qutraj_run.sol[j,0,:,0]
            else:
                expect=np.array([[np.array([0.+0.j]*nstep)]*n_e_ops]*ntraj)
                for j in range(n_e_ops):
                    expect[:,j,:] = qt.qutraj_run.sol[j,:,:,0]
            return expect

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
        # how to return solution
        qt.qutraj_run.return_kets = self.return_kets
        qt.qutraj_run.rho_return_sparse = self.rho_return_sparse

        #run
        qt.qutraj_run.evolve(self.states,instanceno)

        # construct Odedata instance
        sol = Odedata()
        if (self.states):
            sol.states = get_states(size(self.tlist),ntraj)
        else:
            sol.expect = get_expect(size(self.tlist),size(self.e_ops))
        # put to queue
        queue.put(sol)
        #self.sols[instanceno] = sol
        #queue.put('STOP')
        #deallocate stuff
        #finalize()
        return

#
# Misc. converison functions
#

def _realarray_to_fortran(a):
    datad = np.array(a,dtype=wpr)
    return datad

def _qobj_to_fortranfull(A):
    datad = np.array(A.data.toarray(),dtype=wpc)
    return datad

def _qobj_to_fortrancsr(A):
    data = np.array(A.data.data,dtype=wpc)
    indices = np.array(A.data.indices)
    indptr = np.array(A.data.indptr)
    m = A.data.shape[0]
    k = A.data.shape[1]
    return data,indices+1,indptr+1,m,k

