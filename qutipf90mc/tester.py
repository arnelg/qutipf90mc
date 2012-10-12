import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import qutipf90mc as mcf90
import time

def ptracetest():
    gamma = 1.
    neq = 2
    psi0 = qt.basis(neq,neq-1)
    psi0 = qt.tensor(psi0,psi0)
    H = qt.tensor(qt.sigmax(),qt.sigmay())
    c1 = np.sqrt(gamma)*qt.sigmax()
    e1 = np.sqrt(gamma)*qt.sigmaz()
    c_ops = [qt.tensor(c1,c1)]
    e_ops = [qt.tensor(e1,e1),qt.tensor(c1,c1)]
    #e_ops = []
    tlist = np.linspace(0,10,100)
    ntraj = 2000
    ptrace_sel = [0]
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,
            ptrace_sel=ptrace_sel,calc_entropy=True)
    #sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj)
    #sol_me = qt.mesolve(H,psi0,tlist,c_ops,e_ops)
    #exp_f90 = qt.expect(e1,sol_f90.states)
    #exp_me = qt.expect(qt.tensor(e1,qt.qeye(neq)),sol_me.states)
    #plt.figure()
    #plt.plot(tlist,exp_me)
    #plt.plot(tlist,exp_f90)
    #return sol_f90#,sol_me

def test():
    gamma = 1.
    neq = 2
    psi0 = qt.basis(neq,neq-1)
    #a = qt.destroy(neq)
    #ad = a.dag()
    #H = ad*a
    #c_ops = [gamma*a]
    #e_ops = [ad*a]
    H = qt.sigmax()
    c_ops = [np.sqrt(gamma)*qt.sigmax()]
    #c_ops = []
    e_ops = [qt.sigmam()*qt.sigmap(),qt.sigmap()*qt.sigmam()]
    #e_ops = []

    # Times
    T = 2.0
    dt = 0.1
    nstep = int(T/dt)
    tlist = np.linspace(0,T,nstep)

    ntraj=100

    # set options
    opts = qt.Odeoptions()
    opts.num_cpus=2
    #opts.mc_avg = True
    #opts.gui=False
    #opts.max_step=1000
    #opts.atol =
    #opts.rtol =

    sol_f90 = qt.Odedata()
    start_time = time.time()
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve_f90 solutiton took", time.time()-start_time, "s"

    sol_me = qt.Odedata()
    start_time = time.time()
    sol_me = qt.mesolve(H,psi0,tlist,c_ops,e_ops,options=opts)
    print "mesolve solutiton took", time.time()-start_time, "s"

    sol_mc = qt.Odedata()
    start_time = time.time()
    sol_mc = qt.mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve solutiton took", time.time()-start_time, "s"

    if (e_ops == []):
        e_ops = [qt.sigmam()*qt.sigmap(),qt.sigmap()*qt.sigmam()]
        sol_f90expect = [np.array([0.+0.j]*nstep)]*len(e_ops)
        sol_mcexpect = [np.array([0.+0.j]*nstep)]*len(e_ops)
        sol_meexpect = [np.array([0.+0.j]*nstep)]*len(e_ops)
        for i in range(len(e_ops)):
            if (not opts.mc_avg):
                sol_f90expect[i] = sum([qt.expect(e_ops[i],
                    sol_f90.states[j]) for j in range(ntraj)])/ntraj
                sol_mcexpect[i] = sum([qt.expect(e_ops[i],
                    sol_mc.states[j]) for j in range(ntraj)])/ntraj
            else:
                sol_f90expect[i] = qt.expect(e_ops[i],sol_f90.states)
                sol_mcexpect[i] = qt.expect(e_ops[i],sol_mc.states)
            sol_meexpect[i] = qt.expect(e_ops[i],sol_me.states)
    elif (not opts.mc_avg):
        sol_f90expect = sum(sol_f90.expect,0)/ntraj
        sol_mcexpect = sum(sol_f90.expect,0)/ntraj
        sol_meexpect = sol_me.expect
    else:
        sol_f90expect = sol_f90.expect
        sol_mcexpect = sol_mc.expect
        sol_meexpect = sol_me.expect

    plt.figure()
    for i in range(len(e_ops)):
        plt.plot(tlist,sol_f90expect[i],'b')
        plt.plot(tlist,sol_mcexpect[i],'g')
        plt.plot(tlist,sol_meexpect[i],'k')

    return sol_f90, sol_mc

def rundemo(no,fig=False):
    """ Run a demo adapted from qutip

    Parameters
    ----------
    no : int
        Demo number. Available demos are:
        24 - Dynamics of a Heisenberg spin chain
        30 - MC Cavity+Qubig
        31 - Coupled Oscillators
        33 - Trilinear Hamiltonian
        34 - Visualize MC Dissipation
    fig : boolean
        Open a new pyplot.figure before running demo? Default is False.

    See also qutipf90mc.alldemos.
    """
    import qutipf90mc.examples as examples
    print 'running demo #',str(no),'from qutip'
    raw_input('press a key to continue')
    funcstr = 'examples.ex_'+str(no)+'.run()'
    ex_code = compile(funcstr,'<string>','exec')
    if (fig): plt.figure()
    eval(ex_code)

def alldemos():
    """
    Run all available demos. See also qutipf90mc.rundemo.
    """
    rundemo(24,fig=True)
    rundemo(30,fig=True)
    rundemo(31,fig=True)
    rundemo(33)
    rundemo(34)

if __name__ == '__main__':
    ptracetest()
