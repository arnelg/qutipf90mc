import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import qutipf90mc as mcf90
import time

def test():
    gamma = 1
    neq = 2
    psi0 = qt.basis(2,0)
    #a = destroy(neq)
    #ad = a.dag()
    #H = ad*a
    #c_ops = [gamma*a]
    #e_ops = [ad*a]
    H = qt.sigmax()
    c_ops = [np.sqrt(gamma)*qt.sigmax()]
    #e_ops = [qt.sigmam()*qt.sigmap(),qt.sigmap()*qt.sigmam()]
    e_ops = []

    #neq = 15
    #psi0 = qt.tensor(qt.basis(neq,0),qt.basis(neq,0),qt.basis(neq,0))
    #H = qt.tensor(qt.qeye(neq),qt.qeye(neq),qt.qeye(neq))
    #c_ops = []
    #e_ops = []

    #a0=qt.tensor(qt.destroy(neq),qt.qeye(neq),qt.qeye(neq))
    #a1=qt.tensor(qt.qeye(neq),qt.destroy(neq),qt.qeye(neq))
    #a2=qt.tensor(qt.qeye(neq),qt.qeye(neq),qt.destroy(neq))

    ##initial state: coherent mode 0 & vacuum for modes #1 & #2
    #alpha=np.sqrt(7)#initial coherent state param for mode 0
    #psi0=qt.tensor(qt.coherent(neq,alpha),qt.basis(neq,0),qt.basis(neq,0))

    ##trilinear Hamiltonian
    #H=1.0j*(a0*a1.dag()*a2.dag()-a0.dag()*a1*a2)


    # Times
    T = 1.0
    dt = 0.1
    nstep = int(T/dt)
    tlist = np.linspace(0,T,nstep)

    #tlist=np.linspace(0,2.5,50)

    ntraj=1

    # set options
    opts = qt.Odeoptions()
    opts.num_cpus=1
    opts.mc_avg = True
    #opts.gui=False
    #opts.max_step=1000
    #opts.atol =
    #opts.rtol =

    sol_f90 = qt.Odedata()
    start_time = time.time()
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve_f90 solutiton took", time.time()-start_time, "s"

    start_time = time.time()
    sol_me = qt.mesolve(H,psi0,tlist,c_ops,e_ops,options=opts)
    print "mesolve solutiton took", time.time()-start_time, "s"

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

def testdemos():
    import qutipf90mc.examples as examples
    print 'running demo #30 from qutip'
    raw_input('press a key to continue')
    plt.figure()
    ex_code = compile('examples.ex_30.run()','<string>','exec')
    eval(ex_code)
    print 'running demo #31 from qutip'
    raw_input('press a key to continue')
    plt.figure()
    ex_code = compile('examples.ex_31.run()','<string>','exec')
    eval(ex_code)
    print 'running demo #33 from qutip'
    raw_input('press a key to continue')
    #plt.figure()
    ex_code = compile('examples.ex_33.run()','<string>','exec')
    eval(ex_code)

if __name__ == '__main__':
    test()
