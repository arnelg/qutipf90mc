import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import qutipf90mc as mcf90
import time

def test():
    gamma = 1
    neq = 2
    psi0 = qt.basis(neq,0)
    #psi0 = Qobj([[1],[1]])
    #a = destroy(neq)
    #ad = a.dag()
    #H = ad*a
    #c_ops = [gamma*a]
    #e_ops = [ad*a]
    H = qt.sigmax()
    c_ops = [np.sqrt(gamma)*qt.sigmax()]
    e_ops = [qt.sigmam()*qt.sigmap(),qt.sigmap()*qt.sigmam()]
    #e_ops = [sigmam()*sigmap()]

    # Times
    T = 10.0
    dt = 0.1
    nstep = int(T/dt)
    tlist = np.linspace(0,T,nstep)

    ntraj=1000

    # set options
    opts = qt.Odeoptions()
    #opts.num_cpus=1
    #opts.gui=False
    #opts.max_step=1000
    #opts.atol =
    #opts.rtol =

    start_time = time.time()
    sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve_f90 solutiton took", time.time()-start_time, "s"

    start_time = time.time()
    sol_me = qt.mesolve(H,psi0,tlist,c_ops,e_ops,options=opts)
    print "mesolve solutiton took", time.time()-start_time, "s"

    start_time = time.time()
    sol_mc = qt.mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
    print "mcsolve solutiton took", time.time()-start_time, "s"

    plt.figure()
    for i in range(len(e_ops)):
        plt.plot(tlist,sol_f90.expect[i],'b')
        plt.plot(tlist,sol_mc.expect[i],'g')
        plt.plot(tlist,sol_me.expect[i],'k')

    return sol_f90

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

if __name__ == '__main__':
    test()
