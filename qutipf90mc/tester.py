import numpy as np
import matplotlib.pyplot as plt
from qutip import *
#import qutraj_run as qt
import mcsolve_f90 as mcf90
import time

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

#init_tlist(tlist)
#init_psi0(psi0)
#init_hamilt(H,c_ops)
#init_c_ops(c_ops)
#init_e_ops(e_ops)

# set options
opts = Odeoptions()
#opts.num_cpus=1
#opts.gui=False
#opts.max_step=1000
#atol = opts.atol
#rtol = opts.rtol
#
#qt.qutraj_run.ntraj = ntraj
#
#qt.qutraj_run.init_odedata(neq,atol,rtol,mf=10)
#
#qt.qutraj_run.evolve()

start_time = time.time()
sol_f90 = mcf90.mcsolve_f90(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
##sol_f90 = mcsolve_f90(H,psi0,tlist,[],e_ops,ntraj=ntraj,options=opts)
##sol_f90 = mcsolve_f90(H,psi0,tlist,c_ops,[],ntraj=ntraj,options=opts)
print "solutiton took", time.time()-start_time, "s"

start_time = time.time()
##sol= mesolve(H,psi0,tlist,[],e_ops)
##sol= mesolve(H,psi0,tlist,c_ops,[])
sol_me = mesolve(H,psi0,tlist,c_ops,e_ops)
#print "solutiton took", time.time()-start_time, "s"
#start_time = time.time()
#sol_mc = mcsolve(H,psi0,tlist,c_ops,e_ops,ntraj=ntraj,options=opts)
#sol = mcsolve(H,psi0,tlist,c_ops,[],ntraj=ntraj)
#sol_mc = mcsolve(H,psi0,tlist,c_ops,[],ntraj=ntraj,options=opts)
print "solutiton took", time.time()-start_time, "s"

#H_eff = H - 0.5j*sum([c_ops[i].dag()*c_ops[i]
#        for i in range(len(c_ops))])
##sol_me = mesolve(H_eff,psi0,tlist,[],e_ops)
#
#states = sol_f90.states
### avg over trajectories
##avg = zeros(size(tlist))
##for i in range(ntraj):
##    avg = avg+real(expect(e_ops[0],states[i]))
##avg = avg/ntraj
#
#sol_bare = qt.qutraj_run.sol
#
plt.figure()
for i in range(len(e_ops)):
    #plt.plot(tlist,avg)
    plt.plot(tlist,sol_f90.expect[i])
    #plt.plot(tlist,sol_mc.expect[i],'--')
    plt.plot(tlist,sol_me.expect[i],'--')
#
#plt.plot(tlist,sol_f90.expect[0]+sol_f90.expect[1])
#plt.plot(tlist,sol.expect[0]+sol.expect[1],'--')


#finalize()
