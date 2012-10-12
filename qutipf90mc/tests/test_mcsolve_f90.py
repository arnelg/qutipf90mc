#
# Unit tests for mcsolve_f90
# Adapdted from qutip/tests/test_mcsolve.py
#

from qutip import *
from qutip.odechecks import _ode_checks
from numpy import allclose
from numpy.testing import assert_equal
from numpy.testing.decorators import skipif
from qutipf90mc import mcsolve_f90
import unittest
#find Cython if it exists
try:
    import Cython
except:
    Cython_found=0
else:
    Cython_found=1

kappa=0.2
def sqrt_kappa(t,args):
    return sqrt(kappa)

def sqrt_kappa2(t,args):
    return sqrt(kappa*exp(-t))

def const_H1_coeff(t,args):
    return 0.0

#average error for failure
mc_error=5e-2 #5% for ntraj=500

def test_MCNoCollExpt():
    "Monte-carlo: Constant H with no collapse ops (expect)"
    error=1e-8
    N=10 #number of basis states to consider
    a=destroy(N)
    H=a.dag()*a
    psi0=basis(N,9) #initial state
    kappa=0.2 #coupling to oscillator
    c_op_list=[]
    tlist=linspace(0,10,100)
    mcdata=mcsolve_f90(H,psi0,tlist,c_op_list,[a.dag()*a],options=Odeoptions(gui=False))
    expt=mcdata.expect[0]
    actual_answer=9.0*ones(len(tlist))
    diff=mean(abs(actual_answer-expt)/actual_answer)
    assert_equal(diff<error,True)


def test_MCNoCollStates():
    "Monte-carlo: Constant H with no collapse ops (states)"
    error=1e-8
    N=10 #number of basis states to consider
    a=destroy(N)
    H=a.dag()*a
    psi0=basis(N,9) #initial state
    kappa=0.2 #coupling to oscillator
    c_op_list=[]
    tlist=linspace(0,10,100)
    mcdata=mcsolve_f90(H,psi0,tlist,c_op_list,[],options=Odeoptions(gui=False))
    states=mcdata.states
    expt=expect(a.dag()*a,states)
    actual_answer=9.0*ones(len(tlist))
    diff=mean(abs(actual_answer-expt)/actual_answer)
    assert_equal(diff<error,True)

def test_MCSimpleConst():
    "Monte-carlo: Constant H with constant collapse"
    N=10 #number of basis states to consider
    a=destroy(N)
    H=a.dag()*a
    psi0=basis(N,9) #initial state
    kappa=0.2 #coupling to oscillator
    c_op_list=[sqrt(kappa)*a]
    tlist=linspace(0,10,100)
    mcdata=mcsolve_f90(H,psi0,tlist,c_op_list,[a.dag()*a],options=Odeoptions(gui=False))
    expt=mcdata.expect[0]
    actual_answer=9.0*exp(-kappa*tlist)
    avg_diff=mean(abs(actual_answer-expt)/actual_answer)
    assert_equal(avg_diff<mc_error,True)
