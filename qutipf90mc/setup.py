#!/usr/bin/env python
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('qutipf90mc', parent_package, top_path)

    blas_opt = get_info('blas_opt',notfound_action=2)

    config.add_library('zvode',
            sources=[join('zvode','*.f')])

    libs = [
            'zvode',
            ]

    # Remove libraries key from blas_opt
    if 'libraries' in blas_opt: # key doesn't exist on OS X ...
        libs.extend(blas_opt['libraries'])
    newblas = {}
    for key in blas_opt.keys():
        if key == 'libraries':
            continue
        newblas[key] = blas_opt[key]


    config.add_extension('qutraj_run',
                         sources=[
                             'qutraj_run.pyf',
                             'qutraj_precision.f90',
                             'mt19937.f90',
                             'qutraj_general.f90',
                             'qutraj_hilbert.f90',
                             'qutraj_run.f90',
                             ],
                         libraries=libs,
                         **newblas)

    config.add_subpackage('examples')

    return config

if (__name__ == '__main__'):
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
    #setup(packages=['qutipf90mc'])

