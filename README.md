# qutipf90mc

A "wave-function monte carlo" solver written in Fortran 90/95, with a
python interface trough f2py. The program is intended to be used with
the [QuTiP](https://github.com/qutip/qutip)  python package.

## Features:

* Usage (almost, see missing features) identical to QuTiP's `mcsolve`
* Uses sparse (compressed row format) matrices for operators
* Uses zvode to integrate in time
* Time evolution algorithm from QuTiP to find correct times for jumps.
* Automatic parallelization via Python's multiprocessing module.

Missing features:
* Does not accept list as "ntraj" argument.
* Only solves prolbems without explicit time-dependence.


Dependencies:

* QuTiP v.3.1.0 or higher
* A fortran compiler and the BLAS library. BLAS comes with many fortran compilers, such as gfortran. Unoffifical gfortran binaries can be found here https://gcc.gnu.org/wiki/GFortranBinaries

## Platforms:

Has been tested on OSX with anaconda python 3.5, gfortran 6.3 and qutip
v4.2.0.

## Installation:

1. Download code with
```shell
git clone https://github.com/arnelg/qutipf90mc.git
```

2. Enter directory and install
```shell
python setup.py install
```

Or, if you prefer to install locally:
```shell
python setup.py build_ext --inplace
```

## Testing and usage:

Test the installation by leaving the directory, starting python and entering
```python
import qutipf90mc
```

To compare the speed of `mcsolve_f90` vs. `mcsolve` for a decaying
system with Hilbert space dimension `dim`, and `ntraj` trajectories, run on a single CPU.
```python
qutipf90mc.compare.run(dim,ntraj)
```

For more info
```python
help(qutipf90mc.mcsolve_f90)
```

## License

You are free to use this software, with or without modification, provided that the conditions listed in the LICENSE.txt file are satisfied.

## Contributors

A. L. Grismmo, P. D. Nation, and J. R. Johansson
