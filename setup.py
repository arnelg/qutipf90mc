#!/usr/bin/env python
"""
qutipf90mc doc
"""

DOCLINES = __doc__.split('\n')

from distutils.core import Command
import os
from glob import glob
from os.path import join
from os.path import splitext, basename, join as pjoin
from unittest import TextTestRunner, TestLoader

#------ clean command for removing .pyc, .o, .mod and .so files --------#

class CleanCommand(Command):
    user_options = [("all", "a", "All")]

    def initialize_options(self):
        self._clean_me_pyc = []
        self._clean_me_o = []
        self._clean_me_so = []
        self._clean_me_mod = []
        self._clean_me_other = []
        self.all = None
        for root, dirs, files in os.walk('.'):
            for f in files:
                if f.endswith('.pyc'):
                    self._clean_me_pyc.append(pjoin(root, f))
                if f.endswith('.o'):
                    self._clean_me_o.append(pjoin(root, f))
                if f.endswith('.so'):
                    self._clean_me_so.append(pjoin(root, f))
                if f.endswith('.mod'):
                    self._clean_me_mod.append(pjoin(root, f))
                if f=='qutraj_run-f2pywrappers2.f90':
                    self._clean_me_other.append(pjoin(root, f))
                if f=='qutraj_runmodule.c':
                    self._clean_me_other.append(pjoin(root, f))

    def finalize_options(self):
        pass

    def run(self):
        pyc_rm=0; o_rm=0; so_rm=0; mod_rm=0
        for clean_me in self._clean_me_pyc:
            try:
                os.unlink(clean_me)
            except:
                pyc_rm+=1
        for clean_me in self._clean_me_o:
            try:
                os.unlink(clean_me)
            except:
                o_rm+=1
        for clean_me in self._clean_me_so:
            try:
                os.unlink(clean_me)
            except:
                so_rm+=1
        for clean_me in self._clean_me_mod:
            try:
                os.unlink(clean_me)
            except:
                mod_rm+=1
        if pyc_rm>0:
            print("Could not remove "+str(pyc_rm)+" pyc files.")
        else:
            print("Removed all .pyc files.")
        if o_rm>0:
            print("Could not remove "+str(o_rm)+" .o files.")
        else:
            print("Removed all .o files.")
        if so_rm>0:
            print("Could not remove "+str(so_rm)+" .so files.")
        else:
            print("Removed all .so files.")
        if mod_rm>0:
            print("Could not remove "+str(mod_rm)+" .mod files.")
        else:
            print("Removed all .mod files.")
        for clean_me in self._clean_me_other:
            try:
                os.unlink(clean_me)
                print "Removed", clean_me
            except:
                print "Could not remote", clean_me

#--------- configuration -------------#

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration(None, parent_package, top_path)
    config.add_subpackage('qutipf90mc')
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    return config


def setup_package():
    from numpy.distutils.core import setup
    try:
        setup(
            name = 'qutipf90mc',
            version = '0.1',
            maintainer = "Arne L. Grimsmo",
            maintainer_email = "arne.grimsmo@gmail.com",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "",
            download_url = "https://github.com/arnelg/qutipf90mc",
            license = 'GPL3',
            #classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
            platforms = ["Linux", "Mac OS-X"],
            cmdclass = { 'clean': CleanCommand},
            configuration=configuration )
    finally:
        #del sys.path[0]
        #os.chdir(old_path)
        pass
    return


if __name__ == '__main__':
    setup_package()
