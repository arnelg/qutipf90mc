#!/usr/bin/env python
"""
qutipf90mc doc
"""

DOCLINES = __doc__.split('\n')

from os.path import join

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
            maintainer = "Arne L. Grimsmo",
            maintainer_email = "arne.grimsmo@gmail.com",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "",
            download_url = "https://github.com/arnelg/qutipf90mc",
            license = 'GPL3',
            #classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
            platforms = ["Linux", "Mac OS-X"],
            configuration=configuration )
    finally:
        #del sys.path[0]
        #os.chdir(old_path)
        pass
    return


if __name__ == '__main__':
    setup_package()
