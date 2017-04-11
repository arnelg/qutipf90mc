#!/usr/bin/env python
"""
qutipf90mc doc
"""


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
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
            name='qutipf90mc',
            version='0.2',
            maintainer="Arne L. Grimsmo",
            maintainer_email="arne.grimsmo@gmail.com",
            configuration=configuration)
    finally:
        pass
    return


if __name__ == '__main__':
    setup_package()
