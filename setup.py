#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
        name="qutipf90mc",
        version="0.2",
        packages=find_packages(),
        author="Arne L. Grimsmo",
        license="MIT",
        install_requires=["numpy", "scipy"],
)
