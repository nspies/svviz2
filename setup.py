#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext # Cython should be installed via pysam
import pysam
import numpy
import sys

if sys.version_info < (3, 3):
    raise RuntimeError("Python version >= 3.3 required.")

def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(name="genosv",
      version=get_version(open('genosv/__init__.py').read()),
      description="genosv",
      author="Noah Spies",
      packages=find_packages(),
      ext_modules = [Extension("genosv.remap._mapq",
                              sources=["genosv/remap/_mapq.pyx"],
                              include_dirs=pysam.get_include()+[numpy.get_include()],
                              define_macros=pysam.get_defines())],
      cmdclass={'build_ext': build_ext},
      entry_points={"console_scripts": ["genosv = genosv.app.main:main"]},
      install_requires=["cython", "pysam>=0.10", "numpy", "pyfaidx", "tqdm", "pandas", "numpy"], 
     )
