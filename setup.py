#!/usr/bin/env python

from setuptools import setup, find_packages
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
      entry_points={"console_scripts": ["genosv = genosv.app.main:main"]},
      install_requires=["pysam>=0.10", "numpy", "pyfaidx", "tqdm"], 
     )
