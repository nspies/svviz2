#!/usr/bin/env python

from setuptools import setup, find_packages
import sys

if sys.version_info < (3, 3):
    raise RuntimeError("Python version >= 3.3 required.")


setup(name="genosv",
      version="0.1",
      description="genosv",
      author="Noah Spies",
      packages=find_packages(),
      entry_points={"console_scripts": ["genosv = genosv.app.main:main"]},
      install_requires=["pysam>=0.10", "numpy", "pyfaidx", "tqdm"], 
     )
