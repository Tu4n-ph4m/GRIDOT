#!/usr/bin/env python

from setuptools import setup, find_packages
import gridot

setup(
    name="gridot_packages",
    version=gridnet.__version__,
    description="GRN inferences with Granger Graph and Optimal transport",
    # url="https://github.com/alexw16/gridnet",
    license="GPLv3",
    packages=find_packages(),
  
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "torch",
        "scikit-learn",
        "scanpy"
    ],
)
