#!/usr/bin/env python

from setuptools import setup, find_packages
import gridot,schema  # or import gridot.__version__ if it's part of the module

setup(
    name="gridot_packages",
    version=gridot.__version__,  # Ensure __version__ is defined in gridot
    description="GRN inferences with Granger Graph and Optimal transport",
    license="GPLv3",
    packages=find_packages(include=["gridot", "schema"]),  # Include relevant packages
    include_package_data=True,  # Optional if you have non-Python files to include
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas",
        "scipy",
        "seaborn",
        "cvxopt",
        "tables",
        "tqdm",
        "scikit-learn",
        "umap-learn",
        "legacy-api-wrap",
        "setuptools_scm",
        "packaging",
        "sinfo",
        "scanpy"
    ],
)
