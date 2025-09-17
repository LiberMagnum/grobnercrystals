#!/usr/bin/env sage
from setuptools import setup # type: ignore
from GCS import __version__

setup(
    name='GCS',
    version=__version__,
    description='A SageMath package for working with Gröbner crystal structures',

    url='https://github.com/LiberMagnum/GCS',
    author='Abigail Price',
    author_email='price29@illinois.edu',

    py_modules=['GCS'],
    install_requires=[
        'sage-package',
        'numpy',
    ],

    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)