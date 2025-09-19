#!/usr/bin/env sage
from setuptools import setup # type: ignore

setup(
    name='grobnercrystals',
    version='1.1.2',
    description='A SageMath package for working with Gröbner crystal structures',

    url='https://github.com/LiberMagnum/grobnercrystals',
    author='Abigail Price',
    author_email='price29@illinois.edu',

    py_modules=['grobnercrystals'],
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