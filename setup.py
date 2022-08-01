#!/usr/bin/env python
import os
from setuptools import find_packages, setup

# Package meta-data.
NAME = 'cliff'
DESCRIPTION = 'A Python API for calculating epistasis and ruggness of mutation dataset.'
URL = 'https://github.com/cutecutecat/cliff'
EMAIL = 'starkind1997@gmail.com'
AUTHOR = 'Henery Chen'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '1.0'
REQUIRED = [
    "pandas>=1.3.4",
    "networkx>=2.6.3",
    "matplotlib>=3.4.3",
    "numpy>=1.20.3",
    "joblib>=1.1.0",
    "tqdm>=4.62.2",
    "click>=8.0.3"
]

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

REQUIRED = []
with open(os.path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    REQUIRED = f.read().splitlines()
print(REQUIRED)

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('test',)),
    install_requires=REQUIRED,
    include_package_data=True,
    license='MIT',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='protein ruggness epistasis high-order genetics genotype-phenotype-maps',
    entry_points={
        "console_scripts": [
            "cliff = cliff.client:cli",
        ],
    },
)
