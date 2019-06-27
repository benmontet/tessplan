#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import sys
from setuptools import setup

sys.path.insert(0, "tessplan")


long_description = \
    """
tessplan is a tool to schedule observations of TESS objects of interest.
It enables the user to determine visible transits of a TOI from any observatory
at any future time, and can determine in a given time window what planets will 
be transiting. It retrieves information from ExoFOP-TESS, so a user can use that
information (such as a TFOP Working Group prioritization) to select possible
planet candidates to observe. 
"""


setup(
    name='tessplan',
    version='0.0.1rc3',
    license='MIT',
    author='Benjamin Montet',
    author_email='bmontet@uchicago.edu',
    packages=[
        'tessplan',
        ],
    include_package_data=True,
    url='http://github.com/benmontet/tessplan',
    description='A tool for scheduling follow-up of transiting TESS Objects of Interest',
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={'': ['README.md', 'LICENSE']},
    install_requires=[
        'numpy', 'astroquery', 'astropy',
        'astroplan', 'pytz',
        'requests', 'tqdm', 
        'beautifulsoup4'],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.0',
        ],
    )
