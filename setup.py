# -*- coding: utf-8 -*-

from setuptools import setup
from os import path
import sys
#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
from numpy.distutils.core import Extension as Extension
from numpy.distutils.core import setup as setup

import numpy
from Cython.Build import cythonize

sys.path.insert(0, "piszkespipe")
from version import __version__

# Load requirements
requirements = None
with open('requirements.txt') as file:
    requirements = file.read().splitlines()

# If Python3: Add "README.md" to setup.
# Useful for PyPI. Irrelevant for users using Python2.
try:
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
except:
    long_description = ' '

# Command-line tools
entry_points = {'console_scripts': [
    'piszkespipe = piszkespipe:piszkespipe_from_commandline'
]}

extensions = Extension("piszkespipe.Marsh",
                         sources=["piszkespipe/OptExtract/Marsh.c"],
                         include_dirs=[numpy.get_include(),'/opt/local/include','/usr/local/include'],
                         library_dirs=['/opt/local/lib'],
                         libraries=['m', 'gsl', 'gslcblas'])

"""
flib = Extension('piszkespipe.CCF',
              sources=['piszkespipe/CCF/CCF.f90','piszkespipe/CCF/CCF.pyf'],
              swig_opts=['-c', '--fcompiler=gnu95','--f77flags="-ffixed-line-length-none"'],
              define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])
"""

flib = Extension("piszkespipe.CCF",
             sources=["piszkespipe/CCF/CCF_pyc.pyx", "piszkespipe/CCF/CCF.c"],
             include_dirs=[numpy.get_include(),'/opt/local/include'],
             library_dirs=['/opt/local/lib'],
             libraries=['m'],
             define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])

ext_modules = [extensions,flib]

desc='A pipeline for reducing echelle spectra obtained in Piszkesteto.'

setup(name='piszkespipe',
      version=__version__,
      description=desc,
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Attila Bodi',
      author_email='astrobatty@gmail.com',
      url='https://github.com/astrobatty/piszkespipe/',
      classifiers=['Intended Audience :: Science/Research', 'Topic :: Scientific/Engineering :: Astronomy', 'Development Status :: 5 - Production/Stable', 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)', 'Operating System :: OS Independent', 'Programming Language :: Python'],
      packages=['piszkespipe','piszkespipe.Correlation','piszkespipe.GLOBALutils'],
      package_data={"piszkespipe": ["wavcals/*.dat", "data/xc_masks/*", "Correlation/*.dat"]},
      include_package_data=True,
      install_requires=requirements,
      entry_points=entry_points,
      ext_modules = cythonize(ext_modules, compiler_directives={'language_level': 3})
     )
