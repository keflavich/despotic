#!/usr/bin/env python

import os
import os.path as osp
import shutil
import sys
from distutils.core import setup, Command
import subprocess
import numpy
from Cython.Build import cythonize
if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup

with open('README.md') as file:
    long_description = file.read()

with open('CHANGES') as file:
    long_description += file.read()

with open('REQUIREMENTS') as file:
    requirements = file.readlines()

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py


from despotic import __version__ as version
tagname = "despotic_%s" % (version)

download_url = "https://bitbucket.org/krumholz/despotic/get/dfc8dc4b047d.zip"

class PyTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        errno = subprocess.call([sys.executable, 
                                 osp.join('tests', 'run_tests.py')])
        raise SystemExit(errno)

class RunExamples(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        errno = subprocess.call([sys.executable, 
                                 osp.join('examples', 'run_tests.py')])
        raise SystemExit(errno)


import os 
# despotic/cloudfiles is a link back to ./cloudfiles
if not osp.exists(osp.join('despotic', 'cloudfiles')):
    os.symlink(osp.join('..', 'cloudfiles'),
               osp.join('despotic', 'cloudfiles'))

setup(name='DESPOTIC',
      version=version,
      description='a Python / numPy / sciPy package to perform calculations related to line emission and thermal behavior in cold interstellar clouds.',
      long_description=long_description,
      author=['Mark Krumholz'],
      author_email=['mark.krumholz@gmail.com'], 
      url='https://bitbucket.org/krumholz/despotic/',
      download_url=download_url,
      packages=['despotic','despotic.chemistry'],
      package_dir={'despotic':'despotic'}, 
      package_data={'despotic':[osp.join('cloudfiles','*.desp')]},
      requires=requirements,
      cmdclass={'build_py': build_py, 'test': PyTest, 'examples': RunExamples},
      classifiers=[
                   "Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
                  ],
      ext_modules = cythonize(osp.join('despotic', 'collPartner_helper.pyx')),
      include_dirs = [numpy.get_include()]
     )

