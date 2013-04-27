#!/usr/bin/env python

import os
import shutil
import sys
from distutils.core import setup, Command
import subprocess
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

download_url = "https://despotic.googlecode.com/svn/trunk/dist/DESPOTIC-%s.tar.gz" % version

class PyTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        errno = subprocess.call([sys.executable, 'tests/run_tests.py'])
        raise SystemExit(errno)

import os 
# despotic/cloudfiles is a link back to ./cloudfiles
if not os.path.exists('despotic/cloudfiles'):
    os.symlink('../cloudfiles','despotic/cloudfiles')

setup(name='DESPOTIC',
      version=version,
      description='a Python / numPy / sciPy package to perform calculations related to line emission and thermal behavior in cold interstellar clouds.',
      long_description=long_description,
      author=['Mark Krumholz'],
      author_email=['mkrumhol@ucsc.edu'], 
      url='https://code.google.com/p/despotic/',
      download_url=download_url,
      packages=['despotic','despotic.chemistry'],
      package_dir={'despotic':'despotic'}, 
      package_data={'despotic':['cloudfiles/*.desp']},
      requires=requirements,
      cmdclass={'build_py': build_py, 'test': PyTest},
      classifiers=[
                   "Development Status :: 3 - Alpha",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
                  ],
      
     )

