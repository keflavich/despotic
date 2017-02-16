.. highlight:: rest

Installing DESPOTIC
===================

DESPOTIC is distributed in two ways. 

Installing from git
-------------------

The full source distribution, including example Python programs using it and sample cloud descriptor files, is available from `bitbucket <https://bitbucket.org/krumholz/despotic/>`_, and can be obtained via git by doing::

  git clone git@bitbucket.org:krumholz/despotic.git

The package can then be installed by doing::

  python setup.py install

in the ``despotic`` directory. You may need to preface this with
``sudo`` if you want to install in a way that is accessible for all
users on the system.

Installing from pip
-------------------

DESPOTIC is also available through the python package index. If you
have pip installed, you can just type::

  pip install despotic

As with the setup via ``git``, you may need to preface this with
``sudo`` to install globally.

Setting Up the Environment
--------------------------

DESPOTIC creates a local cache of atomic and molecular data files. You
can specify the location of this cache by setting the
``$DESPOTIC_HOME`` environment variable; the data cache will be
created in ``$DESPOTIC_HOME/LAMDA`` If you are using a ``bash``-like
shell, the syntax to set the location of ``$DESPOTIC_HOME`` is::

   export DESPOTIC_HOME = /path/to/despotic

while for a ``csh``-like shell, it is::

   setenv DESPOTIC_HOME /path/to/despotic

If ``$DESPOTIC_HOME`` is not set, then DESPOTIC will attempt to create
the cache in the directory from which it is run.

Requirements and Dependencies
-----------------------------

DESPOTIC requires

* `scipy >= 0.11.0`
* `cython >= 0.20.x`
* `matplotlib >= 1.3.x`

.. _ssec-winds-installation:
  
Installing the ``despotic.winds`` Module
----------------------------------------

The ``despotic.winds`` module (see :ref:`sec-winds`) is written in
C++, and makes use of the `GNU scientific library
<http://www.gnu.org/software/gsl/>`_. Users who wish to make use of
this package will need to compile the C++ portions of the
software. This process should be automatic for users who have an
installed C++ compiler and standard build tools (``make``), and who
have the `GSL <http://www.gnu.org/software/gsl/>`_ headers and library
installed in their ``CXX_INCLUDE_PATH`` and ``LIBRARY_PATH``,
respectively.

If this does not work, users should be able to build the required C++
library by changing into the ``despotic/winds`` directory and doing
``make``. Users can edit the ``Makefile`` in ``despotic/winds`` as
necessary to tailor the build for their environment.

The ``winds`` module is not required for the operation of the rest of
the code, so users who are not planning to make use of this capability
can safely ignore it.
