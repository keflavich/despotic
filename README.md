DESPOTIC
--------

DESPOTIC is a radiative transfer code hosted at [https://sites.google.com/a/ucsc.edu/krumholz/codes/despotic](https://sites.google.com/a/ucsc.edu/krumholz/codes/despotic).

Here's the README from that page:

The DESPOTIC Code

DESPOTIC is a library to Derive the Energetics and SPectra of Optically Thick
Interstellar Clouds. Its capabilities include calculating spectral line
luminosities from clouds of specified physical properties and compositions,
predicting line profiles, calculating rates of heating and cooling due to a
wide range of processes, determining clouds' equilibrium gas and dust
temperatures, and calculating time-dependent thermal evolution of clouds.

DESPOTIC is implemented as a Python package, and is released under the GPL. The
physical model, equations, solved, and some tests and example applications are
described in 
[Krumholz (2014)](http://adsabs.harvard.edu/abs/2014MNRAS.437.1662K),
and an extensive [User's Guide](https://docs.google.com/a/ucsc.edu/viewer?a=v&pid=sites&srcid=dWNzYy5lZHV8a3J1bWhvbHp8Z3g6YzBkMDk3OWI3NzA3OTU0) is available, and is included with the download.


Downloading DESPOTIC
====================

DESPOTIC is distributed in two ways. The full source distribution, including
example Python programs using it and sample cloud descriptor files, is hosted
at [google code](https://code.google.com/p/despotic/). It may be
obtained from the SVN repository by typing:

`svn checkout http://despotic.googlecode.com/svn/trunk/ despotic-read-only`

This will create a directory called despotic containing the distribution. The
User's Guide is located in the despotic/doc directory.

Despotic is also available through the [python package index](https://pypi.python.org/pypi/DESPOTIC). Just type

`pip install despotic`

The tarball may also be downloaded from PyPi and installed using the
standard Python distutils method.
