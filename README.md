DESPOTIC
--------

DESPOTIC is a radiative transfer code hosted at http://www.ucolick.org/~krumholz/codes/despotic/

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
described in <a href="http://arxiv.org/abs/1304.2404">Krumholz (2013)</a>, and an
extensive <a href="doc/UsersGuide.pdf">User's Guide</a> is available, and is
included with the download.


Downloading DESPOTIC
====================

DESPOTIC is distributed in two ways. The full source distribution, including
example Python programs using it and sample cloud descriptor files, is hosted
at <a href="https://code.google.com/p/despotic/">google code</a>. It may be
obtained from the SVN repository by typing:

`svn checkout http://despotic.googlecode.com/svn/trunk/ despotic`

This will create a directory called despotic containing the distribution. The
User's Guide is located in the despotic/doc directory.

Alternately, users who just want the DESPOTIC code itself and not the
associated examples and documentation can obtain it from the <a
href="https://pypi.python.org/">Python Package Index</a>. For those with <a
href="http://www.pip-installer.org">pip</a> or <a
href="http://peak.telecommunity.com/DevCenter/EasyInstall">EasyInstall</a>
installed, installation is as simple as typing

`pip install despotic`

The tarball may also be downloaded from the <a
href="https://pypi.python.org/pypi/DESPOTIC/">DESPOTIC page</a> in the PYPI and
installed using the standard Python distutils method.
