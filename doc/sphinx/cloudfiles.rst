.. highlight:: cloudfiles

.. _sec-cloudfiles:

Cloud Files
===========

File Structure
--------------

DESPOTIC cloud files contain descriptions of clouds that can be read
by the ``cloud`` or ``zonedcloud`` classes, using either the class
constructor or the read method; see :ref:`sec-full` for
details. This section contains a description of the format for these
files. It is recommended but not required that cloud files have names
that end in the extension ``.desp``.

Each line of a cloud file must be blank, contain a comment starting
with the character ``#``, or contain a key-value pair formatted as::

  key = value

The line may also contain comments after ``value``, again beginning
with ``#``. Any content after ``#`` is treated as a comment and is
ignored. DESPOTIC keys are case-insensitive, and whitespace around
keys and values are ignored. For the full list of keys, see
:ref:`tab-cloudfiles`. All quantities must
be in CGS units. Key-value pairs may be placed in any order, with the
exception of the key ``H2opr``, which provides a means of specify the
ratio of ortho-to-para-:math:`\mathrm{H}_2`, instead of directly
setting the ortho-and para-:math:`\mathrm{H}_2` abundances. If this
key is specified, it must precede the key ``xH2``, which gives the
total :math:`\mathrm{H}_2` abundance including both ortho- and para-
species. Not all keys are required to be present. If left unspecified,
most quantities default to a fiducial Milky Way value (if a reasonable
one exists, e.g., for the gas-dust coupling constant and
ISRF strength) or to 0 (if it does not, e.g., for densities and
chemical abundances).

.. _tab-cloudfiles:

.. table:: Cloud file keys and their meanings

   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Key        | Units                                               | Description                                        |
   +============+=====================================================+====================================================+
   |                                                                                                                       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Physical Properties                                                                                                   |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | nH         | :math:`\mathrm{cm}^{-3}`                            | Volume density of H nuclei                         |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | colDen     | :math:`\mathrm{cm}^{-2}`                            | Column density of H nuclei, averaged over area     |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | sigmaNT    | cm :math:`\mathrm{s}^{-1}`                          | Non-thermal velocity dispersion                    |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Tg         | K                                                   | Gas temperature                                    |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Td         | K                                                   | Dust temperature                                   |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   |                                                                                                                       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Dust Properties                                                                                                       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | alphaGD    | erg :math:`\mathrm{cm}^3` :math:`\mathrm{K}^{-3/2}` | Gas-dust collisional coupling coefficient          |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | sigmaD10   | :math:`\mathrm{cm}^2` :math:`\mathrm{H}^{-1}`       | Dust cross section per H to 10 K thermal radiation |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | sigmaDPE   | :math:`\mathrm{cm}^2` :math:`\mathrm{H}^{-1}`       | Dust cross section per H to 8-13.6 eV photons      |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | sigmaDISRF | :math:`\mathrm{cm}^2` :math:`\mathrm{H}^{-1}`       | Dust cross section per H averaged over ISRF        |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | betaDust   | Dimensionless                                       | Dust IR spectral index                             |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Zdust      | Dimensionless                                       | Dust abundance normalized to Milky Way value       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   |                                                                                                                       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Radiation Field Properties                                                                                            |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | TCMB       | K                                                   | CMB temperature                                    |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | TradDust   | K                                                   | Temperature of the dust-trapped IR radiation field |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | ionRate    | :math:`\mathrm{s}^{-1}`                             | Primary cosmic ray / x-ray ionization rate         |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | chi        | Dimensionless                                       | ISRF strength, normalized to Solar neighborhood    |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   |                                                                                                                       |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | Chemical composition                                                                                                  |
   +------------+-----------------------------------------------------+----------------------------------------------------+
   | emitter    |                                                     | See :ref:`ssec-emitters`                           |
   +------------+-----------------------------------------------------+----------------------------------------------------+

.. _ssec-emitters:

Emitters
--------

The ``emitter`` key is more complex than most, and requires special
mention. Lines describing emitters follow the format::

  emitter = name abundance [noextrap] [energySkip] [file:FILE] [url:URL]

Here the brackets indicate optional items, and the optional items may
appear in any order, but must be after the two mandatory ones.

The first mandatory item, ``name``, gives name of the emitting
molecule or atom. Note that molecule and atom names are case
sensitive, in the sense that DESPOTIC will not assume that ``co`` and
``CO`` describe the same species. Any string is acceptable for
``name``, but if the file or URL containing the data for that species
is not explicitly specified, the name is used to guess the
corresponding file name in the Leiden Atomic and Molecular Database
(LAMDA) -- see :ref:`sec-data`. It is therefore generally advisable
to name a species following LAMDA convention, which is that molecules
are specified by their chemical formula, with a number specifying the
atomic weight preceding the if the species is not the most common
isotope. Thus LAMDA refers to :math:`^{28}\mathrm{Si}^{16}\mathrm{O}`
(the molecule composed of the most common isotopes) as ``sio``,
:math:`^{29}\mathrm{Si}^{16}\mathrm{O}` as ``29sio``, and
:math:`^{12}\mathrm{C}^{18}\mathrm{O}` as ``c18o``. The automatic
search for files in LAMDA also includes common variants of
the file name used in LAMDA. The actual file name from
which DESPOTIC reads data for a given emitter is stored in
the emitterData class -- see :ref:`sec-full`. 

The second mandatory item, ``abundance``, gives the abundance of that
species relative to H nuclei. For example, an abundance of 0.1 would
indicate that there is 1 of that species per 10 H nuclei.

The optional items ``noextrap`` and ``energySkip`` change how DESPOTIC
performs computations with that species. If ``noextrap`` is set,
DESPOTIC will raise an error if any attempt is made to calculate a
collision rate coefficient between that species and one of the bulk
components (H, He, etc.) that is outside the range tabulated in the
data file. If not, DESPOTIC will instead handle temperatures outside
the tabulated range by substituting the closest temperature inside the
tabulated range. Note that this behavior can be altered within a
DESPOTIC program by using the ``extrap`` property of the
``emitterData`` class -- see :ref:`sec-fulldoc`. 

The optional item ``energySkip`` specifies that a species should be
ignored when computing heating and cooling rates via the
``cloud.dEdt`` method, and by extension whenever thermal equilibria or
thermal evolution are computed for that cloud. However, line emission
from that species can still be computed using the ``cloud.lineLum``
method. This option is therefore useful for including species for
which the line emission is an interesting observable, but which are
irrelevant to the thermal balance and thus can be omitted when
calculating cloud thermal properties in order to save computational
time.

Finally, the optional items ``file:FILE`` and ``url:URL`` specify
locations of atomic and molecular data files, either on the local file
system or on the web. This capability is useful in part because some
LAMDA files do not follow the usual naming convention, or because for
some species LAMDA provides more than one version of the data for that
species (e.g., two versions of the data file for atomic C exist, one
with only the low-lying IR levels, and another including the
higher-energy UV levels). File specifications must be of the form
``file:FILE`` with ``FILE`` replaced by a file name, which can include
both absolute and relative paths. If no path or a relative path is
given, DESPOTIC searches for the file first in the current directory,
and then in the directory ``$DESPOTIC_HOME/LAMDA``, where
``$DESPOTIC_HOME`` is an environment variable. If it is not specified,
DESPOTIC just looks for a directory called LAMDA relative to the
current directory.

The ``url:URL`` option can be used to specify the location of a file
on the web, usually somewhere on the LAMDA website. It must be
specified as ``url:URL``, where ``URL`` is replaced by an absolute or
relative URL. If an absolute URL is given, DESPOTIC attempts to
download the file from that location. If a relative URL is given,
DESPOTIC attempts to download the file from at
``http://$DESPOTIC_LAMDAURL/datafiles/URL``, where
``$DESPOTIC_LAMDAURL`` is an environment variable. If this environment
variable is not specified, DESPOTIC searches for the file at
``http://home.strw.leidenuniv.nl/~moldata/URL``.
