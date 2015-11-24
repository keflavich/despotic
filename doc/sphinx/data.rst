.. highlight:: data

.. _sec-data:

Atomic and Molecular Data
=========================

DESPOTIC requires atomic and molecular data to work. This section
describes how it handles these data, both on disk and in its internal
workings.

The Local Database
------------------

DESPOTIC uses atomic and molecular data in the format specified by the
`Leiden Atomic and Molecular Database
<http://home.strw.leidenuniv.nl/~moldata/>`_. The user can manually
supply the required data files, but the more common use case is to
access the data directly from LAMDA. When emitter data is required,
DESPOTIC will attempt to guess the name of the required data file and
download it automatically -- see :ref:`ssec-emitters`. When DESPOTIC
downloads a file from `LAMDA
<http://home.strw.leidenuniv.nl/~moldata/>`_, it caches a local copy
for future use. The next time the same emitter is used, unless
DESPOTIC is given an explicit URL from which the file should be
fetched, it will use the local copy instead of re-downloading the file
from LAMDA. (However, see :ref:`ssec-database-updates`.)

The location of the database is up to the user, and is specified
through the environment variable ``$DESPOTIC_HOME``. If this
environment variable is set, LAMDA files will be places in the
directory ``$DESPOTIC_HOME/LAMDA``, and that is the default location
that will be searched when a file is needed. If the environment
variable ``$DESPOTIC_HOME`` is not set, DESPOTIC looks for files in a
subdirectory LAMDA of the current working directory, and caches files
in that directory if they are downloaded. It is recommended that users
set a ``$DESPOTIC_HOME`` environment variable when working with
DESPOTIC, so as to avoid downloading and caching multiple copies of
LAMDA for different projects in different directories.

.. _ssec-database-updates:

Keeping the Atomic and Molecular Data Up to Date
------------------------------------------------

The data in LAMDA are updated regularly as new calculations or
laboratory experiments are published. Some of these updates add new
species, but some also provide improved data on species that are
already in the database. DESPOTIC attempts to ensure that its locally
cached data are up to date by putting an expiration date on them. By
default, if DESPOTIC discovers that a given data file is more than six
months old, it will re-download that file from LAMDA. This behavior
can be overridden by manually specifying a file name, either in the
cloud file (see :ref:`sec-cloudfiles`) or when invoking
the ``cloud.addEmitter`` or ``emitter.__init__`` methods. Users
can also force updates of the local database more frequently using the
``refreshLamda`` function -- see :ref:`sec-full`.

.. _ssec-database-internal:

DESPOTIC's Internal Model for Atomic and Molecular Data
-------------------------------------------------------

When it is running, DESPOTIC maintains a list of emitting species for
which data have been read within the ``emitter`` module (see
:ref:`sec-full`). Whenever a new emitter is created, either for an
existing cloud, for a new cloud being created, or as a free-standing
object of the emitter class, DESPOTIC checks the emitter name against
the central list. If the name is found in the list, DESPOTIC will
simply use the stored data for that object rather than re-reading the
file containing the data. This is done as an efficiency measure, and
also to ensure consistency between emitters of the same species
associated with different clouds. However, this model has some
important consequences of which the user should be aware.

1. Since data on level structure, collision rates, etc. (everything
   stored in the ``emitterData`` class -- see :ref:`sec-full`) is
   shared between all emitters of the same name, and any alterations
   made to the data for one emitter will affect all others of the same
   name.
2. It is not possible to have two emitters of the same name but with
   different data. Should a user desire to achieve this for some
   reason (e.g., to compare results computed using an older LAMDA file
   and a newer one), the way to achieve this is to give the two
   emitters different names, such as ``co_ver1`` and ``co_ver2``.
3. Maintenance of a central emitter list affects how deepcopy and
   pickling operations operate on emitters. See
   :ref:`sec-full` for details.
