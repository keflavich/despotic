.. highlight:: rest

Functional Guide to DESPOTIC Capabilities
=========================================

This section covers the most common tasks for which DESPOTIC can be
used. They are not intended to provide a comprehensive overview of the
library's capabilities, and users who wish to understand every
available option should refer to :ref:`sec-full`. The routines used in
this section are all described in full detail there. For all of these
operations, the user should first have imported the basic DESPOTIC
class ``cloud`` by doing::

  from despotic import cloud

In the examples below we will also assume that ``matplotlib`` and
``numpy`` have both been imported, via::

  import matplotlib.pyplot as plt
  import numpy as np

Unit Conventions
----------------

In this section, and in general when using DESPOTIC, there are two
important conventions to keep in mind.

1. All quantities are in CGS units unless otherwise specified. The
   main exceptions are quantites that are normalized to Solar or Solar
   neighborhood values (e.g. metallicity and interstellar radiation
   field strength) and quantities where the conventional unit is
   non-CGS (e.g. integrated brightness temperatures are expressed in
   the usual units of K km/s).
2. All rates are expressed per H nucleus, rather than per unit mass or
   per unit volume. This includes chemical abundances. Thus for
   example a heating rate of :math:`\Gamma=10^{-26}` should be
   understood as :math:`10^{-26}` erg/s/H nucleus. An abundance
   :math:`x_{\mathrm{He}}=0.1` should be understood as 1 He atom per
   10 H nuclei.

.. _ssec-line-emission:

Line Emission
-------------

The most basic task in DESPOTIC is computing the line emission
emerging from a cloud of specified physical properties. The first step
in any computation of line emission is to create a cloud object with
the desired properties. This is often most easily done by creating a
DESPOTIC cloud file (see :ref:`sec-cloudfiles`), but the user can also
create a cloud with the desired properties manually. The properties
that are important for line emission are the gas volume density,
column density, temperature, velocity dispersion, and chemical
composition; these are stored in the cloud class and the composition
class within it. For example, the following code snippet::

  mycloud = cloud()
  mycloud.nH = 1.0e2
  mycloud.colDen = 1.0e22
  mycloud.sigmaNT = 2.0e5
  mycloud.Tg = 10.0
  mycloud.comp.xoH2 = 0.1
  mycloud.comp.xpH2 = 0.4
  mycloud.comp.xHe = 0.1

creates a cloud with all its parameters set to default values, a
volume density of H nuclei :math:`n_{\mathrm{H}} =
10^2\,\mathrm{cm}^{-3}`, a column density of H nuclei
:math:`N_{\mathrm{H}} = 10^{22}\,\mathrm{cm}^{-2}`, a non-thermal
velocity dispersion of :math:`\sigma_{\mathrm{NT}} = 2.0 \times
10^5\,\mathrm{cm}\,\mathrm{s}^{-1}`, a gas temperature of :math:`T_g =
10\,\mathrm{K}`, and a composition that is 0.1
ortho-:math:`\mathrm{H}_2` molecules per H nucleus, 0.4
para-:math:`\mathrm{H}_2` molecules, and 0.1 He atoms per H nucleus.

The next step is to specify the emitting species whose line emission
is to be computed. As with the physical properties of the cloud, this
is often most easily done by creating a cloud file. However, it can
also be done manually by using the cloud.addEmitter method, which
allows users to add emitting species to clouds. The following code
snippet adds CO as an emitting species, at an abundance of
:math:`10^{-4}` CO molecules per H nucleus::

   mycloud.addEmitter("CO", 1.0e-4)

The first argument is the name of the emitting species, and the second
is the abundance. The requisite molecular data file will be read from
disk if it is available, or automatically downloaded from the Leiden
Atomic and Molecular Database if it is not (see :ref:`sec-data`).

Once an emitter has been added, only a single call is required to
calculate the luminosity of its lines::

  lines = mycloud.lineLum("CO")

The argument is just the name of the species whose line emission
should be computed. Note that it must match the name of an existing
emitter, and that emitter names are case-sensitive. The value returned
by this procedure, which is stored in the variable ``lines``, is a
``list`` of ``dict``, with each ``dict`` describing the properties of
a single line. Lines are ordered by frequency, from lowest to
highest. Each ``dict`` contains the following keys-value pairs

* ``freq`` is the line frequency in Hz
* ``intIntensity`` is the frequency-integrated intensity of the line
  after subtracting off the CMB contribution, in
  :math:`\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{sr}^{-1}`
* ``intTB`` is the velocity-integrated brightness temperature (again
  subtracting off the CMB) in
  :math:`\mathrm{K}\,\mathrm{km}\,\mathrm{s}^{-1}`
* ``lumPerH`` is the rate of energy emission in the line per H nucleus
  in the cloud, in :math:`\mathrm{erg}\,\mathrm{s}^{-1}`. 

This is a partial list of what the ``dict`` contains; see
:ref:`sec-full` for a complete listing.

Once the data been obtained, the user can do what he or she wishes
with them. For example, to plot velocity-integrated brightness
temperature versus line frequency, the user might do::

  freq = [l["freq"] for l in lines]
  TB = [l["intTB"] for l in lines]
  plt.plot(freq, TB, "o")

Heating and Cooling Rates
-------------------------

To use DESPOTIC's capability to calculate heating and cooling rates,
in addition to the quantities specified for a calculation of line
emission one must also add the quantities describing the dust and the
radiation field. As before, this is most easily accomplished by
creating a DESPOTIC cloud file (see :ref:`sec-cloudfiles`), but the
data can also be input manually. The code snippet below does so::

  mycloud.dust.alphaGD   = 3.2e-34    # Dust-gas coupling coefficient
  mycloud.dust.sigma10   = 2.0e-25    # Cross section for 10K thermal radiation
  mycloud.dust.sigmaPE   = 1.0e-21    # Cross section for photoelectric heating
  mycloud.dust.sigmaISRF = 3.0e-22    # Cross section to the ISRF
  mycloud.dust.beta      = 2.0        # Dust spectral index
  mycloud.dust.Zd        = 1.0        # Abundance relative to Milky Way
  mycloud.Td             = 10.0       # Dust temperature
  mycloud.rad.TCMB       = 2.73       # CMB temperature
  mycloud.rad.TradDust   = 0.0        # IR radiation field seen by the dust
  mycloud.rad.ionRate    = 2.0e-17    # Primary ionization rate
  mycloud.rad.chi        = 1.0        # ISRF normalized to Solar neighborhood

These quantities specify the dust-gas coupling constant, the dust
cross section to 10 K thermal radiation, the dust cross section to the
8 - 13.6 eV photons the dominate photoelectric heating, the dust cross
section to the broader interstellar radiation field responsible for
heating the dust, the dust spectral index, the dust abundance relative
to the Milky Way value, the dust temperature, the cosmic microwave
background temperature, the infrared radiation field that heats the
dust, the primary ionization rate due to cosmic rays and x-rays, and
the ISRF strength normalized to the Solar neighborhood value. All of
the numerical values shown in the code snippet above are in fact the
defaults, and so none of the above commands are strictly
necessary. However, it is generally wise to set quantities explicitly
rather than relying on default values.

Once these data have been input, one may compute all the heating and
cooling terms that DESPOTIC includes using the ``cloud.dEdt`` routine::

  rates = mycloud.dEdt()

This call returns a dict which contains the instantaneous rates of
heating and cooling. The entries in the dict are: ``GammaPE``, the gas
photoelectric heating rate, ``GammaCR``, the gas heating rate due to
cosmic ray and X-ray ionization, ``GammaGrav``, the gas heating rate due
to gravitational compression, ``GammaDustISRF``, the dust heating rate
due to the ISRF, ``GammaDustCMB``, the dust heating rate due to the CMB,
``GammaDustIR``, the dust heating rate due to the IR field,
``GammaDustLine``, the dust heating rate due to absorption of line
photons, ``PsiGD``, the gas-dust energy exchange rate (positive means
gas heating, dust cooling), ``LambdaDust``, the dust cooling rate via
thermal emission, and ``LambdaLine``, the gas cooling rate via line
emission. This last quantity is itself a dict, with one entry per
emitting species and the dictionary keys corresponding to the emitter
names. Thus in the above example, one could see the cooling rate via
CO emission by doing::

  print rates["LambdaLine"]["CO"]

.. _ssec-temp-eq:

Temperature Equilibria
----------------------

Computing the equilibrium temperature requires exactly the same
quantities as computing the heating and cooling rates; indeed, the
process of computing the equilibrium temperature simply amounts to
searching for values of :math:`T_g` and :math:`T_d` such that the sum
of the heating and cooling rates returned by ``cloud.dEdt`` are zero. One
may perform this calculation using the ``cloud.setTempEq`` method::

  mycloud.setTempEq()

This routine iterates to find the equilibrium gas and dust
temperatures, and returns True if the iteration converges. After this
call, the computed dust and gas temperatures may simply be read off::

  print mycloud.Td, mycloud.Tg

The ``cloud.setTempEq`` routine determines the dust and gas
temperatures simultaneously. However, there are many situations where
it is preferable to solve for only one of these two, while leaving the
other fixed. This may be accomplished by the calls::

  mycloud.setDustTempEq()
  mycloud.setGasTempEq()

These routines, respectively, set ``mycloud.Td`` while leaving
``mycloud.Tg`` fixed, or vice-versa. Solving for one temperature at a
time is often faster, and if dust-gas coupling is known to be
negligible will produce nearly identical results as solving for the
two together.

.. _ssec-temp-evol:

Time-Dependent Temperature Evolution
------------------------------------

To perform computations of time-dependent temperature evolution,
DESPOTIC provides the method ``cloud.tempEvol``. In its most basic
form, this routine simply accepts an argument specifying the amount of
time for which the cloud is to be integrated, and returning the
temperature as a function of time during this evolution (note that
executing this command may take a few minutes, depending on your
processor)::

  mycloud.Tg = 50.0         # Start the cloud out of equilibrium
  tFinal = 20 * 3.16e10     # 20 kyr
  Tg, t = mycloud.tempEvol(tFinal)

The two values returned are arrays, the second of which gives a series
of 100 equally-spaced times between 0 and ``tFinal``, and the first of
which gives the temperatures at those times. The number of output
times, the spacing between them, and their exact values may all be
controlled by optional arguments -- see :ref:`fulldoc` for details. At
the end of this evolution, the cloud temperature ``mycloud.Tg``
will be changed to its value at the end of 20 kyr of evolution, and
the dust temperature ``mycloud.Tg`` will be set to its thermal
equilibrium value at that cloud temperature.

If one wishes to examine the intermediate states in more detail, one
may also request that the full state of the cloud be saved at every
time::

  clouds, t = mycloud.tempEvol(tFinal, fullOutput=True)

The ``fullOutput`` optional argument, if ``True``, causes the routine
to return a full copy of the state of the cloud at each output time,
instead of just the gas temperature ``Tg``. In this case, ``clouds``
is a sequence of 100 ``cloud`` objects, and one may interrogate their
states (e.g. calculating their line emission) using the usual
routines.

.. _ssec-chem-eq:

Chemical Equilibria
-------------------

DESPOTIC can also compute the chemical state of clouds from a chemical
network. Full details on chemical networks are given in
:ref:`sec-chemistry`, but for this example we will use a simple network
that DESPOTIC ships with, that of `Nelson & Langer (1999, ApJ,
524, 923) <http://adsabs.harvard.edu/abs/1999ApJ...524..923N>`_. This
network computes the chemistry of carbon and oxygen in a region where
the hydrogen is fully molecular. For more details see
:ref:`sssec-NL99`.

To perform computations with this network, one must first import the
class that defines it::

  from despotic.chemistry import NL99

One can set the equilibrium abundances of a cloud to the equilibrium
values determined by the network via the command::

  mycloud.setChemEq(network=NL99)

The argument ``network`` specifies that the calculation should use the
``NL99`` class. This call sets the abundances of all the emitters that
are included in the network to their equilibrium values. In this case,
the network includes CO, and thus it sets the CO abundance to a new
value::

  print mycloud.emitters["CO"].abundance

One can also see the abundances of all the species included in the
network, including those that do not correspond to emitters in the
cloud, by printing the chemical network property ``abundances``::

  print mycloud.chemnetwork.abundances

Once the chemical network is associated with the cloud, subsequent
calls to ``setChemEq`` need not include the ``network``
keyword. DESPOTIC assumes that all subsequent chemical calculations
are to be performed with the same chemical network unless it is
explicitly told otherwise via a call to ``setChemEq`` or
``chemEvol`` (see :ref:`ssec-chem-time`) that specifies a different
chemical network.

Simultaneous Chemical and Thermal Equilibria
--------------------------------------------

The ``setChemEq`` routine (see :ref:`ssec-chem-eq`) called with no
extra arguments leaves the gas temperature fixed. However, it is also
possible to compute a simultaneous equilibrium for the temperature and
the thermal state. To do so, we first import a chemical network to
be used, in this case the Nelson & Langer (1999) network (see
:ref:`ssec-chem-eq`)::

  from despotic.chemistry import NL99

We then call ``cloud.setChemEq`` with an optional keyword
``evolveTemp`` ::

  mycloud.setChemEq(network=NL99, evolveTemp=``iterate``)

The ``network`` keyword specifies that the computation should use the
NL99 network, while ``evolveTemp`` specifies how to handle the
simultaneous thermal and chemical equilibrium calculation. The options
available are

* ``fixed``: gas temperature is held fixed
* ``iterate``: calculation iterates between computing chemical and
  gas thermal equilibria, i.e., chemical equilibrium is computed at
  fixed temperature, equilibrium gas temperature (see
  :ref:`ssec-temp-eq`) is computed for fixed abundances, and the
  process is repeated until the temperature and abundances converge;
  dust temperature is held fixed
* ``iterateDust``: same as ``iterate``, except the dust temperature is
  iterated as well
* ``gasEq``: gas temperature is always set to its instantaneous
  equilibrium value as the chemical state is evolved toward
  equilibrium; dust temperature is held fixed
* ``fullEq``: same as ``gasEq``, except that both gas and dust
  temperatures are set to their instantaneous equilibrium values
* ``evol``: chemical state and gas temperature are evolved in time
  together, while dust temperature is always set to its instantaneous
  equilibrium value; evolution stops once gas temperature and
  abundances stop changing significantly

Note that, while in general the different evolution methods will
converge to the same answer, there is no guarantee that they will do
in systems where multiple equilibria exist.

.. _ssec-chem-time:

Time-Dependent Chemical Evolution
---------------------------------

DESPOTIC can also calculate time-dependent chemical evolution. This is
accomplished through the method cloud.chemEvol. At with
``cloud.tempEvol`` (see :ref:`ssec-temp-evol`), this routine accepts
an argument specifying the amount of time for which the cloud is to be
integrated, and returning the chemical abundances as a function of
time during this evolution::

  mycloud.rad.ionRate = 2.0e-16 # Raise the ionization rate a lot
  tFinal = 0.5 * 3.16e13 # 0.5 Myr
  abd, t = mycloud.chemEvol(tFinal, network=NL99)

Note that the ``network=NL99`` option may be omitted if one has
previously assigned that network to the cloud (for example by
executing the examples in :ref:`ssec-chem-eq`).

The output quantity abd here is an object of class ``abundanceDict``,
which is a specialized dict for handling chemical abundances -- see
:ref:`ssec-abundanceDict`. One can examine the abundances of specific
species just by giving their chemical names. For example, to see
the time-dependent evolution of the abundances of CO, C, and
:math:`\mathrm{C}^+`, one could do::

  plt.plot(t, abd["CO"])
  plt.plot(t, abd["C"])
  plt.plot(t, abd["C+"])

As with ``setChemEq``, this routine modifies the abundances of
emitters in the cloud to the values they achieve at the end of the
evolution, so to see the final CO abundance one could do::

  print mycloud.emitters["CO"].abundance

Multi-Zone Clouds
-----------------

While most DESPOTIC functionality is provided through the ``cloud``
class, which represents a single cloud, it is sometimes useful to have
a cloud that contains zones of different optical depths. This
functionality is provided through the ``zonedcloud`` class. A
``zonedcloud`` is just a collection of ``cloud`` objects that are
characterized by having different column densities (and optionally
volume densities), and on which all the operations listed above can be
performed in a batch fashion.

One can create a ``zonedcloud`` in much the same way as a ``cloud``,
but reading from an input file::

  from despotic import zonedcloud
  zc = zonedcloud(fileName="cloudfiles/MilkyWayGMC.desp")

A ``zonedcloud`` is characterized by column densities for each of its
zones, which can be accessed through the ``colDen`` property::

  print zc.colDen

The column densities of all zones, and the number of zones, can be
controlled when the ``zonedcloud`` is created using the keywords
``nZone`` and ``colDen``; see :ref:`ssec-full` for the full list of
keywords.

Once a ``zonedcloud`` exists, all of the functions described above in
this section are available for it, and will be applied zone by
zone. For example, one can do::

  zc.setTempEq()
  print zc.Tg

to set and then print the temperature in each zone. Commands the
report observable quantities or abundances will return
appropriately-weighted sums over the entire cloud. For example::

  zc.lineLum('co')[0]

returns a dict describing the :math:`J=1\rightarrow 0` line of CO. The
quantities ``intTB`` and ``intIntensity`` that are part of the dict
and contain the velocity-integrated brightness temperature and
frequency-integrated intensity, respectively (see
:ref:`ssec-line-emission`), are sums over all zones, while ones
like ``Tex`` (the excitation temperature) that do not make sense to
sum are returned as an array giving zone-by-zone values.


Computing Line Profiles
-----------------------

Line profile computation operates somewhat differently then the
previous examples, because it is provided through a stand-alone
procedure rather than through the cloud class. This procedure is
called lineProfLTE, and may be imported directly from the DESPOTIC
package. The routine also requires emitter data stored in an
``emitterData`` object. The first step in a line profile calculation is
therefore to import these two objects into the python environment::

  from despotic import lineProfLTE
  from despotic import emitterData

The second step is to read in the emitter data. The interface to read
emitter data is essentially identical to the one used to add an
emitter to a cloud. One simply declares an ``emitterData`` object,
giving the name of the emitter as an argument::

  csData = emitterData(’CS’) # Reads emitter data for the CS molecule

Alternately, emitter data may be obtained from a ``cloud``, since clouds
store emitter data for all their emitters. Using the examples from the
previous sections::

  coData = mycloud.emitters["CO"].data

copies the emitter data for CO to the variable ``coData``.

The third step is to specify the radius of the cloud, and the profiles
of any quantities within the cloud that are to change with radius,
including density, temperature, radial velocity, and non-thermal
velocity dispersion. Each of these can be constant, but the most
interesting applications are when one or more of them are not, in
which case they must be defined by functions. These function each take
a single argument, the radius in units where the outer radius of the
cloud is unity, and return a single floating point value, giving the
quantity in question in CGS units. For example, to compute line
profiles through a cloud of spatially-varying temperature and infall
velocity, one might define the functions::

  R = 0.02 * 3.09e18 # 0.2 pc
  def TProf(r):
      return 8.0 + 12.0*np.exp(-r**2/(2.0*0.5**2))
  def vProf(r):
      return -4.0e4*r

The first function sets a temperature that varies from 20 K in the
center of close to 8 K at the outer edge, and the second defines a
velocity that varies from 0 in the center to :math:`-0.4` km
:math:`\mathrm{s}^{-1}` (where negative indicates infall) at the outer
edge. Similar functions can be defined by density and non-thermal
velocity dispersion if the user so desires. Alternately, the user can
simply define them as constants::

  ncs = 0.1       # CS density 0.1 cm^-3
  sigmaNT = 2.0e4 # Non-thermal velocity dispersion 0.2 km s^-1

The final step is to use the ``lineProfLTE`` routine to compute the
brightness temperature versus velocity::

  TB, v = lineProfLTE(cs, 2, 1, R, ncs, TProf, vProf, sigmaNT).

Here the first argument is the emitter data, the second and third are
the upper and lower quantum states between which the line is to be
computed (ordered by energy, with ground state = 0), followed by the
cloud radius, the volume density, the temperature, the velocity, and
the non-thermal velocity dispersion. Each of these quantities can be
either a float or a callable function of one variable, as in the
example above. If it is a float, that quantity is taken to be
constant, independent of radius. This routine returns two arrays, the
first of which is the brightness temperature and the second of which
is the velocity at which that brightness temperature is computed,
relative to line center. These can be examined in any of the usual
numpy ways, for example by plotting them::

  plt.plot(v, TB)

By default the velocity is sampled at 100 values. The routine attempts
to guess a reasonable range of velocities based on the input values of
radial velocity and velocity dispersion, but these defaults may be
overridden by the optional argument ``vLim``, which is a sequence of
two values giving the lower and upper limits on the velocity::

  TB, v = lineProfLTE(cs, 2, 1, R, ncs, TProf, vProf, sigmaNT,
                      vLim=[-2e5,2e5]).

A variety of other optional arguments can be used to control the
velocities at which the brightness temperature is computed. It is also
possible to compute line profiles at positions offset from the
geometric center of the cloud, using the optional argument offset --
see :ref:`ssec-full`.

Escape Probability Geometries
-----------------------------

DESPOTIC supports three possible geometries that can be used when
computing escape probabilities, and which are controlled by the
``escapeProbGeom`` optional argument. This argument is accepted by all
DESPOTIC functions that use the escape probability formalism,
including all those involving computation of line emission. This
optional argument, if included, must be set equal to one of the three
strings ``sphere`` (the default), ``slab``, or ``LVG``. These choices
correspond to spherical geometry, slab geometry, and the large
velocity gradient approximation, respectively.
