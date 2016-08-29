.. _sec-chemistry:

Chemistry and Chemical Networks
===============================

The chemistry capabilities in DESPOTIC are located in the
``despotic.chemistry`` module.

.. _ssec-operations:

Operations on Chemical Networks
-------------------------------

The central object in a chemical calculation is a chemical network,
which consists of a set of chemical species and a set of reactions
that change the concentrations of each of them. Formally, a network is
defined by a set of :math:`N` species with abundances
:math:`\mathbf{x}` (defined as number density of each species per H
nucleus), and a function :math:`d\mathbf{x}/dt = \mathbf{f}(\mathbf{x},
\mathbf{p})` that gives the time rate of change of the abundances. The
vector :math:`\mathbf{p}` specifies a set of parameters on which the
reaction rates depend, and it may include quantities such as the
ambient radiation field (for photoreactions), the density, the
temperature, etc. The numerical implementation of chemical networks is
described in :ref:`ssec-chemnetworks`. Given any chemical network,
DESPOTIC is capable of two operations:

* Given an initial set of abundances at :math:`\mathbf{x}(t_0)` at
  time :math:`t_0`, compute the abundances at some later time
  :math:`t_1`. This is simply a matter of numerically integrating the
  ordinary differential equation :math:`d\mathbf{x}/dt =
  \mathbf{f}(\mathbf{x},\mathbf{p})` from :math:`t_0` to :math:`t_1`,
  where :math:`\mathbf{f}` is a known function. In DESPOTIC, this
  capability is implemented by the routine ``chemEvol`` in the
  ``despotic.chemistry.chemEvol`` module. The ``cloud`` class
  contains a wrapper around this routine, which allows it to be
  called to operate on a specific instance of ``cloud`` or
  ``zonedcloud``. 

* Find an equilibrium set of abundances
  :math:`\mathbf{x}_{\mathrm{eq}}` such that :math:`d\mathbf{x}/dt =
  \mathbf{f}(\mathbf{x}_{\mathrm{eq}}, \mathbf{p}) = 0`. Note that
  these is in general no guarantee that such an equilibrium exists, or
  that it is unique, and there are no general techniques for
  identifying such equilibria for arbitrary vector functions
  :math:`\mathbf{f}`. DESPOTIC handles this problem by explicitly
  integrating the ODE :math:`d\mathbf{x}/dt =
  \mathbf{f}(\mathbf{x},\mathbf{p})` until :math:`\mathbf{x}` reaches
  constant values (within some tolerance) or until a specified maximum
  time is reached. In DESPOTIC, this capability is implemented by
  the routine ``setChemEq`` in the ``despotic.chemistry.setChemEq``
  module. The ``cloud`` and ``zonedcloud`` classed contains a wrapper
  around this routine, which allows it to be called to operate on a
  specific instance of ``cloud`` or ``zonedcloud``.

.. _ssec-predefined-networks:

Predefined Chemical Networks
----------------------------

DESPOTIC ships with two predefined chemical networks, described below.

.. _sssec-NL99:

``NL99``
^^^^^^^^

The ``NL99`` network implements the reduced C-O network introduced by
`Nelson & Langer (1999, Astrophysical Journal, 524, 923; hereafter
NL99) <http://adsabs.harvard.edu/abs/1999ApJ...524..923N>`_. Readers
are refereed to that paper for a full description of the network and
the physical approximations on which it relies. To summarize briefly
here, the network is intended to capture the chemistry of carbon and
oxygen as it occurs at moderate densities and low temperatures in
:math:`\mathrm{H}_2` dominated clouds. It includes the species C,
:math:`\mathrm{C}^+`, CHx, CO, :math:`\mathrm{HCO}^+`,
:math:`\mathrm{H}^+_3`, :math:`\mathrm{He}^+`, O, OHx, M, and
:math:`\mathrm{M}^+`. Several of these are "super-species" that
agglomerate several distinct species with similar reaction rates and
pathways, including CHx (where :math:`x = 1-4`), OHx (where :math:`x =
1-2`), and M and :math:`\mathrm{M}^+` (which are stand-ins for metals
such as iron and nickel). The network involves two-body reactions
among these species, as well as photochemical reactions induces by UV
from the ISRF and reactions initiated by cosmic ray ionizations. In
addition to the initial abundances of the various species, the network
depends on the ISRF, the ionization rate, and the total abundances of
C and O nuclei.

In implementing the NL99 network in DESPOTIC there are a three design
choices to be made. First, photochemical and ionization reactions
depend on the UV radiation field strength and the ionization
rate. When performing computations on a cloud, DESPOTIC takes these
from the parameters ``chi`` and ``ionRate`` that are part of the radiation
class attached to the cloud.

Second, photochemical reactions depend on the amount of shielding
against the ISRF provided by dust, and, in the case of the reaction
:math:`\mathrm{CO}\rightarrow\mathrm{C}+\mathrm{O}`, line shielding by
CO and :math:`\mathrm{H}_2`. Following its usual approximation for
implementing such shielding in a one-zone model, DESPOTIC takes the
relevant column density to be :math:`N_{\mathrm{H}}/2`, where
:math:`N_\mathrm{H}` is the column density colDen of the cloud, so
that the typical amount of shielding is assumed to correspond to
half the area-averaged column density. For the dust shielding,
NL99 express the shielding in terms of the V-band extinction
:math:`A_V`; unless instructed otherwise, DESPOTIC computes this via

.. math::

   A_V = 0.4 \sigma_{\mathrm{PE}}(N_{\mathrm{H}}/2).

This ratio of V-band to 100 nm extinction is intermediate
between the values expected for Milky Way and SMC dust opacity curves,
as discussed in `Krumholz, Leroy, & McKee (2011, Astrophysical
Journal, 731, 25)
<http://adsabs.harvard.edu/abs/2011ApJ...731...25K>`_. However, the
user may override this choice. For line shielding, DESPOTIC computes
the :math:`\mathrm{H}_2` and CO column densities

.. math::

   N_{\mathrm{H}_2} = x_{\mathrm{H}_2} N_{\mathrm{H}}/2

   N_\mathrm{CO} = x_\mathrm{CO} N_\mathrm{H} /2,

which amounts to assuming that the CO and :math:`\mathrm{H}_2` are
uniformly distributed. Note that the NL99 network explicitly assumes
:math:`x_{\mathrm{H}_2} = 0.5`, as no reactions involving atomic H are
included -- see :ref:`sssec-NL99-GC` for a network that does include
hydrogen chemistry. These column densities are then used to find a
shielding factor by interpolating the tabulated values of `van Dishoeck
& Black (1988, Astrophysical Journal, 334, 771)
<http://adsabs.harvard.edu/abs/1988ApJ...334..771V>`_.

The third choice is how to map between the species included in the
chemistry network and the actual emitting species that are required
to compute line emission, cooling, etc. This is non-trivial both
because the chemical network includes super-species, because the
chemical network does not distinguish between ortho- and para-
sub-species while the rest of DESPOTIC does, and because the network
does not distinguish between different isotopomers of the same
species, while the rest of DESPOTIC does. This does not create
problems in mapping from cloud emitter abundances to the chemical
network, since the abundances can simply be summed, but it does create
a question about how to map output chemical abundances computed by
the network into the abundances of emitters that can be operated on by
the remainder of DESPOTIC. In order to handle this issue, DESPOTIC
makes the following assumptions:

1. OHx is divided evenly between OH and :math:`\mathrm{OH}_2`
2. The ratio of ortho- to para- for all species is the same as that of
   :math:`\mathrm{H}_2`
3. The abundances ratios of all isotopomers of a species remain fixed
   as reactions occur, so, for example, the final abundance ratio of
   :math:`\mathrm{C}^{18}\mathrm{O}` to 
   :math:`\mathrm{C}^{16}\mathrm{O}` as computed by the chemical
   network is always the same as the initial one. 


.. _sssec-NL99-GC:

``NL99_GC``
^^^^^^^^^^^

The ``NL99_GC`` network is an extension of the ``NL99`` network to
handle hydrogen chemistry as well as carbon and oxygen chemistry,
following the recipe described in `Glover & Clark (2012, MNRAS,
421, 116) <http://adsabs.harvard.edu/abs/2012MNRAS.421..116G>`_. This
network effectively uses ``NL99`` for carbon and oxygen (with some
updates to the rate coefficients) and the network of `Glover & Mac Low
(2007, ApJS, 169, 239)
<http://adsabs.harvard.edu/abs/2007ApJS..169..239G>`_
for the hydrogen chemistry. Self-shielding of hydrogen is handled via
the approximate analytic shielding function of `Draine & Bertoldi
(1996, ApJ, 468, 269)
<http://adsabs.harvard.edu/abs/1996ApJ...468..269D>`_. Effective
column densities for shielding are computed as in ``NL99``.
Calculations of hydrogen chemistry are assumed to leave the ratio of
ortho- to para- :math:`\mathrm{H}_2` unchanged from the initial value,
or set it to 0.25 if the initial value is undefined (e.g., because the
calculation started with a composition that was all H and no
:math:`\mathrm{H}_2`). All other assumptions are completely analogous
to ``NL99``.


.. _ssec-chemnetworks:

Defining New Chemical Networks: the ``chemNetwork`` Class
---------------------------------------------------------

DESPOTIC implements chemical networks through the abstract base class
``chemNetwork``, which is defined by module
``despotic.chemistry.chemNetwork`` -- see :ref:`sec-fulldoc` for full
details. This class defines the required elements that all chemistry
networks must contain; users who wish to implement their own chemistry
networks must derive them from this class, and must override the class
methods listed below. Users are encouraged to examine the two
:ref:`ssec-predefined-networks` for examples of how to derive a
chemical network class from ``chemNetwork``.

For any class ``cn`` derived from ``chemNetwork``, the user is
required to define the following non-callable traits:

* ``cn.specList``: a list of strings that describes the chemical
  species included in the network. The names in ``cn.specList`` can be
  arbitrary, and are not used for any purpose other than providing
  human-readable labels on outputs. However, it is often convenient to
  match the names to the names of emitters, as this makes it
  convenient to add the emitters back to the cloud later.
* ``cn.x``: a numpy array of rank 1, with each element specifying the
  abundance of a particular species in the network. The number of
  elements in the array must match the length of ``cn.specList``. As
  with all abundances in DESPOTIC, abundances must be specified
  relative to H nuclei, so that if, for example, ``x[3]`` is ``0.1``,
  this means that there is 1 particle of species 3 per 10 H nuclei.
* ``cn.cloud``: an instance of the ``cloud`` class to which the
  chemical network is attached. This can be ``None``, as chemical
  networks can be free-standing at not attached to specified instances
  of ``cloud``. However, much of the functionality of chemical
  networks is based around integration with the ``cloud`` class.

In addition, a class ``cn`` derived from ``chemNetwork`` must define
the following callable attributes:

* ``cn.__init__(self, cloud=None, info=None)``: this is the
  initialization method. It must accept two keyword arguments. The
  first, ``cloud``, is an instanced of the ``cloud`` class to which
  this chemical network will be attached. This routine should set
  ``cn.cloud`` equal to the input ``cloud`` instance, and it may also
  extract information from the input ``cloud`` instances in order to
  initialize ``cn.x`` or any other required quantities. The second
  keyword argument, ``info``, is a dict containing any additional
  initialization parameters that the chemical network can be or must
  be passed upon instantiation.
* ``cn.dxdt(self, xin, time)``: this is a method that computes the
  time derivative of the abundances. Given an input numpy array of
  abundances ``xin`` (which is the same shape as ``cn.x``) and a time
  ``time``, it must return a numpy array giving the time derivative of
  all abundances in units of :math:`\mathrm{s}^{-1}.`
* ``cn.applyAbundances(self, addEmitters=False)``: this is a method to
  take the abundances stored in the chemical network and use them to
  update the abundances of the corresponding emitters in the ``cloud``
  instance associated with this chemical network. This method is
  called at the end of every chemical calculation, and is responsible
  for copying information from the chemical network back into the
  cloud. The optional Boolean argument ``addEmitters``, if ``True``,
  specifies that the method should attempt to not only alter the
  abundances of any emitters associated with the cloud, it should also
  attempt to add new emitters that are included in the chemical
  network, and whose abundances are to be determined from it. It is up
  to the user whether or not to honor this request and implement this
  behavior. Failure to do so will not prevent the chemical network
  from operating, but should at least be warned so that other users
  are not surprised.

Finally, the following is an optional attribute of the derived class
``cn``:

* ``cn.abundances``: this is a property that returns the abundance
  information defined in ``cn.x`` as a object of class
  :ref:`sssec-abundanceDict`. This class is a utility class that
  provides an interface to the chemical abundances in a network
  that operates like a Python dict. The property abundances is defined
  in the base ``chemNetwork`` class, so it is available by inheritance
  regardless of whether the user defines ``cn.abundances``. However,
  the user may find it convenient to override
  ``chemNetwork.abundances`` to provide more information, e.g., to
  provide abundances of species that are not explicitly evolved in the
  network, and are instead derived via conservation relations. Both
  the :ref:`sssec-NL99` and :ref:`sssec-NL99-GC` classes do this.

Once a chemical network class that defines the above methods has been
defined, that class can be passed as an argument associated with the
``network`` keyword to the ``cloud.setChemEq`` and ``cloud.chemEvol``
methods (and their equivalents in ``zonedcloud``), and these methods
will automatically perform chemical calculations using the input
network.

In setting up chemical networks, it often convenient to make use of
the :ref:`ssec-chemhelpers` that are provided.


.. _ssec-chemhelpers:

Helper Modules for Chemical Networks
------------------------------------

The modules below, implemented in ``despotic.chemistry``, are intended
as helpers for defining and working with chemical networks.

.. _sssec-abundanceDict:

``abundanceDict``
^^^^^^^^^^^^^^^^^

The ``abundanceDict`` class provides a dict-like interface to numpy
arrays containing species abundances. The motivation for its existence
is that it is desirable to keep lists of chemical species abundances
in a numpy array for speed of computation, but for humans it is much
more convenient to be able to query and print chemical species by
name, with a dict-like interface. The ``abundanceDict`` class overlays
a dict-like structure on a numpy array, making it possible to combine
the speed of a numpy array with the convenience of a dict.

Usage is simple. To define an ``abundanceDict`` one simply provides a
set of chemical species names and a numpy array containing abundances
of those species::

  from despotic.chemistry import abundanceDict

  # specList = list of strings containing species names
  # x = numpy array containing species abundances
  abd = abundanceDict(specList, x)

Once created, an ``abundanceDict`` can be handled much like a dict, in
which the species names in ``specList`` are the keys::

  # Print abundance of the species CO
  print abd.abundances["CO"]

  # Set the C+ abundance to 1e-10
  abd["C+"] = 1e-10

These operations will alter the underlying numpy array appropriately;
the array may also be accessed directly, as
``abundanceDict.x``. However, the performance penalty for referring to
objects by name rather than by index in ``abundanceDict.x`` is
negligibly small in almost all cases, so it is recommended to use
keys, as this makes for much more human-readable code.

Mathematical operations on the abundances will pass through and be
applied to the underlying numpy array::

  # Double the abundance of every species
  abd = 2*abd

  # Add the abundances of two different abundanceDict's, abd1 and abd2
  abd3 = abd1 + abd2

Finally, ``abundanceDict`` instances differ from regular dict objects
in that, once they are created, their species lists are immutable;
species abundances can be changed, but species cannot be added or
removed.

For full details on ``abundanceDict``, see :ref:`sec-fulldoc`.


.. _sssec-reactions:

``reactions``
^^^^^^^^^^^^^

The ``reactions`` module provides a generic way to describe chemical
reactions and compute their rates. It contains one basic class and
several specialized derived classes to compute particualr types of
reactions. The goal of all the classes in ``reactions`` is to allow
users to describe the reactions in a chemical network in
human-readable forms, but then compute their rates using efficient
numpy operations at high speed.

.. _ssssec-reaction-matrix:

``reaction_matrix``
"""""""""""""""""""

The most general class in ``reactions`` is called
``reaction_matrix``. To instantiate it, one provides a list of species
and a list of reactions between them::

  from despotic.chemistry.reactions import reaction_matrix
  rm = reaction_matrix(specList, reactions)

Here ``specList`` is a list of strings giving the names of the
chemical species, while ``reactions`` is a list of dict's, one dict
per reaction, describing each reaction. Each dict in ``reactions``
contains two keys: ``reactions["spec"]`` gives the list of species
involved in the reaction, and ``reactions["stoich"]`` gives the
corresponding stoichiometric factors, with reactants on the left hand
side having negative factors (indicating that they are destroyed by
the reaction) and those on the right having positive factors
(indicating that they are created). Thus for example to include the
reactions

.. math::

  \mathrm{C} + \mathrm{O} & \rightarrow \mathrm{CO} \\
  \mathrm{H} + \mathrm{H} & \rightarrow \mathrm{H}_2

in a network, one could define the reactions as::

  reactions \
     = [ { "spec" : ["C", "O", "CO"], "stoich" : [-1, -1, 1] },
         { "spec" : ["H", "H2"], "stoich" : [-2, 1] } ]

Once a ``reaction_matrix`` has been defined, one can compute the rates
of change of all species using the method ``reaction_matrix.dxdt``::

  dxdt = rm.dxdt(x, n, ratecoef)

Here ``x`` is a numpy array giving the instantaneous abundances of all
species, ``n`` is the number density of H nuclei, and ``ratecoeff`` is
a numpy array giving the rate coefficient for all reactions. The
quantity returned is the instantaneous rate of change of all
abundances ``x``. The density-dependence of the reaction rates implied
by the stoichiometric factors in the reaction list is computed
automatically.


``cr_reactions``
""""""""""""""""

The ``cr_reactions`` class is a specialized class derived from
:ref:`ssssec-reaction-matrix` that is intended to handle cosmic
ray-induced reactions -- those where the rate is proportional to
the cosmic ray ionization rate.

Instantiation of a ``cr_reactions`` object takes the same two
arguments as :ref:`ssssec-reaction-matrix`, but ``reactions``, the
list of reactions to be passed, is altered in that each reaction takes
an additional key, ``rate``, that gives the proportionality between
the reaction rate and the cosmic ray ionization rate. Thus for example
if we wished to include the reactions

.. math::

   \mathrm{CR} + \mathrm{H} & \rightarrow \mathrm{H}^+ +
   \mathrm{e}^- \\
   \mathrm{CR} + \mathrm{H}_2 & \rightarrow \mathrm{H}_2^+ +
   \mathrm{e}^-

in a network, with the former occurring at a rate per particle equal to
the rate of primary cosmic ray ionizations, and the other at a rate
per particle that is twice the rate of primary ionizations, we would
define::

  from despotic.chemistry.reactions import cr_reactions
  specList = ["H", "H2", "H+", "H2+", "e-"]
  reactions \
     = [ { "spec" : ["H", "H+", "e-"], "stoich" : [-1, 1, 1],
           "rate" : 1.0 },
	 { "spec" : ["H2", "H2+", "e-"], "stoich" : [-1, 1, 1],
           "rate" : 2.0 } ]
  cr = cr_reactions(specList, reactions)

Once a ``cr_reactions`` object is instantiated, it can be used to
compute the rates of change all species using the
``cr_reactions.dxdt`` routine::

  dxdt = cr.dxdt(x, n, ionrate)

Here ``x`` is a numpy array giving the current abundances, ``n`` is
the number density of H nuclei, and ``ionrate`` is the cosmic ray
primary ionization rate.

Note that, in most cases, the variable ``n`` will not be used, because
it is needed only if there is more than one reactant on the left hand
side of a reaction. It is provided to enable the case where a cosmic
ray ionization is followed immediately by another reaction, and the
network combines the two steps for speed. For example, the
:ref:`sssec-NL99` network combines the two reactions

.. math::

   \mathrm{CR} + \mathrm{H}_2 & \rightarrow \mathrm{H}_2^+ +
   \mathrm{e}^- \\
   \mathrm{H}_2^+ + \mathrm{H}_2 & \rightarrow \mathrm{H}_3^+ +
   \mathrm{H} + \mathrm{e}^-

into the single super-reaction

.. math::

   \mathrm{CR} + 2\mathrm{H}_2 \rightarrow \mathrm{H}_3^+ +
   \mathrm{H} + \mathrm{e}^-

and the rate of this super-reaction does depend on density.

``photoreactions``
""""""""""""""""""

The ``photoreactions`` class is a specialized class derived from
:ref:`ssssec-reaction-matrix` that is intended to handle reactions
that are driven by FUV photons, and thus have rates proportional to
the interstellar radiation field (ISRF) strength, modified by dust and
gas shielding, and possibly also by line shielding.

Instantiation of ``photoreactions`` takes the same two arguments as
:ref:`ssssec-reaction-matrix`: a list of species, and a list of
reactions, each element of which is a dict giving the reactants and
stoichiometric factors. The dict also contains three additional keys
that are specific to photoreactions: ``rate`` gives the reaction rate
per reactant per second in an ISRF equal to the unshielded Solar
neighborhood value (formally, with strength characterized by
:math:`\chi = 1` -- see :ref:`tab-cloudfiles`). The key ``av_fac``
describes the optical depth per unit visual extinction provided by
dust for the photons driving the reaction. That is, reaction rates
will be reduced by a factor :math:`\exp(-\mathrm{av\_fac}\times
A_V)`. Finally, the key ``shield_fac``, which is optional, can be set
equal to a callable that describes the rate by which the reaction rate
is reduced due to line shielding. This function can take any number of
argument -- see below. Thus the final reaction rate per reactant will
be equal to

.. math::

   \mathrm{reaction\, rate} = \chi \times \mathrm{rate} \times
   \mathrm{shield\_fac} \times \exp(-\mathrm{av\_fac}\times A_V)

As an example, suppose we wished to include the reactions

.. math::

   \gamma + \mathrm{C} & \rightarrow \mathrm{C}^+ + \mathrm{e}^- \\
   \gamma + \mathrm{CO} & \rightarrow \mathrm{C} + \mathrm{O} \\

The first reaction occurs at a rate per C atom :math:`5.1\times
10^{-10}\,\mathrm{s}^{-1}` in unshielded space, and the optical depth
to the photons producing the reaction is :math:`3.0\times A_V`. There
is no significant self-shielding. The second reaction occurs at a rate
per CO molecule :math:`1.7\times 10^{-10}\,\mathrm{s}^{-1}` in
unshielded space, with a dust optical depth equal to :math:`1.7\times
A_V`, and with an additional function describing line shielding called
``fShield_CO``. We could create a ``photoeactions`` object containing
this reaction by doing::

  from despotic.chemistry.reactions import photoreactions
  specList = ["C", "C+", "e-", "O", "CO"]
  reactions \
     = [ { "spec" : ["C", "C+", "e-"], "stoich" : [-1, 1, 1],
           "rate" : 5.1e-10, "av_fac" : 3.0 },
         { "spec" : ["CO", "C", "O"], "stoich" : [-1, 1, 1],
           "rate" : 1.7e-10, "av_fac" : 1.7,
	   "shield_fac" : fShield_CO } ]
  phr = photoreactions(specList, reactions)

Once a photoreactions object exists, the rates of change of all
species due to the included photoreactions can be computed using the
``photoreactions.dxdt`` routine. Usage is as follows::

  dxdt = phr.dxdt(x, n, chi, AV,
                  shield_args=[[f1_arg1, f1_arg2, ...],
		               [f2_arg1, f2_arg2, ...],
			       ... ],
                  shield_kw_args = [ { "f1_kw1" : val1,
		                       "f1_kw2" : val2, ... },
				     { "f2_kw1" : val1,
				       "f2_kw2" : val2, ... },
				     ... ])
				     
Here ``x`` is a numpy array giving the current abundances, ``n`` is
the number density of H nuclei, ``chi`` is the unshielded ISRF
strength, and ``AV`` is the dust visual extinction. The optional
arguments ``shield_args`` and ``shield_kw_args`` are used to pass
positional arguments and keyword arguments, respectively, to the
shielding functions. Each one is a list. The first entry in the list
for ``shield_args`` is a list containing the positional arguments to
be passed to the first shielding function provided in the
``reactions`` dict passed to the constructor. The second entry in
``shield_args`` is a list containing positional arguments to be passed
to the second shielding function, and so forth. Similarly, the
first entry in the ``shield_kw_args`` list is a dict containing the
keywords and their values to be passed to the first shielding
function, etc. The length of the ``shield_args`` and
``shield_kw_args`` must be equal to the number of shielding functions
provided in ``reactions``, but elements of the lists can be set to
``None`` if a particular shielding function does not take positional
or keyword arguments. Moreover, both ``shield_args`` and
``shield_kw_args`` are optional. If they are omitted, then all
shielding functions are invoked with no positional or keyword
arguments, respectively.

Thus in the above example, suppose that the function ``fShield_CO``
returns the factor by which the reaction rate is reduced by CO line
shielding. It takes two arguments: ``NH2`` and ``NCO``, giving the
column densities of :math:`\mathrm{H}_2` and :math:`\mathrm{CO}`,
respectively. It also accepts a keyword argument ``order`` that
specifies the order of interpolation to be used in computing the
shielding function from a table. In this case one could compute the
rates of change of the abundances via::

  dxdt = phr.dxdt(x, n, chi, AV,
                  shield_args = [[NH2, NCO]],
		  shield_kw_args = [ {"order" : 2} ])

This would call ``fShield_CO`` as::

  fShield_CO(NH2, NCO, order=2)

