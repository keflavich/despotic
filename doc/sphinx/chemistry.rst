.. highlight:: chemistry

.. _sec-chemistry:

Chemistry and Chemical Networks
===============================

The chemistry capabilities in DESPOTIC are located in the
``despotic.chemistry`` module.

.. _sec-operations:

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

Available Chemical Networks
---------------------------

DESPOTIC ships with two predefined chemical networks, described below.

.. _sssec-NL99:

``NL99``
~~~~~~~~

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
a question about how to map output chemical abun- dances computed by
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
~~~~~~~~~~~

.. _ssec-chemnetworks:

Defining New Chemical Networks: the ``chemNetwork`` Class
---------------------------------------------------------

Helper Modules for Chemical Networks
------------------------------------

``abundanceDict``
~~~~~~~~~~~~~~~~~

``reactions``
~~~~~~~~~~~~~
