.. _sec-winds:

Winds
=====

The wind modelling capability in DESPOTIC is located in the
``despotic.winds`` module. The unofficial name for this portion of the
code is ROCKETSTAR -- as named by the primary author's 5-year old.

The physical model used in the ``despotic.winds`` module is described
in Krumholz et al. (2017), and users are encouraged to read that paper
to understand the physical basis of the calculations. Example scripts
using this code can be found in the repository
`https://bitbucket.org/krumholz/despotic_winds/
<https://bitbucket.org/krumholz/despotic_winds/>`_ associated with
this paper.

.. _ssec-wind-pwind:

The ``pwind`` Class
-------------------

The front end to the ``despotic.winds`` is the ``pwind`` class. This
class requires that the user specify the generalized Eddington ratio
:math:`\Gamma` and the Mach number :math:`\mathcal{M}` for the
wind. By default this will create a spherical, ideal wind, which is
fully defined by these two parameters. In addition, the the user can
specify a number of other proprties:

* The rate of cloud expansion, specified by the keyword
  ``expansion``. Valid values are ``area``, ``intermediate``, and
  ``solid angle``.

* The gravitational potential, specified by the keyword
  ``potential``. Valid values are ``point`` and ``isothermal``.

* The geometry of the wind, specified by the keyword
  ``geometry``. Valid values are ``sphere``, ``cone``, and
  ``cone sheath``. For ``cone`` the user must also specify the cone
  tilt (via the keyword ``phi``) and the cone opening angle (via the
  keyword ``theta``). For ``cone sheath`` the user must also specify
  the inner opening angle (via the keyword ``theta_in``).

* The driving mechanism, specified by the keyword ``driver``. Valid
  values are ``ideal``, ``radiation``, and ``hot``. For ``radiation``,
  the user must specify the optical depth at the mean surface density
  (via the keyword ``tau0``). For ``hot``, the user must specify the
  dimensionless hot gas velocity (via the keyword ``uh``).

The are a range of other keywords the affect the behavior of ``pwind``
objects. See :ref:`sssec-full-pwind` for a full listing.

An example is::

  import numpy as np
  from despotic.winds import pwind

  Gamma = 0.2
  Mach = 50.
  tau0 = 100.
  phi = np.pi/4.0
  theta = np.pi/2.0

  pw = pwind(Gamma, Mach, driver='radiation', potential='isothermal',
             expansion='solid angle', geometry='cone',
	     tau0=tau0, theta=theta, phi=phi)

This creates a ``pwind`` object that represents a radiatively-driven
wind in an isothermal potential, with clouds maintaining constant
solid angle. The wind is confined to a cone that is tipped by
:math:`45^\circ` relative to the vertical, with a :math:`90^\circ`
opening angle. The wind is characterized by a generalized Eddington
ratio :math:`\Gamma = 0.2` and a Mach number :math:`\mathcal{M} =
50`, and the optical depth at the mean surface density of the launch
region is :math:`\tau_0 = 100`.

Caculations Using ``pwind``
---------------------------

The ``pwind`` class defines a series of methods that can be used to
compute the observable properties of the specified wind. There are
four basic types of observables that can be computed:

* ``pwind.tau``: this method computes absorption optical depths. The
  user must specify the dimensionless velocity or velocities ``u`` at
  which to compute the absorption, as well as the dimensionless
  transition strength :math:`t_X/t_w` for the wind. This can be
  specified either directly, via the keyword ``tXtw``, or computed
  from an input oscillator strength ``Omega``, wavelength
  ``wl``, abundance ``abd``, and wind mass removal timescale
  ``tw``. The keyword ``correlated`` specifies whether the wind should
  be treated as correlated or uncorrelated. The keyword ``u_trans``
  specifies that the transition in question is a multiplet, with
  the individual transitions occurring at dimensionless velocities
  given by ``u_trans``.

* ``pwind.Xi``: this method returns the optically thin emission line
  shape function :math:`\Xi`, as defined in the Krumholz et al. (2017)
  paper. The user must specify the dimensionless velocity or
  velocities ``u``.

* ``pwind.temp_LTE``: this returns the brightness or antenna
  temperature for a species in LTE. The user must provide the
  dimensionless velocity or velocities ``u`` and the wind kinetic
  temperature ``T``. In addition, the user must provide the dimensionless
  transition strength :math:`t_X/t_w` for the wind. This can be
  specified either directly, via the keyword ``tXtw``, or computed
  on the fly in one of two ways. First, the user can specify
  ``Omega``, ``wl``, ``abd`` and ``tw``, exactly as for
  ``pwind.tau``. Second, the user can provide a DESPOTIC ``emitter``
  object (see :ref:`sssec-full-emitter`) and a wind removal timescale
  ``tw``. Finally, the keyword ``correlated`` specifies whether the
  wind should be treated as correlated or uncorrelated.

* ``pwind.intTA_LTE``: this computes the velocity-integrated antenna
  temperature for an emitting species in LTE. The user must provide
  the velocity scale ``v0`` for the wind launching region and the wind
  kinetic temperature ``T``. All other parameters are as for
  ``pwind.temp_LTE``.

All of these routines accept the keywords ``varpi`` and ``varpi_t``
which specify the dimensionless axial and transverse position of the
line of sight (:math:`\varpi_a` and :math:`\varpi_t` in the
terminology of Krumholz et al. 2017). In addition, all routines except
``pwind.Xi`` accept the keywords ``fj`` and ``boltzfac``, which
specify the fractional level population for the lower level and the
Boltzmann factor between the two levels of the transition,
respectively.

The are a range of other keywords that affect the behavior of these
computation routines. See :ref:`sssec-full-pwind` for a full
listing.

