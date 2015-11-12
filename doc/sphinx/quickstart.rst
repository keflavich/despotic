.. highlight:: rest

Quickstart
==========

Introduction
------------

DESPOTIC is a tool to Derive the Energetics and SPectra of
Optically Thick Interstellar Clouds. It can perform a
variety of calculations regarding the chemical and thermal state of
interstellar clouds, and predict their observable line
emission. DESPOTIC treats clouds in a simple one-zone model, and is
intended to allow rapid, interactive exploration of parameter space.

In this Quickstart, we will walk through a basic interactive python
session using DESPOTIC. This will work equally well from an ipython
shell or in an ipython notebook. 

The ``cloud`` Class
-------------------

The basic object in DESPOTIC, which provides an interface to most of
its functionality, is the class ``cloud``. This class stores the basic
properties of an interstellar cloud, and provides methods to perform
calculations on those properties. The first step in most DESPOTIC
sessions is to import this class:: 

   from despotic import cloud

The next step is generally to input data describing a cloud upon which
computations are to be performed. The input data describe the cloud's
physical properties (density, temperature, etc.), the bulk composition
of the cloud, what emitting species it contains, and the radiation
field around it. While it is possible to enter the data manually, it
is usually easier to read the data from a file, using the
:ref:`sec-cloudfiles` format. For this Quickstart, we'll use one of
the configuration files that ship with DESPOTIC and that are
included in the ``cloudfiles`` subdirectory of the DESPOTIC
distribution. To create a cloud whose properties are as given in
a particular cloud file, we simply invoke the constructor with
the ``fileName`` optional argument set equal to a string containing
the name of the file to be read::

  gmc = cloud(fileName="cloudfiles/MilkyWayGMC.desp", verbose=True)

Note the optional argument ``verbose``, which we have set to
``True``. Most DESPOTIC methods accept the ``verbose`` argument, which
causes them to produce printed output containing a variety of
information. By default DESPOTIC operations are silent. 

Computing Temperatures
----------------------

At this point most of the calculations one could want to do on a cloud
are provided as methods of the ``cloud`` class. One of the most basic is
to set the cloud to its equilibrium dust and gas temperatures. This is
accomplished via the ``setTempEq`` method::

  gmc.setTempEq(verbose=True)

With ``verbose`` set to ``True``, this command will produce variety of
output as it iterates to calculate the equilibrium gas and dust
temperatures, before finally printing ``True``. This illustrates
another feature of DESPOTIC commands: those that iterate return a
value of ``True`` if they converge, and ``False`` if they do not.

To see the gas and dust temperatures to which the cloud has been set,
we can simply print them::

  print gmc.Tg
  print gmc.Td

This shows that DESPOTIC has calculated an equilibrium gas temperature
of 10.2 K, and an equilibrium dust temperature of 14.4 K.

Line Emission
-------------

Next we might wish to compute the CO line emission emerging from the
cloud. We do this with the ``cloud`` method ``lineLum``::

  lines = gmc.lineLum("co")

The argument ``co`` specifies that we are interested in the emission
from the CO molecule. This method returns a ``list`` of ``dict``, each
of which gives information about one of the CO lines. The ``dict``
contains a variety of fields, but one of them is the
velocity-integrated brightness temperature of the line. Again, we can
just print the values we want. The first element in the list is the
:math:`J = 1 \rightarrow 0` line, and the velocity-integrated
brightness temperature is listed as ``intTB`` in the ``dict``. Thus to
get the velocity-integrated brightness temperature of the first line,
we just do::

  print lines[0][’intTB’]

This shows that the velocity-integrated brightness temperature of the
CO :math:`J = 1 \rightarrow 0` line is 79 K km/s.

Heating and Cooling Rates
-------------------------

Finally, we might wish to know the heating and cooling rates produced
by various processes, which lets us determined what sets the thermal
balance in the cloud. This may be computed using the method ``dEdt``,
as follows::

  rates = gmc.dEdt()

This method returns a ``dict`` that contains all the heating and
cooling terms for gas and dust. For example, we can print the rates of
cosmic ray heating and CO line cooling via::

  print rates["GammaCR"]
  print rates["LambdaLine"]["co"]
