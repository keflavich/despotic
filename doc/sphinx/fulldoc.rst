.. _sec-fulldoc:

Full Documentation of All DESPOTIC Classes and Functions
========================================================

despotic classes
----------------

.. _sssec-full-cloud:

``cloud``
^^^^^^^^^

.. autoclass:: despotic.cloud
   :members:

``collPartner``
^^^^^^^^^^^^^^^

.. autoclass:: despotic.collPartner
   :members:

``composition``
^^^^^^^^^^^^^^^

.. autoclass:: despotic.composition
   :members:

``despoticError``
^^^^^^^^^^^^^^^^^

.. autoclass:: despotic.despoticError
   :members:

``dustProp``
^^^^^^^^^^^^

.. autoclass:: despotic.dustProp
   :members:

.. _sssec-full-emitter:

``emitter``
^^^^^^^^^^^

.. autoclass:: despotic.emitter
   :members:

.. _sssec-full-emitterData:

``emitterData``
^^^^^^^^^^^^^^^

.. autoclass:: despotic.emitterData
   :members:

``radiation``
^^^^^^^^^^^^^
.. autoclass:: despotic.radiation
   :members:

.. _sssec-full-zonedcloud:

``zonedcloud``
^^^^^^^^^^^^^^
.. autoclass:: despotic.zonedcloud
   :members:


despotic functions
------------------

``fetchLamda``
^^^^^^^^^^^^^^

.. autofunction:: despotic.fetchLamda

.. _sssec-full-lineProfLTE:

``lineProfLTE``
^^^^^^^^^^^^^^^

.. autofunction:: despotic.lineProfLTE
		  
.. _sssec-full-refreshLamda:

``refreshLamda``
^^^^^^^^^^^^^^^^

.. autofunction:: despotic.refreshLamda

		  
despotic.chemistry classes
--------------------------

.. _sssec-full-abundanceDict:

``abundanceDict``
^^^^^^^^^^^^^^^^^

.. autoclass:: despotic.chemistry.abundanceDict
   :members:

``chemNetwork``
^^^^^^^^^^^^^^^

.. autoclass:: despotic.chemistry.chemNetwork
   :members:

``cr_reactions``
^^^^^^^^^^^^^^^^

.. autoclass:: despotic.chemistry.cr_reactions
   :members:

``NL99``
^^^^^^^^

.. autoclass:: despotic.chemistry.NL99
   :members:

``NL99_GC``
^^^^^^^^^^^

.. autoclass:: despotic.chemistry.NL99_GC
   :members:

``photoreactions``
^^^^^^^^^^^^^^^^^^
.. autoclass:: despotic.chemistry.photoreactions

``reaction_matrix``
^^^^^^^^^^^^^^^^^^^

.. autoclass:: despotic.chemistry.reaction_matrix
   :members:


despotic.chemistry functions
----------------------------

``chemEvol``
^^^^^^^^^^^^

.. autofunction:: despotic.chemistry.chemEvol

``setChemEq``
^^^^^^^^^^^^^

.. autofunction:: despotic.chemistry.setChemEq

Shielding functions
^^^^^^^^^^^^^^^^^^^

.. autofunction:: despotic.chemistry.shielding.fShield_H2_DB
.. autofunction:: despotic.chemistry.shielding.fShield_CO_vDB
		  

despotic.winds classes
----------------------

.. _sssec-full-pwind:

``pwind``
^^^^^^^^^
.. autoclass:: despotic.winds.pwind
   :members:
      
   .. automethod:: despotic.winds.pwind.__init__
		      
despotic.winds functions
------------------------

``sxMach``
^^^^^^^^^^

.. autofunction:: despotic.winds.pwind_util.sxMach

``zetaM``
^^^^^^^^^

.. autofunction:: despotic.winds.pwind_util.zetaM

``zetaA``
^^^^^^^^^

.. autofunction:: despotic.winds.pwind_util.zetaA

``pM``
^^^^^^

.. autofunction:: despotic.winds.pwind_util.pM

``pA``
^^^^^^

.. autofunction:: despotic.winds.pwind_util.pA

``tX``
^^^^^^

.. autofunction:: despotic.winds.pwind_util.tX
