.. toctree::
   :numbered:

Potentials and fields
*********************

.. automodule:: scarce.fields

Decription class
----------------

The challenge for the field determination is that numerical differentiation of the potential amplifies numerical instabilities, 
thus the potential has to be smoothed before differentiation. A convinient interface is provided by the field description class.

.. autoclass:: scarce.fields.Description

Methods
-------

.. autofunction:: scarce.fields.calculate_planar_sensor_w_potential

.. autofunction:: scarce.fields.calculate_3D_sensor_w_potential

.. autofunction:: scarce.fields.get_weighting_potential_analytic

.. autofunction:: scarce.fields.get_weighting_field_analytic

.. autofunction:: scarce.fields.get_potential_planar_analytic_1D

.. autofunction:: scarce.fields.get_electric_field_analytic








