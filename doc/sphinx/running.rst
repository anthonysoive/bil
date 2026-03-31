Running Bil
===========

Basic usage
-----------

.. code-block:: bash

   bil [options] myfile
   bil               # prints available options

Input files
-----------

Input files use the ``.dat`` or ``.conti`` extension. They specify the mesh,
material properties, boundary conditions, initial conditions, and output
requests.

Example:

.. code-block:: bash

   bil examples/M1/richards.dat

Validation
----------

Reference solutions are stored in ``base/`` organised by model. To validate:

.. code-block:: bash

   cd examples/M1
   bil richards.dat
   diff output.txt ../../base/M1/M1.t1

Memory checks
^^^^^^^^^^^^^

.. code-block:: bash

   make memcheck    # runs Valgrind
