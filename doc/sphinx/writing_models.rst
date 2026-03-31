Writing a New Model
===================

Step 1 — Create the file
-------------------------

Create ``src/Models/ModelFiles/M<N>.c`` where ``<N>`` is the next available
model index.

Step 2 — Required boilerplate
------------------------------

.. code-block:: c

   #include <stdio.h>
   #include <stdlib.h>
   #include <string.h>
   #include <math.h>
   #include "CommonModel.h"

   #define TITLE    "My Model"
   #define AUTHORS  "Your Name"

   #include "PredefinedModelMethods.h"

   /* Number of equations, explicit terms, implicit terms */
   #define NEQ  (1)
   #define NVE  (1)
   #define NVI  (3)

Step 3 — Implement required functions
--------------------------------------

``pm(const char *s)``
   Maps a parameter name (from the input file) to an integer index.

``SetModelProp(Model_t *model)``
   Registers the model title, number of equations, etc.

``ComputeInitialState(Element_t *el, double t)``
   Initialises nodal values at ``t=0``.

``ComputeImplicitTerms(Element_t *el, double t, double dt)``
   Updates implicit (storage) terms.

``ComputeMatrix(Element_t *el, double t, double dt, double *k)``
   Fills the element stiffness/tangent matrix.

``ComputeResidu(Element_t *el, double t, double dt, double *r)``
   Fills the element residual vector.

Step 4 — Register the model
-----------------------------

Add the model index to the registry in ``src/Models/``.

Step 5 — Validate
------------------

Create an input file in ``examples/`` and a reference solution in ``base/``,
then run:

.. code-block:: bash

   bil myfile
