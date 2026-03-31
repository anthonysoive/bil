Physics Models
==============

Models live in ``src/Models/ModelFiles/``. Each model is identified by a number
(``M1``, ``M2``, …) or a descriptive name (``BBMGas``, ``Duracem``, …).

Model File Convention
---------------------

Every model file follows this fixed pattern:

.. code-block:: c

   #include "CommonModel.h"
   #define MODELINDEX  <N>
   #define TITLE       "<description>"
   #define AUTHORS     "<name>"
   #include "OldMethods.h"   // or "Methods.h" for newer models

   // pm()   — maps parameter name string → index
   // dm<N>() — reads material parameters from input file
   // qm<N>() — computes fluxes / internal variables
   // sm<N>() — computes source/storage terms
   // k<N>()  — assembles element stiffness matrix
   // c<N>()  — assembles element coupling/capacity matrix

The ``MODELINDEX`` must match the file name number and the registration in the
model registry.

Available Models
----------------

.. list-table::
   :widths: 15 55 30
   :header-rows: 1

   * - Model
     - Description
     - Application domain
   * - M1
     - Richards Equation (1D)
     - Unsaturated flow
   * - M2
     - Richards Equation (generalized)
     - Unsaturated flow
   * - M5
     - Diffusion–advection
     - Transport
   * - M7
     - Poromechanics
     - Coupled HM
   * - M10
     - Cement hydration
     - Cement chemistry
   * - M15
     - Chloride diffusion in concrete
     - Durability
   * - M54
     - Alkali–silica reaction
     - Concrete degradation
   * - BBMGas
     - Barcelona Basic Model with gas
     - Unsaturated soil mechanics
   * - Duracem
     - Durability of cementitious materials
     - Cement chemistry
