Architecture
============

Data Flow
---------

.. code-block:: text

   Input .dat/.conti file
     → Parser (Flex/Bison: src/Parser/Parser.l + Parser.y)
     → DataSet (Mesh, Fields, Materials, BCs, ICs)
     → Models  (physics kernels, one per material model)
     → Modules (Monolithic or SNIA time-stepping algorithm)
     → Solver  (matrix assembly + sparse linear solve)
     → Outputs (GMSH, VTK, result views)

Source Tree
-----------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Path
     - Role
   * - ``src/Main/Entry.c``
     - Program entry point; session management; MPI init/finalize
   * - ``src/Common/``
     - Memory, buffers, exceptions, message utilities
   * - ``src/DataSet/``
     - Mesh, nodes, elements, fields, materials, BCs, ICs
   * - ``src/Parser/``
     - Flex lexer + Bison grammar for input files
   * - ``src/Models/ModelFiles/``
     - ~50 physics model implementations (M1–M69, BBMGas, Duracem, …)
   * - ``src/Models/ConstitutiveLaws/``
     - Shared constitutive law components (damage, elasticity, plasticity)
   * - ``src/Models/DataBases/``
     - Material property databases (cement chemistry, diffusion coefficients, etc.)
   * - ``src/Models/Methods/``
     - FEM, FVM, and MPM discretization kernels
   * - ``src/Modules/ModuleFiles/``
     - ``Monolithic.c`` and ``SNIA.c`` — time-stepping algorithms
   * - ``src/Solver/``
     - Matrix representation, storage formats, solver dispatch
   * - ``src/Solver/ResolutionMethods/``
     - SuperLU (serial/MT/Dist), PETSc KSP, MA38, Crout
   * - ``src/Tools/``
     - Curve/function evaluation, math utilities, expression parser
   * - ``src/Parallelization/``
     - MPI and OpenMP wrappers
   * - ``src/Libraries/autodiff/``
     - C++17 automatic differentiation library

Object-Oriented-in-C Pattern
-----------------------------

Bil uses an OO-in-C idiom: structs with function pointers, typed with the ``_t``
suffix (e.g. ``Model_t``, ``Solver_t``, ``FEM_t``).

Parallelization
---------------

- **MPI**: domain decomposition via ``src/Parallelization/``. Requires ``ENABLE_MPI=ON`` and PETSc or SuperLU_DIST.
- **OpenMP**: element-level parallelism. Requires ``ENABLE_OPENMP=ON`` and SuperLU_MT if multi-threaded solves are needed.
- The three SuperLU variants (serial, MT, Dist) are mutually exclusive.
