Installation
============

Prerequisites
-------------

- CMake >= 3.10
- C/C++ compiler (GCC, Clang, or MSVC)
- Optional: MPI, PETSc, SuperLU, FLTK, OpenMP

Build
-----

.. code-block:: bash

   git clone https://github.com/Universite-Gustave-Eiffel/bil
   cd bil
   mkdir build && cd build
   cmake ..
   make

The binary is placed in ``../bin/``, the shared library in ``../lib/``.

Build Options
^^^^^^^^^^^^^

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE=Release ..   # optimised build
   cmake -DCMAKE_BUILD_TYPE=Debug ..     # debug build

Feature flags (edit ``OPTIONS`` at the project root):

.. list-table::
   :widths: 40 15
   :header-rows: 1

   * - Option
     - Default
   * - ``ENABLE_MPI``
     - ON
   * - ``ENABLE_AUTODIFF``
     - ON
   * - ``ENABLE_FLTK``
     - ON
   * - ``ENABLE_PETSC``
     - ON
   * - ``ENABLE_SUPERLU``
     - ON
   * - ``ENABLE_OPENMP``
     - OFF
