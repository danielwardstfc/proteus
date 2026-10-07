# Liquid Metal MHD Solver with MFEM backend

This example demonstrates incompressible liquid-metal magnetohydrodynamics (LMMHD) using the MFEM finite-element backend in MOOSE.

The solver couples the fluid velocity and pressure with the electric potential and current density in an imposed magnetic field. The formulation and implementation are based on the finite-element approach described in [Li et al. (2019)](https://doi.org/10.1016/j.cma.2019.112552).

For further information on MFEM in MOOSE, see the [MFEM-MOOSE documentation](https://mooseframework.inl.gov/syntax/MFEM/) and [MFEM](https://mfem.org/).
