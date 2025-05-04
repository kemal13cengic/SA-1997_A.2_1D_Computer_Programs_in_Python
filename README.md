# 1D Computational Fluid Dynamics Programs - Fortran to Python Translation

This project is a translation of 1D computer programs from **Fortran** (as found in the SA-1997 Lecture Notes for the course *Introduction to Computational Fluid Dynamics* by I. Demirdžić and S. Muzaferija) to **Python**.

The programs in the appendix are designed to solve 1D problems of:

- Linear and non-linear steady heat conduction (**Program A.2.1 Steady heat conduction**)
- Steady heat conduction-convection (**Program A.2.2 Steady heat conduction-convection**)
- Transient heat conduction (**Program A.2.3 Transient heat conduction**)

These solvers utilize **linear equation solvers** based on either **TDMA** (Tri-Diagonal Matrix Algorithm) or **GSM** (Gauss-Seidel Method).

There are also **arbitrary precision versions** of all three programs, utilizing the `mpmath` library, available in a separate directory.
