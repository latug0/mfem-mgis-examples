# Robin Thermal Convection Test

**Author:** Pierre Dubois

## Overview

This test case validates the implementation of thermal convection through a **Robin boundary condition** in MFEM-MGIS.

The problem consists of a steady-state heat conduction test in a three-dimensional bar with constant thermal conductivity. A prescribed temperature is imposed on one end of the domain, while a Robin (convective) boundary condition is applied on the opposite end:

\[
\mathbf{j} \cdot \mathbf{n} = h \left(T - T_\infty\right),
\]

where \(h\) is the heat transfer coefficient and \(T_\infty\) is the ambient temperature.

Because the analytical solution is known, this example provides a straightforward verification of the implementation by comparing the numerical and exact solutions. It also validates the residual and tangent contributions associated with the Robin boundary condition.