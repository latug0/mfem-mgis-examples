# Thermomechanical Simulation with MFEM-MGIS

## Project Overview

This project implements a nonlinear thermomechanical coupling solver based on the **MFEM-MGIS** library. It is designed to simulate the complex behavior of materials (such as **U3Si2** fuel and **ALFENI** cladding) subjected to intense thermal and mechanical loads.

The code manages a strong coupling (via an `IterativeCouplingScheme`) between two physics:
1. **Heat Transfer:** Resolution of the heat equation, taking into account a time-varying power source term.
2. **Mechanics:** Resolution of the mechanical equilibrium integrating complex behaviors generated via **MFront** (thermal expansion, solid swelling, Norton PRQ viscoplasticity, and plasticity with hardening).

The code is optimized for parallel computing (MPI) and allows the use of advanced iterative solvers (Hypre family) or direct solvers depending on the physics involved.

---

## Prerequisites and Installation

To compile and run this project, you need the **MFEM-MGIS** environment, which relies on **MFEM** (finite elements) and **TFEL/MGIS** (integration of MFront material behaviors).

### Main Dependencies
* **MPI** (OpenMPI, MPICH, etc.) for parallel computing.
* **TFEL/MFront** (compiled with generic interfaces).
* **MFEM** (compiled with MPI, Hypre, and ideally Metis support).
* **MGIS** (MFront Generic Interface Support).

### Installation Guide
To install the library and its dependencies, please refer to the official documentation:
**[MFEM-MGIS Installation Guide](https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html)**

Once the environment is set up, you can compile this thermomechanical project by linking your `CMakeLists.txt` to the `mfem-mgis` installation.

---

## Usage and Command-Line Arguments

The main program accepts several command-line arguments to configure the mesh, MFront behaviors, and linear solvers. 

Here is a summary table of the available options:

| Short Option | Long Option | Description |
| :--- | :--- | :--- |
| `-m` | `--mesh` | Path to the mesh file (e.g., `.msh`, `.vtk`). |
| `-lU` | `--libraryU3SI2` | Path to the compiled MFront library (`.so`) for the U3Si2 material. |
| `-lA` | `--libraryALFENI` | Path to the compiled MFront library (`.so`) for the ALFENI material. |
| `-svTh` | `--solverTh` | Name of the linear solver to use for the heat transfer problem (e.g., `HypreGMRES`). |
| `-pcTh` | `--preconditionnerTh` | Preconditioner associated with the thermal solver (e.g., `HypreBoomerAMG`). |
| `-svMc` | `--solverMc` | Name of the linear solver to use for the mechanics problem (e.g., `HyprePCG`). |
| `-pcMc` | `--preconditionnerMc` | Preconditioner associated with the mechanics solver. |
| `-o` | `--order` | Finite element order (polynomial degree, default is usually 1). |
| `-r` | `--refinement` | Uniform refinement level of the mesh (default: `0`). |
| `-p` | `--post-processing` | Enables (`1`) or disables (`0`) the export of results for ParaView. |
| `-v` | `--verbosity-level` | Verbosity level of the console logs (`0` = minimal, higher levels = increased details). |

### Parallel Execution Example

```bash
mpirun -np 4 ./Thermomechanical \
  -m ../assemblage_hexa.msh \
  -lU src/libU3SI2-generic.so \
  -lA src/libALFENI-generic.so \
  -svTh HypreGMRES -pcTh HypreBoomerAMG \
  -svMc HyprePCG -pcMc HypreBoomerAMG \
  -o 1 -r 0 -p 1 -v 1