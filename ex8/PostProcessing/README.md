# Post-Processing

This directory contains all the tools, scripts (Bash, Python), and data used to analyze the computational performance and physical results of the simulations.

## Directory Structure

The directory is divided into three thematic subdirectories:

### HPC_scaling

Dedicated to analyzing code performance and scalability on a supercomputer.

- **Contents:** Scripts for aggregating computation times, calculating speedup, and plotting Strong Scaling and Weak Scaling curves.

### SolverPreconditioner

Dedicated to profiling and optimizing linear algebra operations.

- **Contents:** Log extraction scripts (e.g., `aggregation_SvPc_mech.sh`), result CSV files, and Python scripts for comparing the impact of different solver/preconditioner combinations (HyprePCG, MUMPS, etc.).

### Physics_postprocessing

Dedicated to the physical and thermomechanical analysis of the simulation results.

- **Contents:** Processing of physical fields (temperature, strain), extraction of MFront state variables (U3Si2 swelling, ALFENI plasticity), and generation of time-evolution plots.

## Prerequisites

To run the scripts contained in these directories, the following environment is generally required:

- Python 3
- Libraries: `pandas`, `matplotlib`, `numpy`
- Standard Bash tools (`grep`, `sed`, `awk`)
