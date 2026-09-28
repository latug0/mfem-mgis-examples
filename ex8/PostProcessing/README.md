# Post-Processing

> This directory is an archive of the scripts written by Julien Rigal during his internship. They are kept for reference and are not necessarily usable as is.

This directory contains all the tools, scripts (Bash, Python), and data used to analyze the computational performance and physical results of the simulations.

## Directory Structure

The directory is divided into two thematic subdirectories:

### HPC_scaling

Dedicated to analyzing code performance and scalability on a supercomputer.

- **Contents:** Scripts for aggregating computation times, calculating speedup, and plotting Strong Scaling and Weak Scaling curves.

### SolverPreconditioner

Dedicated to profiling and optimizing linear algebra operations.

- **Contents:** Log extraction scripts (e.g., `aggregation_SvPc_mech.sh`), result CSV files, and Python scripts for comparing the impact of different solver/preconditioner combinations (HyprePCG, MUMPS, etc.).

## Prerequisites

To run the scripts contained in these directories, the following environment is generally required:

- Python 3
- Libraries: `pandas`, `matplotlib`, `numpy`
- Standard Bash tools (`grep`, `sed`, `awk`)
