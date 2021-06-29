#!/bin/sh
#SBATCH -J mfem
#SBATCH -p skylake
#SBATCH -n 16
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH -A b171
#SBATCH -t 0:49:00
#SBATCH -o hypre.%j.%a.out
#SBATCH -e hypre.%j.%a.err
#SBATCH --exclusive

# chargement des modules
module purge 
source /home/glatu/mod_mgis.sh

ulimit -a

MPIOPT="-report-bindings  --map-by core -bind-to core"
COMMONOPT="" #" --refine 2"
#time mpirun ${MPIOPT} ptest_full_periodic  --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh 2>&1
time mpirun -n 16 ${MPIOPT} ./Ssna303_3d_hypre  ${COMMONOPT} 2>&1
