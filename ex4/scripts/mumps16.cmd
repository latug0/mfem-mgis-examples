#!/bin/sh
#SBATCH -J mfem
#SBATCH -p skylake
#SBATCH -n 16
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-core=1
#SBATCH -A b171
#SBATCH -t 0:29:00
#SBATCH -o ssna303.%j.%a.out
#SBATCH -e ssna303.%j.%a.err
#SBATCH --exclusive

# chargement des modules
module purge 
source /home/glatu/mod_mgis.sh

ulimit -a

MPIOPT="-report-bindings  --map-by core -bind-to core"
COMMONOPT="" #" --refine 2"
#time mpirun ${MPIOPT} ptest_full_periodic  --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh 2>&1
time mpirun -n 16 ${MPIOPT} ./Ssna303_3d -no-up  ${COMMONOPT} 2>&1
