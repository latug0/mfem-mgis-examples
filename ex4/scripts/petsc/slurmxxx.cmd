#!/bin/sh
#SBATCH -J mfem
#SBATCH -p skylake
#SBATCH -n FFF
#SBATCH --ntasks-per-node=CCC
#SBATCH --ntasks-per-core=1
#SBATCH -A b171
#SBATCH -t LLL
#SBATCH -o mmm.%j.%a.out
#SBATCH -e mmm.%j.%a.err
#SBATCH --exclusive

# chargement des modules
module purge 
source /home/glatu/mod_mgis.sh

ulimit -a

MPIOPT="-report-bindings  --map-by core -bind-to core"
COMMONOPT="" #" --refine 2"
#time mpirun ${MPIOPT} ptest_full_periodic  --mesh ${MFEMMGIS_DIR}/tests/cube_2mat_per.mesh 2>&1
time mpirun -n ${SLURM_NTASKS}  ${MPIOPT} EEE  ${COMMONOPT} 2>&1
