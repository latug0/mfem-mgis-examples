#!/bin/bash
#MSUB -r mfem-mgis-ex6-2 # Request name
#MSUB -n 2
#MSUB -c 1 # 1 socket is reserved for 1 MPI process.
#MSUB -T 8000 # Elapsed time limit in seconds
#MSUB -o ex6_%I.o # Standard output. %I is the job id
#MSUB -e ex6_%I.e # Error output. %I is the job id
#MSUB -q milan # Choosing partition of GPU nodes
set -x
#cd ${BRIDGE_MSUB_PWD}


NAME=cas-cible-1
EXE=perf
MESH=many_spheres.msh
WORKDIR=$CCCSCRATCHDIR/$NAME

cp $EXE $WORKDIR/
cp $MESH $WORKDIR/

cd $WORKDIR

ccc_mprun -n 2 ./$EXE --mesh $MESH
