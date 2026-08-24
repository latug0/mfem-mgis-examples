#!/bin/bash

REF=2

PROC_LIST=(
32
64
128
256
512
1024
2048
4096
8192
16384
)

EXEC="../../build/Thermomechanical"

echo "Strong Scaling starting (Refinement: $REF)"

for NPROC in "${PROC_LIST[@]}"; do

    JOB_FILE="ThMc_scaling_${NPROC}.sh"

    echo "-> Generate and submit for $NPROC ranks..."

    cat <<EOF > $JOB_FILE
#!/bin/bash

#MSUB -n $NPROC
#MSUB -T 180
#MSUB -q milan
#MSUB -m work,scratch

# module purge
# module load gnu/12.3.0 mpi cmake/3.29.6

echo "Run starts"
echo "Refinement : $REF"

ccc_mprun $EXEC --mesh "../../partitionning/CompleteMesh${NPROC}/output-mesh${NPROC}." --refinement $REF

echo "Run ends : \$(date)"
EOF

    ccc_msub $JOB_FILE

    sleep 1

done

echo "All tasks are in the queue."
