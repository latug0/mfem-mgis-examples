#!/bin/bash

echo "Generating jobs for the study of mechanics solver/preconditioner..."

BASE_DIR="/ccc/scratch/cont002/den/rigaljul/Projet/ThermoMécanique/SwellingAdded"
EXEC="$BASE_DIR/build/ThermoMecaFancy"
LIB_U="$BASE_DIR/build/src/libU3SI2-generic.so"
LIB_A="$BASE_DIR/build/src/libALFENI-generic.so"

LOG_DIR="logs_SvPc_mech"
JOB_DIR="jobs_scripts"
mkdir -p "$LOG_DIR"
mkdir -p "$JOB_DIR"

MESH="/ccc/scratch/cont002/den/rigaljul/Partitioning/CompleteMesh128/output-mesh128."
SOLVERS=("HyprePCG"
	 "HypreGMRES"
	 "HypreFGMRES"
	 "MINRESSolver"
	 "CGSolver"
	 "MUMPSSolver"
	 "BiCGSTABSolver")

PRECONDS=("HypreBoomerAMG"
	  "HypreDiagScale"
	  "HypreILU"
	  "HypreParaSails")

COUNT=0

for SV_MC in "${SOLVERS[@]}"; do

    if [[ "$SV_MC" == *"MUMPS"* ]]; then
       PRECONDS_MC=("NONE")
    else
       PRECONDS_MC=("${PRECONDS[@]}")
    fi

    for PC_MC in "${PRECONDS_MC[@]}"; do
        ((COUNT++))

        JOB_FILE="$JOB_DIR/job_sweep_${COUNT}.sh"
        LOG_FILE="$LOG_DIR/run_${MESH_NAME}_Th_HypreGMRES_HypreBoomerAMG_Mc_${SV_MC}_${PC_MC}.log"
        echo "-> Submission [$COUNT] : $MESH_NAME | Th: $HypreGMRES+$HypreBoomerAMG | Mc: $SV_MC+$PC_MC"

        cat <<EOF > $JOB_FILE
#!/bin/bash
#MSUB -n 128
#MSUB -T 3600
#MSUB -q milan
#MSUB -m work,scratch
#MSUB -J Swp_${COUNT}
#MSUB -o ${LOG_FILE}.out
#MSUB -e ${LOG_FILE}.err

echo "Exécution de la configuration $COUNT"

ccc_mprun "$EXEC" -m "$MESH" \\
    --libraryU3SI2 "$LIB_U" --libraryALFENI "$LIB_A" \\
    -r 0 \\
    -svTh "HypreGMRES" -pcTh "HypreBoomerAMG" \\
    -svMc "$SV_MC" -pcMc "$PC_MC" > "$LOG_FILE" 2>&1

if [ \${PIPESTATUS[0]} -ne 0 ]; then
    echo "$MESH_NAME,$HypreGMRES,$HypreBoomerAMG,$SV_MC,$PC_MC,FAILED,FAILED,FAILED" > "$CSV_PART"
fi
EOF

        ccc_msub $JOB_FILE
        sleep 0.5
    done
done

echo -e "\n===================================================================="
echo "All submissions ($COUNT jobs) are in the queue !"
echo "===================================================================="
echo "Once all jobs are finished, type this command to aggregate the results :"
echo ""
echo "bash aggregation_SvPc_mech.sh"
echo "===================================================================="
