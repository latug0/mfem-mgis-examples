#!/bin/bash

set -euo pipefail

OUT="aggregation.csv"

# Header
echo "proc,dof,problem_footprint_GB,solution_footprint_GB,\
mesh_ctor_calls,mesh_ctor_min,mesh_ctor_mean,mesh_ctor_max,mesh_ctor_part,\
mesh_load_calls,mesh_load_min,mesh_load_mean,mesh_load_max,mesh_load_part,\
fed_ctor_calls,fed_ctor_min,fed_ctor_mean,fed_ctor_max,fed_ctor_part,\
thermal_calls,thermal_min,thermal_mean,thermal_max,thermal_part,\
u3si2_swelling_calls,u3si2_swelling_min,u3si2_swelling_mean,u3si2_swelling_max,u3si2_swelling_part,\
mechanics_calls,mechanics_min,mechanics_mean,mechanics_max,mechanics_part,\
total" > "$OUT"


# Extract the 5 numerical values from a timer row:
# number Of Calls | min(s) | mean(s) | max(s) | part(%)
extract_timer() {
    local pattern="$1"
    local file="$2"
    local values

    values=$(
        grep -m1 "$pattern" "$file" \
        | grep -oE '[0-9]+(\.[0-9]+)?' \
        | head -5 \
        | paste -sd, \
        || true
    )

    if [[ -z "$values" ]]; then
        echo "0,0,0,0,0"
    else
        echo "$values"
    fi
}


# Main loop
for f in $(ls Thermo_scaling_*.o | sort -V); do

    echo "Processing $f"

    proc=$(echo "$f" | sed -E 's/^Thermo_scaling_([0-9]+).*/\1/')

    dof=$(
        grep -m1 "Number of finite element unknowns" "$f" \
        | grep -oE '[0-9]+$' \
        || echo 0
    )

    problem_mem=$(
        grep -m1 "After_problem_creation" "$f" \
        | grep -oE '[0-9]+(\.[0-9]+)?' \
        | tail -1 \
        || echo 0
    )

    solution_mem=$(
        grep -m1 "After Solving" "$f" \
        | grep -oE '[0-9]+(\.[0-9]+)?' \
        | tail -1 \
        || echo 0
    )


    # Timers from the timetable

    mesh_ctor=$(extract_timer "Mesh::Constructor" "$f")

    mesh_load=$(extract_timer "Mesh::LoadMesh" "$f")

    fed_ctor=$(extract_timer "FED::Constructor" "$f")


    # Physics solved inside IterativeCouplingScheme::computeNextState
    thermal=$(extract_timer "|    |--> Thermal" "$f")

    u3si2_swelling=$(
        extract_timer "U3SI2_SolidSwelling_Tridimensional" "$f"
    )

    mechanics=$(extract_timer "|    |--> Mechanics" "$f")


    # Total execution time

    total=$(
        grep -m1 "> root" "$f" \
        | grep -oE '[0-9]+(\.[0-9]+)?' \
        | sed -n '3p' \
        || true
    )

    total=${total:-0}

    # Write row

    echo "${proc},${dof},${problem_mem},${solution_mem},\
${mesh_ctor},\
${mesh_load},\
${fed_ctor},\
${thermal},\
${u3si2_swelling},\
${mechanics},\
${total}" >> "$OUT"

done


echo
echo "Aggregation written to $OUT"
