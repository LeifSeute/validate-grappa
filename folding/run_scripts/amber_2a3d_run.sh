#! bin/bash

name='2a3d'
BOX_SIZE=6.4

set -e

BASE_ID=19237

JOB_NAMES=("1${name}")

dirname="../results/${name}_amber99sbildn"

cd "$dirname"

for i in "${!JOB_NAMES[@]}"; do
    jobname="${JOB_NAMES[$i]}"
    # add 1 to BASE_ID
    BASE_ID=$((BASE_ID + 1))
    id="$BASE_ID"
    
    echo "Submitting job $i with ID $id and jobname $jobname..."
    sbatch -J "$jobname" "submit_continue.sh" "$id" "$BOX_SIZE"
done