#! bin/bash

name='2a3d'
BOX_SIZE=6.4

set -e

BASE_ID=84509

JOB_NAMES=("1${name}" "2${name}" "3${name}")

dirname="../results/${name}_grappa"

cd "$dirname"

for i in "${!JOB_NAMES[@]}"; do
    jobname="${JOB_NAMES[$i]}"
    # add 1 to BASE_ID
    BASE_ID=$((BASE_ID + 1))
    id="$BASE_ID"
    
    echo "Submitting job $i with ID $id and jobname $jobname..."
    sbatch -J "$jobname" "submit_continue.sh" "$id" "$BOX_SIZE"
done