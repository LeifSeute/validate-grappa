#! bin/bash

name='chignolin'
BOX_SIZE=4

set -e

BASE_ID=34234

JOB_NAMES=("1${name}" "2${name}") # "3${name}")

dirname="../results/${name}_amber99sbildn"

cd "$dirname"

for i in "${!JOB_NAMES[@]}"; do
    jobname="${JOB_NAMES[$i]}"
    # add 1 to BASE_ID
    # BASE_ID=$((BASE_ID + 1))
    id="$BASE_ID"
    
    echo "Submitting job $i with ID $id and jobname $jobname..."
    sbatch -J "$jobname" "submit_continue.sh" "$id" "$BOX_SIZE"
    # hard code the second id:
    BASE_ID=34893
done