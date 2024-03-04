#! bin/bash

name='protein_g'
BOX_SIZE=5.5

set -e

BASE_ID=35879

JOB_NAMES=("1a${name}")

dirname="../results/${name}_amber99sbildn"

cd "$dirname"

# HERE, WE NEED TO CHANGE THE ADD IONS BECAUSE WE WANT TO ADD BOTH NA AND CL. check this:

if grep -q -- "-pname NA -nname CL -neutral -conc 0.1" submit_continue.sh; then
    echo "Script is fine. Will submit jobs."
else
    echo "You add ions with a given concentration. Please replace the genion cmd with"
    echo 'echo "SOL" | gmx genion -s pep_genion.tpr -p pep.top -o pep_ion.gro -pname NA -nname CL -neutral -conc 0.1'
    exit 1
fi


for i in "${!JOB_NAMES[@]}"; do
    jobname="${JOB_NAMES[$i]}"
    # add 1 to BASE_ID
    BASE_ID=$((BASE_ID + 1))
    id="$BASE_ID"
    
    echo "Submitting job $i with ID $id and jobname $jobname..."
    sbatch -J "$jobname" "submit_continue.sh" "$id" "$BOX_SIZE"
done