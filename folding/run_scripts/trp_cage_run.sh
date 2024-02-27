#! bin/bash

name='trp_cage'
BOX_SIZE=3.7

set -e

BASE_ID=98547

JOB_NAMES=("1${name}" "2${name}" "3${name}")

dirname="../results/${name}_grappa"

cd "$dirname"


# HERE, WE NEED TO CHANGE THE ADD IONS BECAUSE WE WANT TO ADD BOTH NA AND CL. check this:

if grep -q -- "-pname NA -nname CL -neutral -conc 0.065" submit_continue.sh; then
    echo "Script is fine. Will submit jobs."
else
    echo "You need to add ions with a given concentration. Please replace the genion cmd with"
    echo 'echo "SOL" | gmx genion -s pep_genion.tpr -p pep.top -o pep_ion.gro -pname NA -nname CL -neutral -conc 0.065'
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