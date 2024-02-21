#!/bin/bash

# This script removes periodic boundary conditions (PBC).
# For a given directory, it will look for the files md.tpr and md.trr, and if they exist, it will remove PBC from the trajectory file and save it as md_noPBC.trr.

set -e

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

# Navigate to the parent directory containing all MD simulation directories
DIR=$1
pushd "$DIR"

# Loop through each child directory in the folder

# Check if the necessary files exist before proceeding
if [[ -f "md.tpr" && -f "md.trr" ]]; then
    # Automatically select the Protein group (usually group number 1) for centering and for output by echoing their numbers into gmx trjconv
    echo -e "Protein\nSystem" | gmx trjconv -s md.tpr -f md.trr -o md_centered.trr -pbc mol -center -ur compact

else
    echo "Required files md.tpr or md.trr not found in $Dir"
fi

# Return to the original directory
popd
