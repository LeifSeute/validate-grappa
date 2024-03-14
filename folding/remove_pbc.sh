#!/bin/bash

# usage ./remove_pbc.sh <dir>
# <dir>: directory containing the md.tpr and md.trr files
# This script removes periodic boundary conditions (PBC).
# For a given directory, it will look for the files md.tpr and md.trr, and if they exist, it will remove PBC from the trajectory file and save it as md_noPBC.trr.

set -e

# Navigate to the parent directory containing all MD simulation directories
DIR=$1
echo removing PBC from $DIR
OUTPUT_GROUP=${2:-"System"} # set to Protein to only write the protein positions in the output file
pushd "$DIR"

# module load GCC/9.3.0
# module load CUDA/11.5.0

# module use /hits/fast/mbm/broszms/software/modules/
# module load gromacs-2023-cuda-11.5

# Loop through each child directory in the folder

# Check if the necessary files exist before proceeding
if [[ -f "md.tpr" && -f "md.trr" ]]; then
    # Automatically select the Protein group (usually group number 1) for centering and for output by echoing their numbers into gmx trjconv

    # if output exists, delete it:
    if [[ -f "md_centered.trr" ]]; then
        rm md_centered.trr
    fi
    
    echo -e "Protein\n$OUTPUT_GROUP" | gmx trjconv -s md.tpr -f md.trr -o md_centered.trr -pbc mol -center -ur compact

else
    echo "ERROR DURING REMOVE_PBC: Required files md.tpr or md.trr not found in $Dir"
fi

# Return to the original directory
popd
