#!/bin/bash

# module load GCC/9.3.0
# module load CUDA/11.5.0

# module use /hits/fast/mbm/broszms/software/modules/
# module load gromacs-2023-cuda-11.5

# Navigate to the parent directory containing all MD simulation directories
PARENT_DIR=$1
pushd "$PARENT_DIR"

# Loop through each child directory in ../mds
for dir in */; do
    echo "Processing directory: $dir"
    pushd "$dir"
    
    # Check if the necessary files exist before proceeding
    if [[ -f "md.tpr" && -f "md.trr" ]]; then
        # Automatically select the Protein group (usually group number 1) for centering and for output by echoing their numbers into gmx trjconv
        echo -e "1\n1" | gmx trjconv -s md.tpr -f md.trr -o md_noPBC.trr -pbc mol -center -ur compact
    else
        echo "Required files md.tpr or md.trr not found in $dir"
    fi
    
    # Return to the parent directory
    popd
done

# Return to the original directory
popd
