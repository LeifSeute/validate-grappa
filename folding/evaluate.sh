#!bin/bash

set -e

NAME=${1:-"chignolin"}

REF_RMSD=0
if [ "$NAME" == "chignolin" ]; then
    REF_RMSD=1
elif [ "$NAME" == "chignolin_old" ]; then
    REF_RMSD=1
elif [ "$NAME" == "trp_cage" ]; then
    REF_RMSD=1.4
elif [ "$NAME" == "trp_cage_from_amber" ]; then
    REF_RMSD=1.4
elif [ "$NAME" == "protein_g" ]; then
    REF_RMSD=1.2
elif [ "$NAME" == "w_domain" ]; then
    REF_RMSD=1.2
else
    echo "Unknown protein name: $NAME"
    exit 1
fi

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

# BASE_DIR="/hits/fast/mbm/seutelf/md_simulation/validate_grappa/folding"
BASE_DIR="."

MD_DIR=$BASE_DIR"/results/${NAME}_grappa/mds"

# iterate over subdirs to process the trajectories:
for subdir in $(ls $MD_DIR); do
    echo "Processing $subdir"
    # bash remove_pbc.sh $MD_DIR/$subdir Protein
done

mkdir -p $BASE_DIR/analysis

python analysis.py --ref $BASE_DIR/references/${NAME}.pdb --trjdir $MD_DIR --ref_rmsd $REF_RMSD --figpath $BASE_DIR/analysis/${NAME}_grappa_ref_rmsd.png


# if Base dir / results/name_amber/mds exists and contains at least one subfolder, process it as well:
MD_DIR=$BASE_DIR"/results/${NAME}_amber99sbildn/mds"
if [ -d $MD_DIR ]; then
    counter=0
    for subdir in $(ls $MD_DIR); do
        echo "Processing $subdir"
        # bash remove_pbc.sh $MD_DIR/$subdir Protein
        counter=$((counter+1))
    done

    if [ $counter -eq 0 ]; then
        echo "No subfolders found in $MD_DIR"
    else
        python analysis.py --ref $BASE_DIR/references/${NAME}.pdb --trjdir $MD_DIR --ref_rmsd $REF_RMSD --figpath $BASE_DIR/analysis/${NAME}_amber99sbildn_ref_rmsd.png --amber
    fi
fi