#!/bin/bash

# create a run folder, then create an unfolded peptide, run a short equilibration md, pick the state with smallest maximum distance between any two c alpha atoms. creates a figure of the picked state in a box, where one can confirm that the peptide is in the box and does not hit the box due to periodic boundary conditions.

FASTA=$1 # path to fasta file
PEPGEN_ENV=${2:-"pepgen"}
GRAPPA_ENV=${3:-"grappa"}
FORCEFIELD=${4:-"amber99sbildn"} # can be amber99sbildn or grappa

# throw upon error:
set -e

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

DIR=$(basename $FASTA .fasta)

# concat the ff_name to the dir name:
DIR=$DIR"_"$FORCEFIELD

mkdir $DIR
mkdir $DIR/folding
mkdir -p $DIR/logfiles
mkdir -p $DIR/mds


cp -r setup_template/* $DIR/
pushd $DIR

bash create_peptide.sh $PEPGEN_ENV $GRAPPA_ENV $FORCEFIELD

# pick a collapsed structure
popd

# remove periodic boundary conditions
bash remove_pbc.sh $DIR/equilibration

python pick_collapsed.py -gro $DIR/equilibration/md.gro -trr $DIR/equilibration/md.trr -o $DIR/folding/pep.gro -fig $DIR

cp $DIR/equilibration/pep.top $DIR/folding/pep.top

echo finished equilibration. take a look at the figures $DIR/collapsed_structure.png

echo if the peptide is in the box, without periodic boundary condition effects that make it look like it is outside the box, then you can proceed with the folding simulations by running sbatch -J jobname submit_continue.sh taskid box_side