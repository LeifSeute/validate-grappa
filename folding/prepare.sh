#!/bin/bash

# create a run folder, then create an unfolded peptide, run a short equilibration md, pick the state with smallest maximum distance between any two c alpha atoms. creates a figure of the picked state in a box, where one can confirm that the peptide is in the box and does not hit the box due to periodic boundary conditions.

FASTA=$1 # path to fasta file
PEPGEN_ENV=$2
GRAPPA_ENV=$3
FORCEFIELD=$4 # can be amber99sbildn or grappa

# throw upon error:
set -e


module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

DIR=$(basename $FASTA .fasta)

# concat the ff_name to the dir name:
DIR=$DIR"_"$FORCEFIELD

mkdir -p results
mkdir results/$DIR
mkdir results/$DIR/folding
mkdir -p results/$DIR/logfiles
mkdir -p results/$DIR/mds
mkdir -p results/$DIR/figures

cp $FASTA results/$DIR/pep.fasta

cp -r setup_template/* results/$DIR/
pushd results/$DIR

bash create_peptide.sh $PEPGEN_ENV $GRAPPA_ENV $FORCEFIELD

# pick a collapsed structure
popd

# remove periodic boundary conditions
bash remove_pbc.sh results/$DIR/equilibration

source ~/.bashrc
conda activate $GRAPPA_ENV
python pick_collapsed.py -gro results/$DIR/equilibration/md.gro -trr results/$DIR/equilibration/md_centered.trr -o results/$DIR/folding/pep.gro -fig results/$DIR/figures

cp results/$DIR/equilibration/pep_unsolvated.top results/$DIR/folding/pep.top

echo finished equilibration. take a look at the figures results/$DIR/figures/collapsed_structure.png

echo if the peptide is in the box, without periodic boundary condition effects that make it look like it is outside the box, then you can proceed with the folding simulations by running sbatch -J jobname submit_continue.sh taskid box_side in results/$DIR