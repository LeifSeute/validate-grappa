#!/bin/bash

PEPGEN_ENV=$1
GRAPPA_ENV=$2
FORCEFIELD=$3 # can be amber99sbildn or grappa

echo "PEPGEN_ENV: $PEPGEN_ENV"
echo "GRAPPA_ENV: $GRAPPA_ENV"
echo "FORCEFIELD: $FORCEFIELD"

# throw upon error:
set -e

# check if forcefield is valid
if [ $FORCEFIELD != "amber99sbildn" ] && [ $FORCEFIELD != "grappa" ]; then
    echo "forcefield must be amber99sbildn or grappa but is $FORCEFIELD"
    exit 1
fi

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

FASTA='pep.fasta'

# name of the fasta file:
DIR=equilibration

source ~/.bashrc

# set up pepgen
conda activate $PEPGEN_ENV
echo $CONDA_PREFIX

echo "fasta filename:"
echo $FASTA
seq=$( tail -n 1 $FASTA )
echo $seq

# run pepgen in uncapped mode:
pepgen $seq $DIR -t

conda activate $GRAPPA_ENV

# create gmx files
cp -r runfiles_initial/* $DIR/
pushd $DIR

# (the 6 1 flags are to select the amber99sbildn forcefield and water model)
printf "6\n1\n "|gmx pdb2gmx -f pep.pdb -o pep.gro -p pep.top -ignh

if [ $FORCEFIELD == "grappa" ]; then
    # parametrize with grappa
    mv pep.top pep_amber99sbildn.top
    grappa_gmx -f pep_amber99sbildn.top -o pep.top -t grappa-1.1.0 -c classical

fi

# copy the topology before solvation and ion addition:
cp pep.top pep_unsolvated.top

# run mds
bash quick_gmx.sh

echo "Done with equilibration! Now pick a structure that fits in the box described in the folding paper SI (different for every protein)."
popd