#!/bin/bash
PDBNAME=$1 #4-letter-code
GRAPPA_ENV=$2
FORCEFIELD=$3 # can be amber99sbildn or grappa

echo "PDBNAME: $PDBNAME"
echo "GRAPPA_ENV: $GRAPPA_ENV"
echo "FORCEFIELD: $FORCEFIELD"

# throw upon error:
set -e

mkdir $PDBNAME
DIR=$PDBNAME/mdrun
mkdir $DIR
cp $PDBNAME.pdb $PDBNAME/mdrun/pep.pdb

# create gmx files
cp -r setup_template/runfiles_initial/* $DIR

# check if forcefield is valid
if [ $FORCEFIELD != "amber99sbildn" ] && [ $FORCEFIELD != "grappa" ]; then
    echo "forcefield must be amber99sbildn or grappa but is $FORCEFIELD"
    exit 1
fi


source ~/.bashrc

conda activate $GRAPPA_ENV

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

echo "Done!"
popd