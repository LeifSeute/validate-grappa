#!/bin/bash


# (the 6 1 flags are to select the amber99sbildn forcefield and water model)
rm folding/pep.top
printf "6\n1\n "|gmx pdb2gmx -f equilibration/pep.pdb -o pep_new.gro -p folding/pep.top -ignh
