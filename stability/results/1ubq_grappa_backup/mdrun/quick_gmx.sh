#!/usr/bin/env bash

set -e

# assumes that there is a file pep.gro and pep.top in the current directory
# equilibrates in a short md, then creates files md.tpr and md.gro

gmx editconf -f pep.gro -o pep_box.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp pep_box.gro -p pep.top -o pep_solv.gro

gmx grompp -f ions.mdp -c pep_solv.gro -p pep.top -o pep_out_genion.tpr
echo "SOL" | gmx genion -s pep_out_genion.tpr -p pep.top -o pep_out_ion.gro -neutral

gmx grompp -f minim.mdp -c pep_out_ion.gro -p pep.top -o pep_out_min.tpr
gmx mdrun -deffnm pep_out_min -v

gmx grompp -f nvt.mdp -c pep_out_min.gro -p pep.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -p pep.top -o npt.tpr
gmx mdrun -v -deffnm npt

gmx grompp -f md.mdp -c npt.gro -t npt.trr -o md.tpr -p pep.top
gmx mdrun -deffnm md -v

echo "q" | gmx make_ndx -f npt.gro
