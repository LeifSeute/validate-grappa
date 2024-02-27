#!/bin/bash

box_side=$1 # in nm

# throw upon error
set -e

rundir='folding'


echo run_dir: $rundir
echo TASK_ID: $TASK_ID
echo box_side: $box_side

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

CYCLE=24

if [ "${rundir: -1}" == "/" ]
then
	rundir=${rundir::-1}
fi

if [ -e mds/${rundir}_${TASK_ID}/md.tpr ];then
cd mds/${rundir}_${TASK_ID}
gmx mdrun -deffnm md -maxh $CYCLE -cpi md.cpt
cd ../..
else
cp -r $rundir mds/${rundir}_${TASK_ID}
cp mdps/* mds/${rundir}_${TASK_ID}
cd mds/${rundir}_${TASK_ID}

gmx editconf -f pep.gro -o pep_box.gro -c -box $box_side $box_side $box_side -bt cubic
gmx solvate -cp pep_box.gro -p pep.top -o pep_solv.gro 

gmx grompp -f ions.mdp -c pep_solv.gro -p pep.top -o pep_genion.tpr 
echo "SOL" | gmx genion -s pep_genion.tpr -p pep.top -o pep_ion.gro -neutral 

gmx grompp -f minim.mdp -c pep_ion.gro -p pep.top -o pep_min.tpr
gmx mdrun -deffnm pep_min 

gmx grompp -f nvt.mdp -c pep_min.gro -p pep.top -o nvt.tpr 
gmx mdrun -deffnm nvt 

gmx grompp -f npt.mdp -c nvt.gro -p pep.top -o npt.tpr 
gmx mdrun -deffnm npt 

gmx grompp -f md.mdp -c npt.gro -p pep.top -o md.tpr -maxwarn 2

gmx mdrun -deffnm md -maxh 1
cd ../..
fi

# printf "1\n0\n " | gmx trjconv -f md.xtc -s md.tpr -o md_all.xtc -center -pbc mol 
# printf "1\n1\n " | gmx trjconv -f md.xtc -s md.tpr -o md_protein.xtc -center -pbc mol 

