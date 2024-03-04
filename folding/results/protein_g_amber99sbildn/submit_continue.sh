#!/bin/bash

#SBATCH -N 1                    
#SBATCH -p cascade.p              
#SBATCH --gres=gpu:1
#SBATCH --mincpus=20 
#SBATCH -t 24:00:00
#SBATCH --job-name="md"   
#SBATCH --output=logfiles/job.o%j
#SBATCH --error=logfiles/job.e%j

## call with sbatch -J <job name> submit_continue.sh <rundir> <nr> <box len>
## make sure directory $1 contains pep.gro, pep.top
## and mdp/ contains mdp files

echo "Job Name:"
echo "$SLURM_JOB_NAME"

# throw upon error
set -e

rundir='folding'

TASK_ID=$1
box_side=$2 # in nm

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
# echo "SOL" | gmx genion -s pep_genion.tpr -p pep.top -o pep_ion.gro -neutral 
echo "SOL" | gmx genion -s pep_genion.tpr -p pep.top -o pep_ion.gro -pname NA -nname CL -neutral -conc 0.1

gmx grompp -f minim.mdp -c pep_ion.gro -p pep.top -o pep_min.tpr
gmx mdrun -deffnm pep_min 

gmx grompp -f nvt.mdp -c pep_min.gro -p pep.top -o nvt.tpr 
gmx mdrun -deffnm nvt 

gmx grompp -f npt.mdp -c nvt.gro -p pep.top -o npt.tpr 
gmx mdrun -deffnm npt 

gmx grompp -f md.mdp -c npt.gro -p pep.top -o md.tpr -maxwarn 2

gmx mdrun -deffnm md -maxh $CYCLE
cd ../..
fi

# printf "1\n0\n " | gmx trjconv -f md.xtc -s md.tpr -o md_all.xtc -center -pbc mol 
# printf "1\n1\n " | gmx trjconv -f md.xtc -s md.tpr -o md_protein.xtc -center -pbc mol 

END=$(date +"%s")
echo "$(((END-START))) seconds ran"
echo "$(((END-START)/3600)) full hours ran"
let "CYCLE--"
if [ $(((END-START)/3600)) -lt $CYCLE ]
then
        echo "last cycle was just $(((END-START)/3600))h long and therefore finito"
else
        sbatch -J $SLURM_JOB_NAME submit_continue.sh ${TASK_ID} ${box_side}
fi
