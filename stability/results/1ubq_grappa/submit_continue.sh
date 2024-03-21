#!/bin/bash

#SBATCH -N 1                    
#SBATCH -p cascade.p              
#SBATCH --gres=gpu:1
#SBATCH --mincpus=20 
#SBATCH -t 24:00:00
#SBATCH --job-name="md"   
#SBATCH --output=logfiles/job.o%j
#SBATCH --error=logfiles/job.e%j


echo "Job Name:"
echo "$SLURM_JOB_NAME"

# throw upon error
set -e

rundir='mdrun'

echo run_dir: $rundir

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5

CYCLE=24

if [ "${rundir: -1}" == "/" ]
then
	rundir=${rundir::-1}
fi

if [ -e ${rundir}/md.tpr ];then
cd ${rundir}
gmx mdrun -deffnm md -maxh $CYCLE -cpi md.cpt
cd ../..
else
echo "no md.tpr found in ${rundir}"
exit 1
fi

END=$(date +"%s")
echo "$(((END-START))) seconds ran"
echo "$(((END-START)/3600)) full hours ran"
let "CYCLE--"
if [ $(((END-START)/3600)) -lt $CYCLE ]
then
        echo "last cycle was just $(((END-START)/3600))h long and therefore finito"
else
        sbatch -J $SLURM_JOB_NAME submit_continue.sh
fi
