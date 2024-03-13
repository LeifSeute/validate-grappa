

GRAPPA_ENV=casc_grappa_gromacs
PDBNAME=1ubq
FORCEFIELD=grappa

module load GCC/9.3.0
module load CUDA/11.5.0

module use /hits/fast/mbm/broszms/software/modules/
module load gromacs-2023-cuda-11.5


bash run.sh $PDBNAME $GRAPPA_ENV $FORCEFIELD

PDBNAME=1rf7

bash run.sh $PDBNAME $GRAPPA_ENV $FORCEFIELD