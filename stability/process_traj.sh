
PDBNAME=$1

set -e

DIR1="results/${PDBNAME}_grappa/mdrun"
DIR2="results/${PDBNAME}_amber99sbildn/mdrun"

bash ../folding/remove_pbc.sh $DIR1 Protein
bash ../folding/remove_pbc.sh $DIR2 Protein