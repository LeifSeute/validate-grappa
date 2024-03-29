#! bin/bash

set -e

cd ..

NAME="chignolin_2"

FASTA="references/$NAME.fasta"
PEPGEN_ENV=pepgen
GRAPPA_ENV=casc_grappa_gromacs

bash prepare.sh $FASTA $PEPGEN_ENV $GRAPPA_ENV amber99sbildn > prepare_$NAME.log