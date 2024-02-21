#! bin/bash

set -e

FASTA="chignolin.fasta"
PEPGEN_ENV=pepgen
GRAPPA_ENV=casc_grappa_gromacs

bash prepare.sh $FASTA $PEPGEN_ENV $GRAPPA_ENV amber99sbildn > prepare_chignolin.log

bash prepare.sh $FASTA $PEPGEN_ENV $GRAPPA_ENV grappa > prepare_chignolin.log