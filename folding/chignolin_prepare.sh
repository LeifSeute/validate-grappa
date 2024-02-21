#! bin/bash

set -e

FASTA="chignolin.fasta"
PEPGEN_ENV={$1:-"pepgen"}
GRAPPA_ENV={$2:-"casc_grappa_gromacs"}

bash prepare.sh $FASTA $PEPGEN_ENV $GRAPPA_ENV amber99sbildn > prepare_chignolin.log

bash prepare.sh $FASTA $PEPGEN_ENV $GRAPPA_ENV grappa > prepare_chignolin.log