#!/bin/bash

pushd ..
NAMES=( "chignolin" "trp_cage" "trp_cage_from_amber") # "w_domain" "protein_g" )

for NAME in "${NAMES[@]}"
do
    echo "Evaluating $NAME"
    bash evaluate.sh $NAME
done

popd