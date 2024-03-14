#!/bin/bash

pushd ..
NAMES=( "chignolin" "trp_cage") # "w_domain" "protein_g" )

for NAME in "${NAMES[@]}"
do
    echo "Evaluating $NAME"
    bash evaluate.sh $NAME
done

popd