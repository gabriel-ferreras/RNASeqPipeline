#! /bin/bash

if [ $# -ne 1 ]
then
    echo ""
    echo "    Usage: bash RNASeqPipeline.sh <params_file>"
    echo ""
    echo "    params_file: imput file with the parameters."
    echo ""
    exit
fi

PARAMS=$1

EXP= $(grep experiment_name: $PARAMS | awk '{ print $2 }')
echo "experiment name = "$EXP