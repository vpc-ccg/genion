#!/bin/bash
set -e
SCRIPT_PATH="$(dirname `which $0`)/scripts"
echo $SCRIPT_PATH
SNAKEFILE=makeref.snake
#Checking if tools exist
command -v snakemake && echo "... Exists :)" >&2 || { echo snakemake missing! >&2; FAILED=1; }

if [ -z $FAILED ]; then 
    snakemake -s "${SCRIPT_PATH}/${SNAKEFILE}" --config SCPA=${SCRIPT_PATH}  "$@";
fi
