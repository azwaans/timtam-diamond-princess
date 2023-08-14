#!/usr/bin/env bash
# -*- mode:sh; -*-
#
# Run MCMC
# ========
#
# This script runs a specified number of chains in parallel using ant
# and mcmc2.7.3
#
# Usage:
# ------
#
# ./script_name NUM_CHAINS CHAIN_LEN
#

# check the number of input arguments
if [[ $# -ne 2 ]]; then
    echo "Invalid number of arguments."
    echo "Usage: $0 NUM_CHAINS CHAIN_LEN"
    exit 1
fi

if [ ! -f ./data/diamond.fasta ]; then
    echo "The file ./data/diamond.fasta does not exist."
    echo "Please contact the author to obtain the file."
    exit 1
fi

# check if the arguments are valid integers
if ! [[ $1 =~ ^[0-9]+$ ]] || ! [[ $2 =~ ^[0-9]+$ ]]; then
    echo "Both arguments should be valid integers."
    exit 1
fi

# ensure the required output directory, out/log-files, exists
if [ ! -d ./out/log-files ]; then
	mkdir -p ./out/log-files
fi

NUM_CHAINS=$1
CHAIN_LEN=$2

BEAST_XML=xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml

IX=1
# use a random seed number 
# ID=$RANDOM
ID=15607
PIDS=() # Array to hold process IDs

while [ $IX -le $NUM_CHAINS ]
do
    ant -DstateFile=out/tmp-timtam-chain-num-"$IX".xml.state -Dseed=$ID  -DchainLength=$CHAIN_LEN -DbeastXML=$BEAST_XML mcmc2.7.3 & PIDS+=($!)
    ((IX++))
    ((ID++))
done

# Wait for all processes to complete
wait "${PIDS[@]}"

echo "All processes completed."
