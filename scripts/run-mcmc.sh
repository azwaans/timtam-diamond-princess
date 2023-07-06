#!/usr/bin/env bash
# -*- mode:sh; -*-

# This script runs a specified number of chains in parallel using ant and mcmc2.7.3

# Usage:
# ./script_name NUM_CHAINS CHAIN_LEN

# check the number of input arguments
if [[ $# -ne 2 ]]; then
    echo "Invalid number of arguments."
    echo "Usage: $0 NUM_CHAINS CHAIN_LEN"
    exit 1
fi

# check if the arguments are valid integers
if ! [[ $1 =~ ^[0-9]+$ ]] || ! [[ $2 =~ ^[0-9]+$ ]]; then
    echo "Both arguments should be valid integers."
    exit 1
fi

NUM_CHAINS=$1
CHAIN_LEN=$2

BEAST_XML=xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml

IX=1687361450300
PIDS=() # Array to hold process IDs

while [ $IX -le $NUM_CHAINS ]
do
    ant -DstateFile=out/tmp-timtam-chain-num-"$IX".xml.state -Dseed=$IX -DchainLength=$CHAIN_LEN -DbeastXML=$BEAST_XML mcmc2.7.3 & PIDS+=($!)
    ((IX++))
done

# Wait for all processes to complete
wait "${PIDS[@]}"

echo "All processes completed."
