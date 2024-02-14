#!/usr/bin/env bash
# -*- mode:sh; -*-

SEED_NR=$(head -n 1 ./out/log-files/run_seed.txt)

echo "Postprocessing the outputs for seed : $SEED_NR"
 
XML_FILE=./xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml
LOG_PARAMS=./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.$SEED_NR.log
LOG_TREES=./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories-diamond.$SEED_NR.trees
OUTPUT_CSV=./out/prevalence-estimate-HKY-weekly-histories.csv

Rscript ./R/combine-ltt-and-histories.R -v -x $XML_FILE -t $LOG_TREES -p $LOG_PARAMS -o $OUTPUT_CSV

echo "Generating plots for seed: $SEED_NR"

OUTPUT_PNG_R=./out/manuscript/r0-estimates.png
OUTPUT_PNG_P=./out/manuscript/prevalence-estimates.png

Rscript ./R/postprocessing-part-4.R -v -x $XML_FILE -t $LOG_TREES -p $LOG_PARAMS -r $OUTPUT_PNG_R -n $OUTPUT_PNG_P -c $OUTPUT_CSV 

echo "All plots generated."