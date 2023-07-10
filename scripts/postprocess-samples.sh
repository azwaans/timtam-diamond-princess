#!/usr/bin/env bash
# -*- mode:sh; -*-


XML_FILE=./xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml
LOG_PARAMS=./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.1.log
LOG_TREES=./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories_diamond.1.trees
OUTPUT_CSV=./out/prevalence-estimate-HKY-weekly-histories.csv

Rscript ./R/combine-ltt-and-histories.R -v -x $XML_FILE -t $LOG_TREES -p $LOG_PARAMS -o $OUTPUT_CSV
Rscript R/postprocessing-part-4.R
