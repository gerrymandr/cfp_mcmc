#!/bin/bash

set -e

make
rm histograms/*
INPUTS=( "cong12.csv" "cong16.csv" "pres08.csv" "pres12.csv" "pres16.csv" "sen10.csv" "sen16.csv")
RESULTS=( "Feb_8_Wes_units.csv" "Feb_11_Wes_units.csv" )
for input in "${INPUTS[@]}"; do
    for result in "${RESULTS[@]}"; do
        echo $input $result
        ./chain -f CurrentRep.txt -n 21 -d 22 -1 68 -M -p 0.1 --histogram --filename_wes_units bill_plans/$result --filename_election_results $input | grep median_mean | tail -1 | awk '{ print $NF }'
        cp median_mean_hist.txt histograms/histogram-$result-$input
    done
done
