#!/bin/bash

set -e
# try for 4 hours
target_time=240

make

data_root="$( cd "$(dirname "$0")" ; pwd -P )"/..

PLANS=( "CurrentRep.txt" "NewRep.txt" )
INPUTS=( "cong12.csv" "cong16.csv" "pres08.csv" "pres12.csv" "pres16.csv" "sen10.csv" "sen16.csv")
RESULTS=( "Feb_8_Wes_units.csv" "Feb_11_Wes_units.csv" )
for plan in "${PLANS[@]}"; do
    for input in "${INPUTS[@]}"; do
        for result in "${RESULTS[@]}"; do
            [[ "$plan" == "CurrentRep.txt" ]] && l1_value=158 || l1_value=69
            folder="${plan}-${input}_${result}"
            echo $data_root/$folder
            mkdir -p $data_root/$folder
            pushd $data_root/$folder
            cp ../chain ./
            ./chain -f $data_root/$plan -n 30 -1 $l1_value -M -p 0.1 -H -E -S \
                    --filename_wes_units $data_root/bill_plans/$result \
                    --filename_election_results $data_root/$input --target_time $target_time \
                    &> out_err.log || echo "errored: $folder"
            popd
        done
    done
done
