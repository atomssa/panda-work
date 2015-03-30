#!/bin/bash

. /vol0/panda/svn/apr13/pandaroot/oct14/build/config.sh

tt=$1
bb=$2

echo
echo

function tag_start_t {
    echo "-----------------------------------------------"
    echo -n "$1 start: "
    date
    start=$(date +%s)
}

function tag_end_t {
    echo -n "$1 end: "
    date
    end=$(date +%s)
    echo "Total $1: $(date -d @$(($(date +%s)-$start)) +"%H:%M:%S hours")"
    echo "-----------------------------------------------"
    echo
}


for i in `seq 0 20`;
do
    echo type:$tt brem:$bb fid:$i
    tag_start_t $i
    #root -l -b -q -w Brem_corr_ele.C"(\"$wdir\")"  #\(\"$wdir\"\)
    root -l -b -q -w anav2.C"(0, $tt, $bb, $i)"
    tag_end_t $i
done
