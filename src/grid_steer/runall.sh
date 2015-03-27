#!/bin/bash

. /vol0/panda/svn/pandaroot-jan14/build/config.sh

nev=10000

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
    #echo "Total $1: $(date -u -d "0 $end seconds - $start seconds" +"%H:%M:%S")"
    echo "Total $1: $(date -d @$(($(date +%s)-$start)) +"%H:%M:%S hours")"
    echo "-----------------------------------------------"
    echo
}

tag_start_t sim
root -l -b -q -w sim_complete.C\($nev\) >> sim.log 2>&1
tag_end_t sim

tag_start_t digi
root -l -b -q -w digi_complete.C >> digi.log 2>&1
tag_end_t digi

tag_start_t sim
root -l -b -q -w reco_complete.C >> reco.log 2>&1
tag_end_t digi

tag_start_t pid
root -l -b -q -w pid_complete.C >> pid.log 2>&1
tag_end_t pid

