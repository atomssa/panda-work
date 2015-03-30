#!/bin/bash

. /vol0/panda/svn/pandaroot-oct14/build/config.sh

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
    echo $i: $wdir
    wdir=/vol0/panda/work/jpsi_pi0/grid.out/psim_oct14_binsong_config2/runall.$i
    tag_start_t $i
    root -b -q -w -l brem.C"(\"$wdir\")"; 
    tag_end_t $i
done

for i in `seq 0 20`; 
do 
    echo $i: $wdir
    wdir=/vol0/panda/work/jpsi_pi0/grid.out/psim_oct14_binsong_config8/runall.$i
    tag_start_t $i
    root -b -q -w -l brem.C"(\"$wdir\")"; 
    tag_end_t $i
done

for i in `seq 0 21`;
do
    for j in `seq 0 5`;
    do 
	wdir=/vol0/panda/work/jpsi_pi0/grid.out/psim_oct14_binsong_configs/runall.$i.$j
	echo $i.$j: $wdir
	tag_start_t $i
	root -l -b -q -w brem.C"(\"$wdir\")"
	tag_end_t $i
    done
done
