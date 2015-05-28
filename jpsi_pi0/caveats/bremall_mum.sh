#!/bin/bash

#. /vol0/panda/svn/pandaroot-oct14/build/config.sh
#. /vol0/panda/svn/apr13/pandaroot/oct14/build/config.sh
. /vol0/panda/svn/jul14p3/pandaroot/oct14/build/config.sh

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
    echo "Total $1: $(date -u -d "0 $end seconds - $start seconds" +"%H:%M:%S")"
    echo "-----------------------------------------------"
    echo
}

for i in `seq 0 24`;
do
    [[ $# == 1 ]] && beg=$1 || beg=0
    [[ $# == 1 ]] && end=$1 || end=5
    [[ $# == 2 ]] && beg=$1 && end=$2
    for j in `seq $beg $end`;
    do
	export ODIR=/vol0/panda/work/jpsi_pi0/grid.out/brem/mum_oct14_binsong_configs/runall.$i.$j
	echo $ODIR
	tag_start_t $i
	root -l -b -q -w brem.C
	tag_end_t $i
    done
done

