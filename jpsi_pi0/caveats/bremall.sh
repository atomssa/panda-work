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

arg=$1

if [[ $arg == 17 || $arg == 2 || $arg == 8 ]]; then
    [[ $arg == 17 ]] && q=esim || q=psim
    #[[ $arg == 17 ]] && n=9 || n=20
    n=20
    for i in `seq 0 $n`;
    do
	export ODIR=/vol0/panda/work/jpsi_pi0/grid.out/${q}_oct14_binsong_config${arg}/runall.$i
	echo $ODIR
	tag_start_t $i
	root -b -q -w -l brem.C
	tag_end_t $i
    done
else
    for i in `seq 0 21`;
    do
	for j in `seq 0 4`;
	do
	    export ODIR=/vol0/panda/work/jpsi_pi0/grid.out/psim_oct14_binsong_configs/runall.$i.$j
	    echo $ODIR
	    tag_start_t $i
	    root -l -b -q -w brem.C
	    tag_end_t $i
	done
    done
fi
