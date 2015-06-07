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

batch=$1
if [[ $batch == 0 || $batch == 1 || $batch == 2 || $batch == 3 ]]; then
    if [[ $2 == "" || $3 == "" ]]; then
	s=0
	e=0
    else
	s=$2
	e=$3
    fi
else
    s=0
    e=0
fi

if [[ $e == 0 ]]; then
    for ifile in `more lists/batch.${batch}.list`;
    do
	tag=batch${batch}_file${ifile}
	logfile=logs/${tag}.log
	echo "doing analysis for $tag, logfile: $logfile"
	tag_start_t $tag
	root -l -b -q -w effhists.C\($batch,$ifile\) >> $logfile 2>&1
	tag_end_t $tag
    done
else
    for ifile in `seq $s $e`;
    do
	tag=batch${batch}_file${ifile}
	logfile=logs/${tag}.log
	echo "doing analysis for $tag, logfile: $logfile"
	tag_start_t $tag
	root -l -b -q -w effhists.C\($batch,$ifile\) >> $logfile 2>&1
	tag_end_t $tag
    done
fi

