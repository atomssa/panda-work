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

_type=$1
_iplab=$2
_brem=1

if [[ $_type == 0 ]]; then
    listf=lists/pi0pipm_iplab${_iplab}.list
else
    listf=lists/pi0jpsi_iplab${_iplab}.short.list
fi

echo listf=$listf

for ifile in `more $listf`; do
    tag=type${_type}_plab${_plab}_brem${_brem}_${ifile}
    logfile=log/anav2/${tag}.log
    echo "doing analysis for $tag, logfile: $logfile"
    nevt=0
    if [[ $_type == 1 ]]; then
	if [[ $_iplab == 0 && $ifile == 0 ]]; then
	    nevt=8040
	    echo nevt=$nevt
	elif [[ $_iplab == 1 && $ifile == 3 ]]; then
	    nevt=3013
	    echo nevt=$nevt
	elif [[ $_iplab == 2 && $ifile == 0 ]]; then
	    nevt=4425
	    echo nevt=$nevt
	else
	    echo "nevt= 0"
	fi
    else
	echo nevt=$nevt
    fi
    tag_start_t $tag
    root -l -b -q -w macros/anav2.C\($_iplab,$_type,$_brem,$ifile,$nevt\) >> $logfile 2>&1
    tag_end_t $tag
done

#
