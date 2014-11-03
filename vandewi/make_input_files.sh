#!/bin/bash

for j in `ls fevrier2014/| sed s/"\/"//`;
do
    k=$(echo $j| awk -F_ '{print $2"_"$4}');
    echo ===========================
    echo $j $k;
    [ -e $k.in ] && rm -f ${k}.in
    touch ${k}.in
    for i in `ls fevrier2014/${j}/${k}_????.dat fevrier2014/${j}/${k}_???_???.dat`;
    do
	s=$(echo $i| awk -F/ '{print $3}' | sed s/.dat//| awk -F_ '{print $2}');
	if [[ $s == "pi0pi0" ]];
	then
	    t="#bar{p}p#rightarrow#pi^{0}#pi^{0}"
	elif [[ $s == "pi0eta" ]];
	then
	    t="#bar{p}p#rightarrow#pi^{0}#eta"
	elif [[ $s == "pi0gamma" ]];
	then
	    t="#bar{p}p#rightarrow#pi^{0}#gamma"
	elif [[ $s == "etaeta" ]];
	then
	    t="#bar{p}p#rightarrow#eta#eta"
	else
	    t="#bar{p}p#rightarrow#pi^{+}#pi{-}"
	fi
	e=$(echo $i| awk -F/ '{print $3}' | sed s/.dat//| awk -F_ '{print $3}');
	l=$(head -1 $i | wc -w);
	echo $k $t $i $l $e >> ${k}.in ;
    done
    echo

done
