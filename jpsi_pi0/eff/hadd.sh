#!/bin/bash

. /vol0/panda/svn/jul14p3/pandaroot/oct14/build/config.sh

if [[ $1 != "elec" && $1 != "posit" && $1 != "pip" && $1 != "pim" && $1 != "fin" ]]; then
    echo "unrecognized option $1"

elif [[ $1 == "fin" ]]; then
    echo  logs/hadd.pi.a.log
    hadd hadd_out/hadd.pi.a.root hadd_out/hadd.pip.a.root hadd_out/hadd.pim.a.root >> logs/hadd.pi.a.log 2>&1
    echo  logs/hadd.pi.b.log
    hadd hadd_out/hadd.pi.b.root hadd_out/hadd.pip.b.root hadd_out/hadd.pim.b.root >> logs/hadd.pi.b.log 2>&1
    echo  logs/hadd.pi.c.log
    hadd hadd_out/hadd.pi.c.root hadd_out/hadd.pip.c.root hadd_out/hadd.pim.c.root >> logs/hadd.pi.c.log 2>&1
    echo  logs/hadd.pi.log
    hadd hadd_out/hadd.pi.root hadd_out/hadd.pi.?.root >> logs/hadd.pi.log 2>&1
    echo  logs/hadd.e.a.log
    hadd hadd_out/hadd.e.a.root hadd_out/hadd.elec.a.root hadd_out/hadd.posit.a.root >> logs/hadd.e.a.log 2>&1
    echo  logs/hadd.e.b.log
    hadd hadd_out/hadd.e.b.root hadd_out/hadd.elec.b.root hadd_out/hadd.posit.b.root >> logs/hadd.e.b.log 2>&1
    echo  logs/hadd.e.c.log
    hadd hadd_out/hadd.e.c.root hadd_out/hadd.elec.c.root hadd_out/hadd.posit.c.root >> logs/hadd.e.c.log 2>&1
    echo  logs/hadd.e.log
    hadd hadd_out/hadd.e.root hadd_out/hadd.e.?.root >> logs/hadd.e.log 2>&1
else

    for i in `ls lists/hadd.${1}.*.list`;
    do
	k=$(echo $i| sed s/lists/hadd_out/ | sed s/\.list/\.root/);
	l=$(echo $i| sed s/lists/logs/ | sed s/\.list/\.log/);
	echo $i "  " $k "  " $l;
	hadd $k $(< $i) >> $l 2>&1;
    done

    for i in `a b c`;
    do
	if [[ $(ls hadd_out/hadd.${1}.${i}*.root| wc -l| awk '{print $1}') > 1 ]]; then
	    echo  logs/hadd.${1}.${i}.log
	    hadd hadd_out/hadd.${1}.${i}.root hadd_out/hadd.${1}.${i}.?.root >> logs/hadd.${1}.${i}.log 2>&1
	fi
    done

fi


