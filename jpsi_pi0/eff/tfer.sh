#!/bin/bash

start=/projet/panda/Ermias/tda/

if [[ $1 == 0 ]]; then
    for i in `seq 500 699`;
    do
	echo $i;
	cd $start/pim_flat;
	scp -r atomssa@ipngrid02:/nfs1/scratch/ermias/tda/pim_flat/runall.${i} .;
	[ -e $start/pim_flat/runall.${i} ] && cd $start/pim_flat/runall.${i};
	[ -e output.tar ] && tar -xvf output.tar;
	[ -e pid_complete.root -a -e simparams.root ] && rm -vrf output.tar;
	cd $start/pip_flat;
	scp -r atomssa@ipngrid02:/nfs1/scratch/ermias/tda/pip_flat/runall.${i} .;
	[ -e $start/pip_flat/runall.${i} ] && cd $start/pip_flat/runall.${i};
	[ -e output.tar ] && tar -xvf output.tar;
	[ -e pid_complete.root -a -e simparams.root ] && rm -vrf output.tar;
    done
elif [[ $1 == 1 ]]; then
    for i in `seq 101 140`;
    do
	echo $i;
	cd $start/elec_flat;
	scp -r atomssa@ipngrid02:/nfs1/scratch/ermias/tda/elec_flat/runall.${i} .;
	[ -e $start/elec_flat/runall.${i} ] && cd $start/elec_flat/runall.${i};
	[ -e output.tar ] && tar -xvf output.tar;
	[ -e pid_complete.root -a -e simparams.root ] && rm -vrf output.tar;
	cd $start/posit_flat;
	scp -r atomssa@ipngrid02:/nfs1/scratch/ermias/tda/posit_flat/runall.${i} .;
	[ -e $start/posit_flat/runall.${i} ] && cd $start/posit_flat/runall.${i};
	[ -e output.tar ] && tar -xvf output.tar;
	[ -e pid_complete.root -a -e simparams.root ] && rm -vrf output.tar;
    done
else
    #to transfer already untarred files on ipngrid02
    _dir=$PWD
    _base=$(basename $_dir)
    echo transfering directory $_dir with basename $_base from $1 to $2

    for i in `seq $1 $2`;
    do
    	echo $i;
    	[ -e $_dir/runall.${i} ] || mkdir $_dir/runall.${i};
    	cd $_dir/runall.${i};
    	[ -e pid_complete.root ] || scp atomssa@ipngrid02:/nfs1/scratch/ermias/tda/$_base/runall.${i}/pid_complete.root .;
    	[ -e simparams.root ] || scp atomssa@ipngrid02:/nfs1/scratch/ermias/tda/$_base/runall.${i}/simparams.root . ;
    	cd $_dir
    done
fi