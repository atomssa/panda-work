#!/bin/bash

# here is how to add in one go
#hadd hadd/elec.a.root $(< lists/hadd.elec.a.list)
# or all of electrons:

#for i in `ls lists/hadd.elec.*.list`; do echo $i; k=$(echo $i| sed s/lists/hadd.hadd_out/ | sed s/list/root/ | sed s/hadd.//); echo $k; hadd $k $(< $i); done

function create_list {

    tag=$1
    grp=$2

    if [[ $grp != "a" && $grp != "b" && $grp != "c" ]]; then
       echo "grp $grp unrecongized "
       return
    fi

    echo tag=$tag, grp=$grp

    if [[ $tag == "elec" ]]; then

	if [[ $grp == "a" ]]; then
	    ls /projet/panda/Ermias/tda/elec_flat/runall.?/effhists.root > lists/hadd.elec.a.list
	    ls /projet/panda/Ermias/tda/elec_flat/runall.??/effhists.root >> lists/hadd.elec.a.list
	    ls /projet/panda/Ermias/tda/elec_flat/runall.100/effhists.root >> lists/hadd.elec.a.list
	elif [[ $grp == "b" ]]; then
	    ls /projet/panda/Ermias/tda/elec_flat/runall.10[1-9]/effhists.root > lists/hadd.elec.b.list
	    ls /projet/panda/Ermias/tda/elec_flat/runall.1[1-4]?/effhists.root >> lists/hadd.elec.b.list
	elif [[ $grp == "c" ]]; then
	    ls /projet/panda/Ermias/jacek/elec/runall.*/effhists.root > lists/hadd.elec.c.list
	else
	    echo unrecognized grp $grp
	fi

    elif [[ $tag == "posit" ]]; then
	echo tag=$tag

	if [[ $grp == "a" ]]; then
	    ls /projet/panda/Ermias/tda/posit_flat/runall.?/effhists.root > lists/hadd.posit.a.list
	    ls /projet/panda/Ermias/tda/posit_flat/runall.??/effhists.root >> lists/hadd.posit.a.list
	    ls /projet/panda/Ermias/tda/posit_flat/runall.100/effhists.root >> lists/hadd.posit.a.list
	elif [[ $grp == "b" ]]; then
	    ls /projet/panda/Ermias/tda/posit_flat/runall.10[1-9]/effhists.root > lists/hadd.posit.b.list
	    ls /projet/panda/Ermias/tda/posit_flat/runall.1[1-4]?/effhists.root >> lists/hadd.posit.b.list
	elif [[ $grp == "c" ]]; then
	    ls /projet/panda/Ermias/jacek/posit/runall.*/effhists.root > lists/hadd.posit.c.list
	else
	    echo unrecognized grp $grp
	fi
    elif [[ $tag == "pip" ]]; then
	echo tag=$tag

	if [[ $grp == "a" ]]; then
	    ls /projet/panda/Ermias/tda/pip_flat/runall.1??/effhists.root > lists/hadd.pip.a.0.list
	    ls /vol0/panda/work/jpsi_pi0/grid.out/tda/pip_flat/runall.2??/effhists.root > lists/hadd.pip.a.1.list
	    ls /projet/panda/Ermias/tda/pip_flat/runall.3??/effhists.root > lists/hadd.pip.a.2.list
	    ls /projet/panda/Ermias/tda/pip_flat/runall.4??/effhists.root > lists/hadd.pip.a.3.list
	elif [[ $grp == "b" ]]; then
	    ls /projet/panda/Ermias/tda/pip_flat/runall.5??/effhists.root > lists/hadd.pip.b.0.list
	    ls /projet/panda/Ermias/tda/pip_flat/runall.6??/effhists.root > lists/hadd.pip.b.1.list
	elif [[ $grp == "c" ]]; then
	    ls /projet/panda/Ermias/jacek/piplus/runall.*/effhists.root > tmp
	    split -d -100 tmp tmp.
	    for i in `ls tmp.*`; do mv $i lists/hadd.pip.c.$(expr 0 + $(echo $i| awk -F. '{print $2}')).list; done
	    rm -vf tmp
	else
	    echo unrecognized grp $grp
	fi

    elif [[ $tag == "pim" ]]; then
	echo tag=$tag

	if [[ $grp == "a" ]]; then
	    ls /projet/panda/Ermias/tda/pim_flat/runall.1??/effhists.root > lists/hadd.pim.a.0.list
	    ls /vol0/panda/work/jpsi_pi0/grid.out/tda/pim_flat/runall.2??/effhists.root > lists/hadd.pim.a.1.list
	    ls /projet/panda/Ermias/tda/pim_flat/runall.3??/effhists.root > lists/hadd.pim.a.2.list
	    ls /projet/panda/Ermias/tda/pim_flat/runall.4??/effhists.root > lists/hadd.pim.a.3.list
	elif [[ $grp == "b" ]]; then
	    ls /projet/panda/Ermias/tda/pim_flat/runall.5??/effhists.root > lists/hadd.pim.b.0.list
	    ls /projet/panda/Ermias/tda/pim_flat/runall.6??/effhists.root > lists/hadd.pim.b.1.list
	elif [[ $grp == "c" ]]; then
	    ls /projet/panda/Ermias/jacek/piminus/runall.*/effhists.root > tmp
	    split -d -100 tmp tmp.
	    for i in `ls tmp.*`; do mv $i lists/hadd.pim.c.$(expr 0 + $(echo $i| awk -F. '{print $2}')).list; done
	    rm -vf tmp
	else
	    echo unrecognized grp $grp
	fi
    else
	echo unrecongnized tag $tag
    fi
}

if [[ $1 == "all" ]]; then
    create_list elec a
    create_list elec b
    create_list elec c
    create_list posit a
    create_list posit b
    create_list posit c
    create_list pip a
    create_list pip b
    create_list pip c
    create_list pim a
    create_list pim b
    create_list pim c
elif [[ $2 == "all" ]]; then
     create_list $1 a
     create_list $1 b
     create_list $1 c
else
    create_list $1 $2
fi
