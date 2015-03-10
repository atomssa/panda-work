#!/bin/bash

echo "utils.sh setting up "

function check_path {
    echo -n "Checking existence of $1... "
    if [[ ! -e $1 ]]; then
	echo "Doesn't exist... Quitting"
	return 1;
    else
	echo "OK!"
	return 0;
    fi
}

export HN=$(hostname)
# if using a differnt machine, change ipnphen01 to its host name!!!
# except if its on the grid. The else clause assumes that any other
# hostnae means stuff is being compiled from the grid workers
if [[ $HN == "ipnphen01" ]]; then
    export SOFT_DIR=/vol0/panda/svn
    export EXT_DIR=$SOFT_DIR/$EXT_VER
    export SIMPATH=$EXT_DIR/install
    export PATCH_DIR=/vol0/panda/work/src/patch-$PR_VER
    export NCORE=16
else
    export SOFT_DIR=/nfs1/panda/ermias/soft
    export EXT_DIR=$SOFT_DIR/$EXT_VER
    export SIMPATH=$EXT_DIR/install
    export PATCH_DIR=$SOFTDIR/panda-work/src/patch-$PR_VER
    export NCORE=4
fi

export PNDROOT_DIR=$SOFT_DIR/$EXT_VER/pandaroot/$PR_VER

echo "utils.sh setting up done "
