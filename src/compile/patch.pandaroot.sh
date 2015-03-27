#!/bin/bash

if [[ $# < 2 ]]; then
    echo "usage: ./patch.pandaroot.sh EXT_VER PR_VER"
    echo "example1(local): ./patch.pandaroot.sh ext-apr13 oct14"
    echo "example1(grid): ./patch.pandaroot.sh apr13 oct14"
    echo "example2(grid): ./patch.pandaroot.sh jul14p3 trunk-26841"
    exit
else
    if [[ $1 == '-n' ]]; then
	DO_COMPILE=1
	export EXT_VER=$2
	export PR_VER=$3
    else
	DO_COMPILE=0
	export EXT_VER=$1
	export PR_VER=$2
    fi
    echo "Patching EXT_VER=$EXT_VER and PR_VER=$PR_VER"
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
UTILS_SCRIPT=$SCRIPT_DIR/utils.sh
if [ ! -e $UTILS_SCRIPT ]; then
    echo "utils.sh not found in the same directory as compile script. quitting"
    exit
fi
. $UTILS_SCRIPT

#input sanitation
if ! check_path $PNDROOT_DIR; then exit; fi
if ! check_path $PATCH_DIR; then exit; fi
if ! check_path $SIMPATH; then exit; fi
if [[ -e  $EXT_DIR/FairRoot ]]; then
    export FAIRROOTPATH=$EXT_DIR/FairRoot/install
    if ! check_path $FAIRROOTPATH; then exit; fi
fi

cd $PATCH_DIR
for SRC_FILE in `find . -type f`;
do
    echo "############################################################################################################"
    #echo "Checking patch candidate $SRC_FILE"
    if [[ -e $PNDROOT_DIR/$SRC_FILE ]]; then
	echo "$SRC_FILE found in $PNDROOT_DIR"
	diff_cnt=$(diff -b $PNDROOT_DIR/$SRC_FILE $SRC_FILE |wc -l)
	if [[ $diff_cnt -gt 0 ]]; then
	    echo "Patch file $SRC_FILE different from $PNDROOT_DIR version. Patching.."
	    nbkp=$(ls $PNDROOT_DIR/$SRC_FILE.org.* 2>/dev/null | wc -l| tr -d '[[:space:]]')
	    git diff -b $PNDROOT_DIR/$SRC_FILE $PATCH_DIR/$SRC_FILE
	    mv -vf $PNDROOT_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE.org.$nbkp
	    cp -vf $PATCH_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE
	else
	    echo "Patch file identical to $PNDROOT_DIR version. Doing nothing"
	fi
    else
	echo "$SRC_FILE not found in $PNDROOT_DIR"
	SRC_FILE_DIR=$(dirname $SRC_FILE)
	echo "Creating parent directory if it doesn't exist"
	[ -e $PNDROOT_DIR/$SRC_FILE_DIR ] || mkdir -vp $PNDROOT_DIR/$SRC_FILE_DIR
	cp -vf $PATCH_DIR/$SRC_FILE $PNDROOT_DIR/$SRC_FILE
    fi
done
echo "############################################################################################################"
echo "Done patching files."

if [[ $DO_COMPILE -eq 0 ]]; then
    echo "Now building ..."
    # build the patch
    cd $PNDROOT_DIR/build
    make -j 4
    if [[ $HN != "ipnphen01" && $HN != "rasalula" ]]; then
	# reset write permissions
	cd $SOFT_DIR/$EXT_VER
	chmod -f g+w .
	chmod -Rf g+w pandaroot
    fi
fi
