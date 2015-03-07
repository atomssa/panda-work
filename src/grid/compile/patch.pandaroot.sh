#!/bin/bash

if [[ $# != 2 ]]; then
    echo usage: ./patch.pandaroot.sh EXT_VER PR_VER
    echo example1: ./patch.pandaroot.sh apr13 oct14
    echo example2: ./patch.pandaroot.sh jul14p3 trunk-26841
    return
else
    export EXT_VER=$1
    export PR_VER=$2
    echo "Patching EXT_VER=$EXT_VER and PR_VER=$PR_VER"
fi

function check_path {
    echo "Checking existence of $1"
    if [[ ! -e $1 ]]; then
	echo "Doesn't exist... Quitting"
	return 1;
    else
	return 0;
    fi
}

#input sanitation
export SOFTDIR=/nfs1/panda/ermias/soft
export PNDROOT_DIR=$SOFTDIR/$EXT_VER/pandaroot/$PR_VER
export PATCH_DIR=$SOFTDIR/patch-$PR_VER

if ! check_path $PNDROOT_DIR; then return; fi
if ! check_path $PATCH_DIR; then return; fi

export SIMPATH=$SOFTDIR/$EXT_VER/install
if ! check_path $SIMPATH; then return; fi

if [[ -e  $SOFTDIR/$EXT_VER/FairRoot ]]; then
    export FAIRROOTPATH=$SOFTDIR/$EXT_VER/FairRoot/install
    if ! check_path $FAIRROOTPATH; then return; fi
fi

cd $PATCH_DIR
for SRC_FILE in `find . -type f`;
do
    echo "############################################################################################################"
    #echo "Checking patch candidate $SRC_FILE"
    if [[ -e $PNDROOT_DIR/$SRC_FILE ]]; then
	echo "$SRC_FILE found in $PNDROOT_DIR"
	diff_cnt=$(diff -b $PNDROOT_DIR/$SRC_FILE $SRC_FILE |wc -l)
	if [[ $diff_cnt > 0 ]]; then
	    echo "Patch file $SRC_FILE different from $PNDROOT_DIR version. Patching.."
	    nbkp=$(ls $PNDROOT_DIR/$SRC_FILE.org.* 2>/dev/null | wc -l)
	    diff -b $PNDROOT_DIR/$SRC_FILE $PATCH_DIR/$SRC_FILE
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

echo "Done patching files. Now to building"

# build the patch
cd $PNDROOT_DIR/build
make -j 4

# reset write permissions
cd $SOFTDIR/$EXT_VER
chmod -f g+w .
chmod -Rf g+w pandaroot
