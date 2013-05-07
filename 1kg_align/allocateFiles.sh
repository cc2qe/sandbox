#!/bin/bash

if [ $# -ne 3 ]
then
    echo usage $0 [sampleList] [batchDirectory] [node]
    exit 1
fi

# Directory and data names
SAMPLELIST=$1
BATCHDIR=$2
ROOTDIR=/scratch/cc2qe/1kg/batch1
NODE=$3

# ---------------------
# STEP 1: Allocate the data to the local drive

# copy the files to the local drive
# Require a lot of memory for this so we don't have tons of jobs writing to drives at once

# make the working directory
MOVE_FILES_CMD="for SAMPLE in \`cat $SAMPLELIST\`
do
    WORKDIR=$ROOTDIR/\$SAMPLE &&
    SAMPLEDIR=$BATCHDIR/\$SAMPLE

    mkdir -p \$WORKDIR &&
    rsync -rv \$SAMPLEDIR/* \$WORKDIR
done"

#MOVE_FILES_CMD="echo MOVE_FILES_CMD"

# set a high priority on the move files so it will finish before anyone starts aligning
#MOVE_FILES_Q=`$QUICK_Q -m 32500mb -d $NODE -t 1 -n move_${NODE} -c " $MOVE_FILES_CMD " -q $QUEUE -p 100`


#for SAMPLE in `cat $SAMPLELIST`
#do
#    `~/code/sandbox/1kg_align/alignDriver.sh $SAMPLE $BATCHDIR/$SAMPLE $NODE $MOVE_FILES_Q`
#done