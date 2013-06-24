#!/bin/bash

if [ $# -lt 1 ]
then
    echo usage $0 [sampleList]
    exit 1
fi

ROOTDIR=`pwd`
SAMPLELIST=$1

for SAMPLE in `cat $SAMPLELIST`
do
    cd $SAMPLE

    LOGFILE=qc.log
    echo -ne "" > $LOGFILE

    for BAM in `ls *.bam`
    do
	echo -ne "checkBam $BAM ... " >> $LOGFILE
	echo `~/code/seque/checkBam.sh $BAM` >> $LOGFILE
    done

    cd $ROOTDIR


done