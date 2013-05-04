#!/bin/bash

for sample in `cat sampleList20`
do
    ls $sample/*_1.filt.fastq.gz
    numRuns=`ls $sample/*_1.filt.fastq.gz | wc -l`
    
    echo numRuns: $numRuns
    
    if [ $numRuns -eq 1 ]
    then
        echo Only one run... renaming
	strand1=`ls $sample/*_1.filt.fastq.gz`
	mv $strand1 $sample/${sample}_1.fq.gz
    else
	if [ $numRuns -gt 1 ]
        then
	    echo Multiple runs... concatenating files
	    zcat $sample/*_1.filt.fastq.gz >
	fi
    fi
    
    
    
done