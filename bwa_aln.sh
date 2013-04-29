#!/bin/bash

set -e

USAGE='\tusage: '$0' <ref> <fastq> <outname> [threads=16]'

if [ $# -lt 3 ] || [ $# -gt 4 ]
then
	echo -e "\n$USAGE\n"
	exit
fi

if [ $# -eq 3 ]
then
	threads=16
else
	threads=$4
fi

filepath=$2
filename=${2##*/}
basename=${filename%.*}
outname=$3
outbase=${outname%.*}
ref=$1

if [ ${outname##*.} != "bam" ]
then
        echo -e "\nOutput file must have .bam extension"
        echo -e "$USAGE\n"
        exit
fi


echo ""
echo Ref: $ref
echo Fastq: $filepath
echo Output: $outname
echo ""

# align with bwa
bwa aln -n 1 -l 15 -k 1 -O 20 -E 10 -t $threads $ref $filepath > $outbase.sai
bwa samse $ref $outbase.sai $filepath | samtools view -Sb - -o $outbase.unsorted.bam

rm $outbase.sai

#sort the file
java -Xmx4g -Djava.io.tmpdir=/mnt/thor_pool1/user_data/cc2qe/tmp -jar /shared/bin/picard-tools-1.40/SortSam.jar I=$outbase.unsorted.bam O=$outbase.bam SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

rm $outbase.unsorted.bam