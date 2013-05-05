#!/bin/bash

if [ $# -ne 3 ]
then
    echo usage $0 [sample] [sampleDirectory] [node]
    exit 1
fi

# Directory and data names
SAMPLE=$1
SAMPLEDIR=$2
ROOTDIR=/scratch/cc2qe/1kg/batch1
WORKDIR=$ROOTDIR/$SAMPLE

# Annotations
REF=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/human_b37_hs37d5.fa
NOVOREF=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/human_b37_hs37d5.k14s1.novoindex
INDELS1=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
INDELS2=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz
DBSNP=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz
INTERVALS=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/output.intervals

# PBS parameters
NODE=$3
QUEUE=primary

# Software paths
NOVOALIGN=/shared/external_bin/novoalign
GATK=/shared/external_bin/GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar
SAMTOOLS=/shared/bin/samtools
PICARD=/mnt/thor_pool1/user_data/cc2qe/software/picard-tools-1.90
QUICK_Q=/mnt/thor_pool1/user_data/cc2qe/code/bin/quick_q

# ---------------------
# STEP 1: Allocate the data to the local drive

# copy the files to the local drive
# Require a lot of memory for this so we don't have tons of jobs writing to drives at once

# make the working directory
MOVE_FILES_CMD="mkdir -p $WORKDIR &&
rsync -rv $SAMPLEDIR/* $WORKDIR"

MOVE_FILES_CMD="echo hi"

#echo $MOVE_FILES_CMD
MOVE_FILES=`$QUICK_Q -m 512mb -d $NODE -t 1 -n moveTest -c " $MOVE_FILES_CMD " -q $QUEUE`

########### MAKE SUR EYO UFIX OIEHEWORIHJWO THE FASTQ FILEPATHS!!!!!!!!!!!
#zcat *_1.filt.fastq.gz | gzip -c > ${WORKDIR}/${SAMPLE}_1.fq.gz
#zcat *_2.filt.fastq.gz | gzip -c > ${WORKDIR}/${SAMPLE}_2.fq.gz


# ---------------------
# STEP 2: Align the fastq files with novoalign
# 12 cores and 16g of memory

ALIGN_CMD="cd $WORKDIR &&
for i in \$(seq 1 \`cat fqlist1 | wc -l\`)
do
    FASTQ1=\`sed -n \${i}p fqlist1\` &&
    FASTQ2=\`sed -n \${i}p fqlist2\` &&
    READGROUP=\`echo \$FASTQ1 | sed 's/_.*//g'\` &&

    RGSTRING=\`cat \${READGROUP}_readgroup.txt\` &&

    time $NOVOALIGN -d $NOVOREF -f \$FASTQ1 \$FASTQ2 \
	-r Random -c 12 -o sam \$RGSTRING | $SAMTOOLS view -Sb - > $SAMPLE.\$READGROUP.novo.bam ;
done"

ALIGN_CMD="echo hi"

echo $ALIGN_CMD
ALIGN=`$QUICK_Q -m 16gb -d $NODE -t 12 -n novo_$SAMPLE -c " $ALIGN_CMD " -q $QUEUE -z "-W depend=afterok:${MOVE_FILES}"`



# ---------------------
# STEP 3: Sort and fix flags on the bam file

# this only requires one core but a decent amount of memory.
SORT_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do

    time $SAMTOOLS view -bu $SAMPLE.\$READGROUP.novo.bam | \
	$SAMTOOLS sort -n -o - samtools_nsort_tmp | \
	$SAMTOOLS fixmate /dev/stdin /dev/stdout | $SAMTOOLS sort -o - samtools_csort_tmp | \
	$SAMTOOLS fillmd -b - $REF > $SAMPLE.\$READGROUP.novo.fixed.bam &&
    
    $SAMTOOLS index $SAMPLE.\$READGROUP.novo.fixed.bam &&
    rm $SAMPLE.\$READGROUP.novo.bam
done"

echo $SORT_CMD
SORT_CMD="echo hi"

SORT=`$QUICK_Q -m 8gb -d $NODE -t 1 -n sort_$SAMPLE -c " $SORT_CMD " -q $QUEUE -z "-W depend=afterok:$ALIGN"`



# ---------------------
# STEP 5: GATK reprocessing

    # make the set of regions for local realignment (don't need to do this step because it is unrelated tot eh alignment. Just need to do it once globally) . This would normally be after MarkDups.
    # java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp -jar $GATK -T RealignerTargetCreator -R $REF -o output.intervals -known $INDELS1 -known $INDELS2

GATK_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do
    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar INPUT=$SAMPLE.\$READGROUP.novo.fixed.bam OUTPUT=$SAMPLE.\$READGROUP.novo.fixed.mkdup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES=1000 CREATE_INDEX=true &&


    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK -T IndelRealigner -R $REF -I $SAMPLE.\$READGROUP.novo.fixed.mkdup.bam -o $SAMPLE.\$READGROUP.novo.realign.fixed.bam -targetIntervals $INTERVALS -known $INDELS1 -known $INDELS2 -LOD 0.4 -model KNOWNS_ONLY &&

    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
	-T CountCovariates \
	-R $REF \
	-I $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
	-recalFile $SAMPLE.\$READGROUP.recal_data.csv \
	-knownSites $DBSNP \
	-l INFO \
	-L '1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y;MT' \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate &&

    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
	-T TableRecalibration \
	-R $REF \
	-recalFile $SAMPLE.\$READGROUP.recal_data.csv \
	-I $SAMPLE.\$READGROUP.realign.fixed.bam \
	-o $SAMPLE.\$READGROUP.recal.bam \
	-l INFO \
	-noOQs \
	--disable_bam_indexing
done"


SORT=`$QUICK_Q -m 8gb -d $NODE -t 1 -n gatk_$SAMPLE -c " $GATK_CMD " -q $QUEUE -z "-W depend=afterok:$SORT"`

exit 0

# -----------------------
# STEP 6: quick calmd

cd $WORKDIR
for READGROUP in `cat rglist`
do
    # Calmd
    $SAMTOOLS calmd -Erb $SAMPLE.$READGROUP.recal.bam $REF > $SAMPLE.$READGROUP.recal.bq.bam
    $SAMTOOLS index $SAMPLE.$READGROUP.recal.bq.bam
done


# -----------------------
# STEP 7: Merging files
# Using Picard instead of samtools because it does a better job of preserving header information

cd $WORKDIR

INPUT_STRING=""
for READGROUP in `cat rglist`
do
    INPUT_STRING+=" I=$SAMPLE.$READGROUP.recal.bq.bam"
done

java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MergeSamFiles.jar $INPUT_STRING O=$SAMPLE.merged.bam SO=coordinate ASSUME_SORTED=true CREATE_INDEX=true


# -----------------------
# STEP 8: Mark duplicates again following the merge.


java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar INPUT=$SAMPLE.merged.bam OUTPUT=$SAMPLE.$READGROUP.novo.fixed.mkdup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES=1000 CREATE_INDEX=true    



# -----------------------
# STEP 9: Reduce reads


# ---------------------
# STEP 10: Move back to hall13 and cleanup.







