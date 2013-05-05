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
MOVE_FILES_Q=`$QUICK_Q -m 512mb -d $NODE -t 1 -n moveTest -c " $MOVE_FILES_CMD " -q $QUEUE`

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
ALIGN_Q=`$QUICK_Q -m 16gb -d $NODE -t 12 -n novo_$SAMPLE -c " $ALIGN_CMD " -q $QUEUE -z "-W depend=afterok:$MOVE_FILES_Q"`



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

SORT_Q=`$QUICK_Q -m 8gb -d $NODE -t 1 -n sort_$SAMPLE -c " $SORT_CMD " -q $QUEUE -z "-W depend=afterok:$ALIGN_Q"`



# ---------------------
# STEP 5: GATK reprocessing

GATK_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do
    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar \
         INPUT=$SAMPLE.\$READGROUP.novo.fixed.bam \
         OUTPUT=$SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
         ASSUME_SORTED=TRUE \
         METRICS_FILE=/dev/null \
         VALIDATION_STRINGENCY=SILENT \
         MAX_FILE_HANDLES=1000 \
         CREATE_INDEX=true &&

    echo 'make the set of regions for local realignment (don't need to do this step because it is unrelated tot eh alignment. Just need to do it once globally).' &&
    echo 'java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp -jar $GATK -T RealignerTargetCreator -R $REF -o output.intervals -known $INDELS1 -known $INDELS2' &&

    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
         -T IndelRealigner \
         -R $REF \
         -I $SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
         -o $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
         -targetIntervals $INTERVALS \
         -known $INDELS1 \
         -known $INDELS2 \
         -LOD 0.4 \
         -model KNOWNS_ONLY &&

    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
        -T BaseRecalibrator \
        -nct 3 \
        -I $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
        -R $REF \
        -knownSites $DBSNP \
        -l INFO \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        -cov ContextCovariate \
        -o $SAMPLE.\$READGROUP.recal_data.grp &&


    java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
        -T PrintReads \
        -R $REF \
        -I $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
        -BQSR $SAMPLE.\$READGROUP.recal_data.grp \
        --disable_bam_indexing \
        -l INFO \
        -o $SAMPLE.\$READGROUP.recal.bam &&

    echo 'cleaning up...' &&
    rm $SAMPLE.\$READGROUP.novo.fixed.bam \
        $SAMPLE.\$READGROUP.novo.fixed.bam.bai \
        $SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
        $SAMPLE.\$READGROUP.novo.fixed.mkdup.bai

done"

GATK_CMD="echo hi"

GATK_Q=`$QUICK_Q -m 8gb -d $NODE -t 3 -n gatk_$SAMPLE -c " $GATK_CMD " -q $QUEUE -z "-W depend=afterok:$SORT_Q"`



# -----------------------
# STEP 6: quick calmd

CALMD_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do
    $SAMTOOLS calmd -Erb $SAMPLE.\$READGROUP.recal.bam $REF > $SAMPLE.\$READGROUP.recal.bq.bam &&
    $SAMTOOLS index $SAMPLE.\$READGROUP.recal.bq.bam
done"

CALMD_Q=`$QUICK_Q -m 512mb -d $NODE -t 1 -n calmd_$SAMPLE -c " $CALMD_CMD " -q $QUEUE -W depend-afterok:$GATK_Q`


# -----------------------
# STEP 7: Merging files
# Using Picard instead of samtools because it does a better job of preserving header information

MERGE_CMD="cd $WORKDIR &&

INPUT_STRING='' &&
for READGROUP in \`cat rglist\`
do
    INPUT_STRING+=' I=$SAMPLE.\$READGROUP.recal.bq.bam'
done &&

java -Xmx4g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MergeSamFiles.jar \$INPUT_STRING O=$SAMPLE.merged.bam SO=coordinate ASSUME_SORTED=true CREATE_INDEX=true"

MERGE_Q=`$QUICK_Q -m 4gb -d $NODE -t 1 -n merge_$SAMPLE -c " $MERGE_CMD " -q $QUEUE -W depend-afterok:$CALMD_Q`


# -----------------------
# STEP 8: Mark duplicates again following the merge.

MKDUP2_CMD="cd $WORKDIR &&

time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar INPUT=$SAMPLE.merged.bam OUTPUT=$SAMPLE.merged.mkdup.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES=1000 CREATE_INDEX=true"

MKDUP2_Q=`$QUICK_Q -m 8gb -d $NODE -t 1 -n mkdup2_$SAMPLE -c " $MKDUP2_cmd " -q $QUEUE -W depend-afterok:$MERGE_Q`


# -----------------------
# STEP 9: Reduce reads


# ---------------------
# STEP 10: Move back to hall13 and cleanup.







