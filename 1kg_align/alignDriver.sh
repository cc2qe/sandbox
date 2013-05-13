#!/bin/bash

if [ $# -lt 3 ]
then
    echo usage $0 [sampleList] [batchDirectory] [node]
    exit 1
fi

# Directory and data names
SAMPLELIST=$1
BATCHDIR=$2
ROOTDIR=/scratch/cc2qe/1kg/batch1

# Annotations
REF=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/human_b37_hs37d5.fa
NOVOREF=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/human_b37_hs37d5.k14s1.novoindex
INDELS1=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
INDELS2=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz
DBSNP=/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz

# PBS parameters
NODE=$3
QUEUE=primary
MOVE_FILES_Q=$4

# Software paths
#NOVOALIGN=/shared/external_bin/novoalign
NOVOALIGN=/shared/external_bin/novocraft-2.08.03/novoalign
GATK=/shared/external_bin/GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar
SAMTOOLS=/shared/bin/samtools
PICARD=/mnt/thor_pool1/user_data/cc2qe/software/picard-tools-1.90
QUICK_Q=/mnt/thor_pool1/user_data/cc2qe/code/bin/quick_q

# ---------------------
# STEP 1: Allocate the data to the local drive

# copy the files to the local drive
# Require a lot of memory for this so we don't have tons of jobs writing to drives at once

MOVE_FILES_CMD="for SAMPLE in \`cat $SAMPLELIST\`
do 
    WORKDIR=$ROOTDIR/\$SAMPLE &&
    SAMPLEDIR=$BATCHDIR/\$SAMPLE &&

    mkdir -p \$WORKDIR &&
    rsync -rv \$SAMPLEDIR/fqlist* \$SAMPLEDIR/rglist \$SAMPLEDIR/*_readgroup.txt \$SAMPLEDIR/*.fastq.gz \$WORKDIR
done"

#MOVE_FILES_CMD="echo MOVE_FILES_CMD"

MOVE_FILES_Q=`$QUICK_Q -m 1gb -d $NODE -t 1 -n move_${NODE} -c " $MOVE_FILES_CMD " -q $QUEUE`


for SAMPLE in `cat $SAMPLELIST`
do
    WORKDIR=$ROOTDIR/$SAMPLE &&
    SAMPLEDIR=$BATCHDIR/$SAMPLE &&


# ---------------------
# STEP 2: Align the fastq files with novoalign
# 8 cores and 15500mb of memory

ALIGN_CMD="cd $WORKDIR &&
for i in \$(seq 1 \`cat fqlist1 | wc -l\`)
do
    FASTQ1=\`sed -n \${i}p fqlist1\` &&
    FASTQ2=\`sed -n \${i}p fqlist2\` &&
    READGROUP=\`echo \$FASTQ1 | sed 's/_.*//g'\` &&

    RGSTRING=\`cat \${READGROUP}_readgroup.txt\` &&

    time $NOVOALIGN -d $NOVOREF -f \$FASTQ1 \$FASTQ2 \
	-r Random -c 8 -o sam \$RGSTRING | $SAMTOOLS view -Sb - > $SAMPLE.\$READGROUP.novo.bam ;
done"

#ALIGN_CMD="echo ALIGN_CMD"

#echo $ALIGN_CMD

# 16gb mem
ALIGN_Q=`$QUICK_Q -d $NODE -t 8 -m 15500mb -n novo_${SAMPLE}_${NODE} -c " $ALIGN_CMD " -q $QUEUE -W depend=afterok:$MOVE_FILES_Q`


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
    
    time $SAMTOOLS index $SAMPLE.\$READGROUP.novo.fixed.bam &&
    rm $SAMPLE.\$READGROUP.novo.bam
done"

# echo $SORT_CMD
#SORT_CMD="echo SORT_CMD"

SORT_Q=`$QUICK_Q -m 8gb -d $NODE -t 1 -n sort_${SAMPLE}_${NODE} -c " $SORT_CMD " -q $QUEUE -z "-W depend=afterok:$ALIGN_Q"`


# ---------------------
# STEP 5: GATK reprocessing

GATK_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do
    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar \
         INPUT=$SAMPLE.\$READGROUP.novo.fixed.bam \
         OUTPUT=$SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
         ASSUME_SORTED=TRUE \
         METRICS_FILE=/dev/null \
         VALIDATION_STRINGENCY=SILENT \
         MAX_FILE_HANDLES=1000 \
         CREATE_INDEX=true &&

    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp -jar $GATK \
        -T RealignerTargetCreator \
        -nt 3 \
        -R $REF \
        -I $SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
        -o $SAMPLE.\$READGROUP.intervals \
        -known $INDELS1 \
        -known $INDELS2 &&

    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
         -T IndelRealigner \
         -R $REF \
         -I $SAMPLE.\$READGROUP.novo.fixed.mkdup.bam \
         -o $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
         -targetIntervals $SAMPLE.\$READGROUP.intervals \
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


    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
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
        $SAMPLE.\$READGROUP.novo.fixed.mkdup.bai \
        $SAMPLE.\$READGROUP.novo.realign.fixed.bam \
        $SAMPLE.\$READGROUP.novo.realign.fixed.bai

done"

#GATK_CMD="echo GATK_CMD"

GATK_Q=`$QUICK_Q -m 8gb -d $NODE -t 3 -n gatk_${SAMPLE}_${NODE} -c " $GATK_CMD " -q $QUEUE -W depend=afterok:$SORT_Q`



# -----------------------
# STEP 6: Samtools calmd

CALMD_CMD="cd $WORKDIR &&
for READGROUP in \`cat rglist\`
do
    time $SAMTOOLS calmd -Erb $SAMPLE.\$READGROUP.recal.bam $REF > $SAMPLE.\$READGROUP.recal.bq.bam &&
    time $SAMTOOLS index $SAMPLE.\$READGROUP.recal.bq.bam &&

    echo 'cleaning up...' &&
    rm $SAMPLE.\$READGROUP.recal.bam
done"

#CALMD_CMD="echo calmd_cmd"

CALMD_Q=`$QUICK_Q -m 1gb -d $NODE -t 1 -n calmd_${SAMPLE}_${NODE} -c " $CALMD_CMD " -q $QUEUE -W depend=afterok:$GATK_Q`



# -----------------------
# STEP 7: Merging files
# Using Picard instead of samtools because it does a better job of preserving header information

# only merge if there are multiple readgroups
MERGE_CMD="
cd $WORKDIR &&
if [ \`cat rglist | wc -l\` -gt 1 ] ;
then
    INPUT_STRING='' &&
    for READGROUP in \`cat rglist\`
    do
        INPUT_STRING+=\" I=$SAMPLE.\$READGROUP.recal.bq.bam\"
    done &&

    time java -Xmx4g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MergeSamFiles.jar \$INPUT_STRING O=$SAMPLE.merged.bam SO=coordinate ASSUME_SORTED=true CREATE_INDEX=true &&

    for READGROUP in \`cat rglist\`
    do
        rm $SAMPLE.\$READGROUP.recal.bq.bam $SAMPLE.\$READGROUP.recal.bq.bam.bai
    done
else
    for READGROUP in \`head -n 1 rglist\`
    do
        mv $SAMPLE.\$READGROUP.recal.bq.bam $SAMPLE.novo.bam &&
        mv $SAMPLE.\$READGROUP.recal.bq.bam.bai $SAMPLE.novo.bai
    done
fi
"

#MERGE_CMD="echo merge_cmd command"

MERGE_Q=`$QUICK_Q -m 4gb -d $NODE -t 1 -n merge_${SAMPLE}_${NODE} -c " $MERGE_CMD " -q $QUEUE -W depend=afterok:$CALMD_Q`


# -----------------------
# STEP 8: Mark duplicates again following the merge.

MKDUP2_CMD="
cd $WORKDIR &&
if [ \`cat rglist | wc -l\` -gt 1 ] ;
then
    time java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $PICARD/MarkDuplicates.jar INPUT=$SAMPLE.merged.bam OUTPUT=$SAMPLE.novo.bam ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES=1000 CREATE_INDEX=true &&

    echo 'clean up files...' &&
    rm $SAMPLE.merged.bam $SAMPLE.merged.bai
fi
"

#MKDUP2_CMD="echo mkdup2 command"

MKDUP2_Q=`$QUICK_Q -m 8gb -d $NODE -t 1 -n mkdup2_${SAMPLE}_${NODE} -c " $MKDUP2_CMD " -q $QUEUE -W depend=afterok:$MERGE_Q`


# -----------------------
# STEP 9: Reduce reads

REDUCE_CMD="cd $WORKDIR &&
time java -Xmx16g -Djava.io.tmpdir=$WORK_DIR/tmp/ -jar $GATK \
    -T ReduceReads \
    -R $REF \
    -I $SAMPLE.novo.bam \
    -o $SAMPLE.novo.reduced.bam"

#REDUCE_CMD="echo reduce command"

REDUCE_Q=`$QUICK_Q -m 16gb -d $NODE -t 1 -n reduce_${SAMPLE}_${NODE} -c " $REDUCE_CMD " -q $QUEUE -W depend=afterok:$MKDUP2_Q`


# ---------------------
# STEP 10: Move back to hall13 and cleanup.

RESTORE_CMD="cd $WORKDIR &&
rsync -rv $SAMPLE.novo.bam $SAMPLE.novo.bai \
    $SAMPLE.novo.reduced.bam $SAMPLE.novo.reduced.bai \
    *.recal_data.grp \
    *.intervals \
    $SAMPLEDIR &&

echo 'removing scratch directory...' &&
rm -r $WORKDIR &&

echo $SAMPLE >> $SAMPLEDIR/../completed.txt"

#RESTORE_CMD="echo RESTORE_CMD"

RESTORE_Q=`$QUICK_Q -m 1gb -d $NODE -t 1 -n restore_${SAMPLE}_${NODE} -c " $RESTORE_CMD " -q $QUEUE -W depend=afterok:$REDUCE_Q`


done






