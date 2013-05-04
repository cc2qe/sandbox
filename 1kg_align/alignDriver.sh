#!/bin/bash

# Directory and data names
ROOTDIR='/scratch/cc2qe/1kg/batch1'
WORKDIR=${ROOTDIR}/$SAMPLE
SAMPLE=$1
SAMPLEDIR=$2

# Annotations
REF='/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/human_b37_hs37d5.k14s1.novoindex'
INDELS1='/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz'
INDELS2='/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz'
DBSNP='/mnt/thor_pool1/user_data/cc2qe/refdata/genomes/b37/annotations/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz'

# PBS parameters
NODE=$2
QUEUE=full

# Software paths
NOVOALIGN=/shared/external_bin/novoalign
GATK=/shared/external_bin/GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar
SAMTOOLS=/shared/bin/samtools

# ---------------------
# STEP 1: Allocate the data to the local drive

# copy the files to the local drive
# Require a lot of memory for this so we don't have tons of jobs writing to drives at once

# make the working directory
mkdir $WORKDIR
rsync -rv $SAMPLEDIR/* $WORKDIR


########### MAKE SUR EYO UFIX OIEHEWORIHJWO THE FASTQ FILEPATHS!!!!!!!!!!!
#zcat *_1.filt.fastq.gz | gzip -c > ${WORKDIR}/${SAMPLE}_1.fq.gz
#zcat *_2.filt.fastq.gz | gzip -c > ${WORKDIR}/${SAMPLE}_2.fq.gz

# change directory to the sample working directory
cd $WORKDIR


# ---------------------
# STEP 2: Align the fastq files with novoalign
# 12 cores and 16g of memory

for i in $(seq 1 `cat fqlist1 | wc -l`)
do
    FASTQ1=`sed -n ${i}p fqlist1`
    FASTQ2=`sed -n ${i}p fqlist2`
    READGROUP=`echo $FASTQ1 | sed 's/_.*//g'`

    # readgroup string
    RGSTRING=`cat ${READGROUP}_readgroup.txt`

    $NOVOALIGN -d $REF -f ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz \
	-r Random -c 12 -o sam $RGSTRING | $SAMTOOLS view -Sb - > $SAMPLE.novo.bam
done

# ---------------------
# STEP 3: Sort and fix flags on the bam file

# this only requires one core but a decent amount of memory.
$SAMTOOLS view -bu $SAMPLE.novo.bam | \
    $SAMTOOLS sort -n -o - samtools_nsort_tmp | \
    $SAMTOOLS fixmate /dev/stdin /dev/stdout | $SAMTOOLS sort -o - samtools_csort_tmp | \
    $SAMTOOLS fillmd -u - $REF > $SAMPLE.novo.fixed.bam

# index the bam file
$SAMTOOLS index $SAMPLE.novo.fixed.bam



# ---------------------
# STEP 5: GATK reprocessing

# make the set of regions for local realignment
java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp -jar $GATK -T RealignerTargetCreator -R $REF -o output.intervals -known $INDELS1 -known $INDELS2

# perform the local realignment
java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK -T IndelRealigner -R $REF -I $SAMPLE.novo.fixed.bam -o $SAMPLE.novo.realign.fixed.bam -targetIntervals output.intervals -known $INDELS1 -known $INDELS2 -LOD 0.4 -model KNOWNS_ONLY -compress 0

# Count covariates
java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
    -T CountCovariates \
    -R $REF \
    -I $SAMPLE.novo.realign.fixed.bam \
    -recalFile $SAMPLE.recal_data.csv \
    -knownSites $DBSNP \
    -l INFO \
    -L '1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;X;Y;MT' \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate

# Base recalibration
java -Xmx8g -Djava.io.tmpdir=$WORKDIR/tmp/ -jar $GATK \
    -T TableRecalibration \
    -R $REF \
    -recalFile $SAMPLE.recal_data.csv \
    -I $SAMPLE.realign.fixed.bam \
    -o $SAMPLE.recal.bam \
    -l INFO \
    -noOQs \
    --disable_bam_indexing \
    -compress 0


# -----------------------
# STEP 6: quick calmd

# Calmd
$SAMTOOLS calmd -Erb $SAMPLE.recal.bam $REF > $SAMPLE.recal.bq.bam


# ----------------------
# STEP 7: Merging files






# ---------------------
# STEP 8: Move back to hall13 and cleanup.







