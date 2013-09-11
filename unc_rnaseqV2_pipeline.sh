#!/bin/bash -e

if [ $# -lt 4 ]
then
    echo usage $0 [SAMPLE_NAME] [FASTQ1] [FASTQ2] [THREADS]
    exit 1
fi

SAMPLE=$1
FASTQ1=$2
FASTQ2=$3
THREADS=$4

echo "Processing $SAMPLE from seq files $FASTQ1, $FASTQ2"

# 0. unzip the fastqs
zcat $FASTQ1 | sed 's/^@.* /@/g' > ${SAMPLE}_1.fastq
zcat $FASTQ2 | sed 's/^@.* /@/g' > ${SAMPLE}_2.fastq

# 1. Format fastq 1 for Mapsplice
## omitted because my files are alread phred33 and i will process them that way
# java -Xmx512M -jar ubu.jar fastq-format --phred33to64 --strip --suffix /1 –in raw_1.fastq --out working/prep_1.fastq > working/mapsplice_prep1.log

# 2. Format fastq 2 for Mapsplice
## omitted because my files are already phred33 and i will process them that way.
# java -Xmx512M -jar ubu.jar fastq-format --phred33to64 --strip --suffix /2 –in raw_2.fastq --out working/prep_2.fastq > working/mapsplice_prep2.log

# 3. Mapsplice
echo "3. Mapsplice"
## example command. this is for Mapsplice 12_07 though, i'm modifying it for Mapsplice 2.0.1.9, which they used in more recent samples.
# python mapsplice_multi_thread.py --fusion --all-chromosomes-files hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-filesdir hg19_M_rCRS/chromosomes --Bowtieidx hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o SAMPLE_BARCODE 2> working/mapsplice.log
## actual command
time python ~/software/MapSplice_multi_threads_2.0.1.9/mapsplice.py --fusion --bam -p $THREADS -c ~/refdata/genomes/unc_tcga_hg19/chromosomes --qual-scale phred33 -x ~/refdata/genomes/unc_tcga_hg19/ebwt/humanchridx_M_rCRS -1 ${SAMPLE}_1.fastq -2 ${SAMPLE}_2.fastq -o $SAMPLE.mapsplice > mapsplice.log 2> mapsplice.log

# 4. Add read groups
echo "4. Add read groups"
## omitting these tags because i don't know them for my samples:
time java -Xmx4G -jar ~/software/picard-tools-1.96/AddOrReplaceReadGroups.jar INPUT=$SAMPLE.mapsplice/alignments.bam OUTPUT=${SAMPLE}.rg.alignments.bam RGSM=${SAMPLE} RGID=${SAMPLE} RGLB=TruSeq RGPL=illumina RGPU=barcode VALIDATION_STRINGENCY=SILENT TMP_DIR=add_rg_tag_tmp > add_rg_tag.log 2> add_rg_tag.log

## omitting this step because it's stupid. and i aligned in phred33 to begin with
# 5. Convert back to phred33
# java -Xmx512M -jar ubu.jar sam-convert --phred64to33 --in working/rg_alignments.bam –out working/phred33_alignments.bam > working/sam_convert.log 2> working/sam_convert.log

# 6. Sort by coordinate
echo "6. Sort by coordinate"
## I'm gonna use picard instead of samtools
# samtools sort ${SAMPLE}_rg_alignments.bam ${SAMPLE}.genome.aln
java -Xmx8g -jar -Djava.io.tmpdir=tmp ~/software/picard-tools-1.96/SortSam.jar I=${SAMPLE}.rg.alignments.bam O=${SAMPLE}.genome.aln.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

# 7. Flagstat
echo "7. Flagstat"
samtools flagstat ${SAMPLE}.genome.aln.bam > ${SAMPLE}.genome.aln.bam.flagstat

# 8. Index
## don't need to do this because I've already indexed it with picard.
# samtools index ${SAMPLE}.genome.aln.bam

# 9. Sort by chromosome, then read id
echo "Sort by chromosome, then read id"
time perl ~/software/ubu/src/perl/sort_bam_by_reference_and_name.pl --input ${SAMPLE}.genome.aln.bam --output ${SAMPLE}.alignments.chromReadSorted.bam --temp-dir . --samtools /shared/bin/samtools > sorted_by_chr_read.log 2> sorted_by_chr_read.log

# 10. Translate to transcriptome coords
echo "10. Translate to transcriptome coords"
time java -Xms8g -Xmx8g -jar ~/software/ubu-1.2/ubu-1.2-jar-with-dependencies.jar sam-xlate --bed ~/refdata/genomes/unc_tcga_hg19/unc_hg19.bed --in ${SAMPLE}.alignments.chromReadSorted.bam --out ${SAMPLE}.transcriptome.aln.bam --order ~/refdata/genomes/unc_tcga_hg19/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse > genome_to_transcriptome.log 2> genome_to_transcriptome.log

# 11. Filter indels, large inserts, zero mapping quality from transcriptome bam
echo "11. Filter indels, large inserts, zero mapping quality from transcriptome bam"
time java -Xmx1g -jar ~/software/ubu-1.2/ubu-1.2-jar-with-dependencies.jar sam-filter --in ${SAMPLE}.transcriptome.aln.bam --out ${SAMPLE}.transcriptome.aln.filtered.bam --strip- indels --max-insert 10000 --mapq 1 > sam_filter.log 2> sam_filter.log

# 12. RSEM
echo "12. RSEM"
## i don't know why they have a --gcr-output-file flag in their command, but it is not a valid rsem parameter so I'm omitting it
time ~/software/rsem-1.1.13/rsem-calculate-expression --paired-end --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.aln.filtered.bam ~/refdata/genomes/unc_tcga_hg19/rsem_ref/hg19_M_rCRS_ref $SAMPLE.rsem > rsem.log 2> rsem.log
## example cmd
# rsem-calculate-expression --gcr-output-file --paired-end --bam --estimate-rspd -p 8 working/transcriptome_alignments_filtered.bam /datastore/tier1data/nextgenseq/seqware-analysis/mapsplice_rsem/rsem_ref/hg19_M_rCRS_ref rsem > working/rsem.log 2> working/rsem.log

# 13. Strip trailing tabs from rsem.isoforms.results
perl ~/software/bin/strip_trailing_tabs.pl --input $SAMPLE.rsem.isoforms.results --temp $SAMPLE.rsem.orig.isoforms.results > trim_isoform_tabs.log 2> trim_isoform_tabs.log

# 14. Prune isoforms from gene quant file
mv $SAMPLE.rsem.genes.results $SAMPLE.rsem.orig.genes.results; sed /^uc0/d $SAMPLE.rsem.orig.genes.results > $SAMPLE.rsem.genes.results

# 15. Normalize gene quant
perl ~/software/bin/quartile_norm.pl -c 2 -q 75 -t 1000 -o $SAMPLE.rsem.genes.normalized_results $SAMPLE.rsem.genes.results

# 16. Normalize isoform quant
perl ~/software/bin/quartile_norm.pl -c 2 -q 75 -t 300 -o $SAMPLE.rsem.isoforms.normalized_results $SAMPLE.rsem.isoforms.results

# 17. Junction counts
java -Xmx512M -jar ~/software/ubu-1.2/ubu-1.2-jar-with-dependencies.jar sam-junc --junctions ~/refdata/genomes/unc_tcga_hg19/splice_junctions.txt --in $SAMPLE.genome.aln.bam --out junction_quantification.txt > junction_quantification.log 2> junction_quantification.log

# 18. Exon counts
coverageBed -split -abam $SAMPLE.genome.aln.bam -b ~/refdata/genomes/unc_tcga_hg19/composite_exons.bed | perl ~/software/bin/normalizeBedToolsExonQuant.pl $SAMPLE.genome.aln.bam ~/refdata/genomes/unc_tcga_hg19/composite_exons.bed > bt.exon_quantification.txt 2> bt_exon_quantification.log

# 19. Cleanup large intermediate output
rm $SAMPLE.mapsplice/alignments.bam $SAMPLE.rg.alignments.bam $SAMPLE.alignments.chromReadSorted.bam $SAMPLE.transcriptome.aln.bam $SAMPLE.transcriptome.aln.filtered.bam ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq