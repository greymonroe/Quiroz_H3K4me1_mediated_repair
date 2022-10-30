#!/bin/bash

PREFIX=$1

module load trimmomatic/0.39

java -jar $TRIMMOMATIC_HOME/trimmomatic-0.39.jar PE -threads 8 -phred33 \
${PREFIX}_L002_R1_001.fastq.gz ${PREFIX}_L002_R3_001.fastq.gz \
1_fastq/${PREFIX}_1.trimmed.fastq.gz 1_fastq/${PREFIX}_1un.trimmed.fastq.gz \
1_fastq/${PREFIX}_2.trimmed.fastq.gz 1_fastq/${PREFIX}_2un.trimmed.fastq.gz \
SLIDINGWINDOW:4:20

# Use the trimmed PAIRED reads for next steps
sbatch bwa_PREFIX.sh $PREFIX

REF=TAIR10_chr_all.fasta
READ1=1_fastq/${PREFIX}_1.trimmed.fastq.gz
READ2=1_fastq/${PREFIX}_2.trimmed.fastq.gz

bwa mem -t 32 -r "@RG\tID:$PREFIX\tSM:$PREFIX\tPL:ILLUMINA" $REF $READ1 $READ2 | samtools sort -n -@5 -o bam/$PREFIX.bam

samtools fixmate -m bam/$PREFIX.bam - | samtools sort -@5 -o bam/$PREFIX.fix.bam

samtools index bam/$PREFIX.fix.bam

samtools markdup -@5 -s bam/$PREFIX.fix.bam - | samtools sort -@5 -o  bam/$PREFIX.fix.markdup.bam

samtools index bam/$PREFIX.fix.markdup.bam


