

conda create -n chipseq_bw_analysis
conda activate chipseq_bw_analysis
conda install -c bioconda ucsc-bigwigtobedgraph bedtools macs2 samtools


samtools faidx reference_genome.fa
cut -f1,2 reference_genome.fa.fai > reference_genome.chrom.sizes

chipseq_data=$1
input_control=$2

conda activate chipseq_bw_analysis


input_control=GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bw
input=$(basename "$input_control" .bw)
bigWigToBedGraph ${input_control} ${input}.bedGraph
bedtools bedtobam -i ${input}.bedGraph -g ../reference_genome.chrom.sizes > ${input}.bam
samtools sort ${input}.bam -o ${input}_sorted.bam
samtools index ${input}_sorted.bam




macs2 callpeak -t ${chip}_sorted.bam -c ${input}_sorted.bam -f BAM -g 374471240 -n ${chip}_chipseq_peaks



macs2 callpeak -t chipseq_data_sorted.bam -c input_control_sorted.bam -f BAM -g 374471240 --broad --broad-cutoff 0.1 -n chipseq_broad_peaks --outdir macs2_output_broad

samtools view ${chip}_sorted.bam

files=(GSM3674680_ChIP_Rice_7days_leaf_H3.bw  GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bw GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bw GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bw GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bw GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bw GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bw)

for file in $files; do 

echo $file

chipseq_data=GSE128434_RAW/$file

chip=$(basename "$chipseq_data" .bw)

bigWigToBedGraph ${chipseq_data} ${chip}.bedGraph

bedtools bedtobam -i ${chip}.bedGraph -g ../reference_genome.chrom.sizes >  ${chip}.bam

samtools sort ${chip}.bam -o ${chip}_sorted.bam

samtools index ${chip}_sorted.bam

macs2 callpeak -t ${chip}_sorted.bam -c ${input}_sorted.bam -f BAM -g 374471240 -n ${chip}_chipseq_peaks --nomodel --extsize 147 --outdir narrow

done


GSM3674680_ChIP_Rice_7days_leaf_H3.bw  GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bw GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bw GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bw GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bw GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bw GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bw
