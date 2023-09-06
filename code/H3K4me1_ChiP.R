library(data.table)
library(seqinr)
library(ggplot2)
library(polymorphology)
source("code/libraries_functions.R")


genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genes$chr<-as.character(genes$chr)


setkey(genes, chr, start, stop)

# data from Niu et al. 2021, Nature Communications:
infiles<-c("data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz", 
           "data/GSM4668652_Col0_rep2_H3K4me1_Input_unique_reads.bedGraph.gz")
chipfiles<-c("data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
             "data/GSM4668651_Col0_rep2_H3K4me1_ChIP_unique_reads.bedGraph.gz")

inputs<-lapply(infiles, function(i){ 
  chip_overlaps(i, genes)
})

chips<-lapply(chipfiles, function(i){ 
  chip_overlaps(bedfile = i, featureobject = genes)
})

in1_total<-chip_total(infiles[1])
in2_total<-chip_total(infiles[2])
chip1_total<-chip_total(chipfiles[1])
chip2_total<-chip_total(chipfiles[2])

genes$chip1_H3K4me1<-rowSums(as.data.table(chips)[,c(1)], na.rm=F)
genes$input1_H3K4me1<-rowSums(as.data.table(inputs)[,c(1)], na.rm=F)
genes$enrich1_H3K4me1<-log2((1+genes$chip1_H3K4me1)/chip1_total) - log2((1+genes$input1_H3K4me1)/in1_total)

genes$chip2_H3K4me1<-rowSums(as.data.table(chips)[,c(2)], na.rm=F)
genes$input2_H3K4me1<-rowSums(as.data.table(inputs)[,c(2)], na.rm=F)
genes$enrich2_H3K4me1<-log2((1+genes$chip2_H3K4me1)/chip2_total) - log2((1+genes$input2_H3K4me1)/in2_total)
genes$enrich_H3K4me1<-rowMeans(genes[,c("enrich1_H3K4me1","enrich2_H3K4me1"), with=F], na.rm=T)

fwrite(genes, "data/A_thal_genes_PDS5_enrich.csv")

