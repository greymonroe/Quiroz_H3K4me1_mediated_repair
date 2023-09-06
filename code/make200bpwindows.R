library(data.table)
library(seqinr)
library(ggplot2)
library(polymorphology)
source("code/libraries_functions.R")


chr<-rbindlist(lapply(list.files("data/Arabidopsis Epigenome Profiles/", full.names = T, pattern="*200*"), function(f){
  chr<-fread(f)
  chr$CHROM<-as.character(chr$chr)
  chr$START<-chr$start
  chr$STOP<-chr$stop
  return(chr)
}))

chr$chr<-as.character(chr)
setkey(chr, chr, start, stop)

# data from Niu et al. 2021, Nature Communications:
infiles<-c("data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz", 
           "data/GSM4668652_Col0_rep2_H3K4me1_Input_unique_reads.bedGraph.gz")
chipfiles<-c("data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
             "data/GSM4668651_Col0_rep2_H3K4me1_ChIP_unique_reads.bedGraph.gz")

inputs<-lapply(infiles, function(i){ 
  chip_overlaps_window(i,  featureobject = chr)
})

chips<-lapply(chipfiles, function(i){ 
  chip_overlaps_window(bedfile = i, featureobject = chr)
})

in1_total<-chip_total(infiles[1])
in2_total<-chip_total(infiles[2])
chip1_total<-chip_total(chipfiles[1])
chip2_total<-chip_total(chipfiles[2])

chr$chip1_H3K4me1<-rowSums(as.data.table(chips)[,c(1)], na.rm=F)
chr$input1_H3K4me1<-rowSums(as.data.table(inputs)[,c(1)], na.rm=F)
chr$enrich1_H3K4me1<-log2((1+chr$chip1_H3K4me1)/chip1_total) - log2((1+chr$input1_H3K4me1)/in1_total)

chr$chip2_H3K4me1<-rowSums(as.data.table(chips)[,c(2)], na.rm=F)
chr$input2_H3K4me1<-rowSums(as.data.table(inputs)[,c(2)], na.rm=F)
chr$enrich2_H3K4me1<-log2((1+chr$chip2_H3K4me1)/chip2_total) - log2((1+chr$input2_H3K4me1)/in2_total)
chr$enrich_H3K4me1<-rowMeans(chr[,c("enrich1_H3K4me1","enrich2_H3K4me1"), with=F], na.rm=T)


# data from Niu et al. 2021, Nature Communications:
infiles<-c("data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph", 
           "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")
chipfiles<-c("data/GSM4668646_RDM15_FLAG_ChIP_rep1_unique_reads.bedGraph",
             "data/GSM4668648_RDM15_FLAG_ChIP_rep2_unique_reads.bedGraph", 
             "data/GSM4711900_RDM15_FLAG_rep1_uniq_more_reads.bedGraph",
             "data/GSM4711901_RDM15_FLAG_rep2_uniq_more_reads.bedGraph")

inputs<-lapply(infiles, function(i){ 
  chip_overlaps_window(i, chr)
})

chips1<-lapply(chipfiles[c(1,3)], function(i){ 
  chip_overlaps_window(bedfile = i, featureobject = chr)
})

chips2<-lapply(chipfiles[c(2,4)], function(i){ 
  chip_overlaps_window(i, chr)
})

in1_total<-chip_total(infiles[1])
in2_total<-chip_total(infiles[2])
chip1_total<-chip_total(chipfiles[1])+chip_total(chipfiles[3])
chip2_total<-chip_total(chipfiles[2])+chip_total(chipfiles[4])

chr$chip1_PDS5C<-rowSums(as.data.table(chips1), na.rm=F)
chr$input1_PDS5C<-rowSums(as.data.table(inputs)[,c(1)], na.rm=F)
chr$enrich1_PDS5C<-log2((1+chr$chip1_PDS5C)/chip1_total) - log2((1+chr$input1_PDS5C)/in1_total)

chr$chip2_PDS5C<-rowSums(as.data.table(chips2), na.rm=F)
chr$input2_PDS5C<-rowSums(as.data.table(inputs)[,c(2)], na.rm=F)
chr$enrich2_PDS5C<-log2((1+chr$chip2_PDS5C)/chip2_total) - log2((1+chr$input2_PDS5C)/in2_total)
chr$enrich_PDS5C<-rowMeans(chr[,c("enrich1_PDS5C","enrich2_PDS5C"), with=F], na.rm=T)


fwrite(chr, "data/200bpwindow_ChIP.csv")




