
setwd("~/Dropbox/Research/rice mutation paper/")

genome<-read.fasta("data/Osativa_204_softmasked.fa.gz")

#load rice data

LoF<-data.table(read.xlsx("data/Genes-Mutated-by-Loss-of-Function-Mutations-1504lines.xlsx", sheet=1, startRow = 4))
gene_annotations_all<-fread("data/all.locus_brief_info.7.0")
gene_annotations_all$length<-gene_annotations_all$stop-gene_annotations_all$start
gene_annotations_all$window_ID<-1:nrow(gene_annotations_all)
gene_annotations_all$ID<-1:nrow(gene_annotations_all)
setkey(gene_annotations_all, chr, start, stop)
TE<-gene_annotations_all[is_TE=="Y" & is_representative=="Y"]
gene_annotations_basic<-gene_annotations_all[is_TE=="N" & is_representative=="Y"]
gene_annotations_basic$direction<-gene_annotations_basic$ori
gene_annotations_basic$type="gene"

snps<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=1, startRow = 4))
snps$CHROM<-snps$Chromosome
snps$POS<-snps$Position
snps$unique<-paste0(snps$CHROM, snps$POS, snps$Single.base.substitution)
snps$chr<-snps$CHROM
snps$start<-snps$POS
snps$stop<-snps$POS
snps<-snps[!is.na(start)]
snps$REF<-substr(snps$Single.base.substitution,1,1)
snps$ALT<-substr(snps$Single.base.substitution,3,3)
snps$context<-contexts(snps, genome)
context_table<-data.table(table(context=snps$context))
context_table$context_only<-substr(context_table$context, 1, 3)
context_table$mut<-paste(substr(context_table$context, 2,2),substr(context_table$context, 6,6), sep=">")
Athal<-fread("data/A_thal_germ_SBS.txt")
setkey(snps, chr, start, stop)

del<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=2, startRow = 2))
del$CHROM<-del$Chromosome
del$POS<-del$Start.position
ins<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=3, startRow = 2))
ins$CHROM<-ins$Chromosome
ins$POS<-ins$Position
indels<-rbind(ins[`Size.(bp)`<10], del[`Size.(bp)`<10], fill=T)
indels$unique<-paste0(indels$CHROM, indels$POS)
indels$start<-indels$POS
indels$stop<-indels$POS
indels$chr<-indels$CHROM
indels$unique<-paste0(indels$CHROM, indels$POS)
indels<-indels[!is.na(start)]
indels$ID<-1:nrow(indels)
setkey(indels, chr, start, stop)

#http://glab.hzau.edu.cn/RiceENCODE/download/Peaks/Nipponbare/Histone/

H3K4me1<-read_encode("data/H3K4me1_seedlings.bed")
H3K4me1_panicles<-read_encode("data/H3K4me1_panicles.bed")
H3K4me1_leaves<-read_encode("data/H3K4me1_mature_leaves.bed")
H3K9me1<-read_encode("data/SRR3213598_H3K9me1_seedlings.bed")
H3K36me3<-read_encode("data/SRR094791_H3K36me3_seedlings.bed")
H3K4me3<-read_encode("data/H3K4me3_seedlings.bed")
H3K9me2<-read_encode("data/H3K9me2_seedlings.bed")
H3K27ac<-read_encode("data/H3K27ac_seedlings.bed")
PII<-read_encode("data/PII_seedlings.bed")
H3K27me3<-read_encode("data/H3K27me3_seedlings.bed")
H3K4ac<-read_encode("data/SRR3213599_H3K4ac_seedlings.bed")
H3K12ac<-read_encode("data/SRR6510886_H4K12ac_seedlings.bed")
H3K9ac<-read_encode("data/SRR6795643_H3K9ac_seedlings.bed")

missense<-fread("data/rice_missense.txt");colnames(missense)<-c("CHROM","POS")
missense$chr<-paste0('Chr',missense$CHROM)
missense$start<-missense$POS
missense$stop<-missense$POS
missense$Mutant.ID<-1:nrow(missense)
setkey(missense, chr, start, stop)

synonymous<-fread("data/rice_synonymous.txt");colnames(synonymous)<-c("CHROM","POS")
synonymous$chr<-paste0('Chr',synonymous$CHROM)
synonymous$start<-synonymous$POS
synonymous$stop<-synonymous$POS
synonymous$Mutant.ID<-1:nrow(synonymous)
setkey(synonymous, chr, start, stop)


ns_s<-fread("data/Mutations_Final_StopCodonAnnotated_IDAdded.csv")
ns_s$CHROM<-ns_s$Chromosome
ns_s$POS<-ns_s$Position
ns_s$unique<-paste0(ns_s$CHROM, ns_s$POS, ns_s$Single.base.substitution)
ns_s$chr<-ns_s$CHROM
ns_s$start<-ns_s$POS
ns_s$Mutant.ID<-ns_s$`Mutant ID`
ns_s$stop<-ns_s$POS
ns_s$gene<-gsub("\\..+", "", ns_s$ID)
ns<-ns_s[MutationType=="non-synonymous"]
s<-ns_s[MutationType=="synonymous"]

setkey(ns_s, chr, start, stop)

