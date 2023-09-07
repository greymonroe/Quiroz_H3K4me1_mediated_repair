library(polymorphology)

files<-list.files("data/strelka_snpeff/download", full.names = T)

PASS3<-fread("data/MSH6_KO_mutations.csv");PASS3$CHROM<-as.character(PASS3$CHROM)

snpeff<-rbindlist(lapply(files, function(file){
  dat<-read.vcfR(file)
  calls<-cbind(data.table(dat@fix),data.table(dat@gt))
  calls$unique<-paste(calls$CHROM, calls$POS, sep = "_")
  calls$passed<-calls$unique %in% PASS3$unique
  return(calls[passed==T])
}))

snpeff$unique2<-paste(snpeff$CHROM, snpeff$POS, snpeff$ALT, sep = "_")
snpeff$synonymous<-grepl("synonymous", snpeff$INFO)
snpeff$nonsynonymous<-grepl("missense", snpeff$INFO)
snpeff$effect<-"non-coding region"
snpeff$effect[which(snpeff$synonymous)]<-"synonymous"
snpeff$effect[which(snpeff$nonsynonymous)]<-"nonsynonymous"

snpeff_cleaned<-snpeff[FILTER=="PASS", .(effect=unique(effect), nonsynonymous=unique(nonsynonymous), synonymous=unique(synonymous)), by=unique2]

sum(snpeff_cleaned$nonsynonymous)
PASS3$unique2<-paste(PASS3$CHROM, PASS3$POS, PASS3$ALT, sep = "_")
PASS3_effect<-merge(PASS3, snpeff_cleaned)
fwrite(PASS3_effect, "data/Strelka2_mutations_SNPeff.csv")


MSH6_genotypes<-fread("data/MSH6_genotypes.csv")

PASS3_effect$line<-MSH6_genotypes$line[match(PASS3_effect$file, MSH6_genotypes$file)]
PASS3_effect$genotype<-MSH6_genotypes$genotype[match(PASS3_effect$file, MSH6_genotypes$file)]
PASS3_effect$trt<-MSH6_genotypes$trt[match(PASS3_effect$file, MSH6_genotypes$file)]

chr<-fread("data/200bpwindow_ChIP.csv")
chr$chr<-as.character(chr$chr)
setkey(chr, chr, start, stop)
setkey(PASS3_effect, CHROM, start, stop)
PASS3_effect$ID
PASS3_effect_overlap<-foverlaps(chr, PASS3_effect)
PASS3_effect_overlap[!is.na(ID),
                  .(enrich_H3K4me1=mean(enrich_H3K4me1), 
                    enrich1_H3K4me1=mean(enrich1_H3K4me1),
                    enrich2_H3K4me1=mean(enrich2_H3K4me1)), 
                  by=ID]
PASS3_effect$enrich_H3K4me1=PASS3_effect_overlap[!is.na(ID),.(enrich_H3K4me1=mean(enrich_H3K4me1)), by=ID]$enrich_H3K4me1
PASS3_effect$enrich1_H3K4me1=PASS3_effect_overlap[!is.na(ID),.(enrich1_H3K4me1=mean(enrich1_H3K4me1)), by=ID]$enrich1_H3K4me1
PASS3_effect$enrich2_H3K4me1=PASS3_effect_overlap[!is.na(ID),.(enrich2_H3K4me1=mean(enrich2_H3K4me1)), by=ID]$enrich2_H3K4me1
PASS3_effect$H3K4me1_enriched<-PASS3_effect$enrich_H3K4me1>0

fwrite(PASS3_effect,"data/Strelka2_mutations.csv")


