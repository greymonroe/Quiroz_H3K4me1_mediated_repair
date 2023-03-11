library(polymorphology)
library(vcfR)
library(seqinr)


genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genome<-read.fasta("data/TAIR10_chr_all.fas.gz")
names(genome)<-gsub("Chr","",names(genome))
bp<- c("A","C","G","T")
##

PASS3<-fread("data/MSH6_KO_mutations.csv");PASS3$CHROM<-as.character(PASS3$CHROM)
PASS3$intron<-features_overlap_mutation(intron, PASS3)
PASS3$exon<-features_overlap_mutation(exon, PASS3)

somatic<-PASS3[type=='somatic']
somatic$CHROM<-as.character(somatic$CHROM)
setkey(somatic, CHROM, start, stop)

somatic[,.(N=.N), by=.(file, trt)][,.(N=mean(N), se=sd(N)/sqrt(.N)), by=.(trt)]

germline<-PASS3[type!='somatic']
germline$CHROM<-as.character(germline$CHROM)
setkey(germline, CHROM, start, stop)

Athal_mut<-fread("data/A_thaliana_germline_mutations.txt")[!is.na(CHROM)]
Athal_mut$CHROM<-as.character(Athal_mut$CHROM)
Athal_mut$context_SBS<-contexts(vars = Athal_mut ,fasta = genome)


pdf("figures/MSH6_KO_SBS_mutonly.pdf", width=1.2, height=2)
contexts_plot(somatic[trt=="WT" ]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations", limits=c(0,170))+ggtitle("wt/wt + wt/msh6")
contexts_plot(somatic[trt=="MSH6"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations", limits=c(0,170))+ggtitle("msh6/msh6")

contexts_plot(Athal_mut[!is.na(context_SBS)]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("germline old\nall")
contexts_plot(germline[trt=="WT" ]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("germline WT\nall")
contexts_plot(germline[trt=="MSH6" ]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("germline MSH6\nall")

contexts_plot(somatic[trt=="WT" & genic==T]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("Wildtype\ngenic")
contexts_plot(somatic[trt=="MSH6"& genic==T]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("MSH6 KO\ngenic")

contexts_plot(somatic[trt=="WT" & genic==F]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("Wildtype\nnon-genic")
contexts_plot(somatic[trt=="MSH6"& genic==F]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("MSH6 KO\nnon-genic")

contexts_plot(somatic[trt=="WT" & loc=="intergenic"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("Wildtype\nintergenic")
contexts_plot(somatic[trt=="MSH6"& loc=="intergenic"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("MSH6 KO\nintergenic")

contexts_plot(somatic[trt=="WT" & loc=="genic"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("Wildtype\ngenic")
contexts_plot(somatic[trt=="MSH6"& loc=="genic"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("MSH6 KO\ngenic")

contexts_plot(somatic[trt=="WT" & loc=="TE"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("Wildtype\nTE")
contexts_plot(somatic[trt=="MSH6"& loc=="TE"]$context_SBS,full=F)[[2]]+scale_y_continuous(name="Mutations")+ggtitle("MSH6 KO\nTE")

dev.off()

pdf("figures/MSH6_MAF_somatic.pdf", width=2, height=2)
somatic$genotype<-factor(ifelse(somatic$trt=="MSH6","msh6/msh6","wt/wt + wt/msh6"), levels=c("wt/wt + wt/msh6","msh6/msh6"))
ggplot(somatic ,aes(x=depth_pct, fill=genotype))+
  geom_histogram(col="black")+
  facet_wrap(~genotype)+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c("dodgerblue","orange"), guide="none")+
  ggtitle("Somatic")+
  scale_x_continuous(name="Alt frequency")+
  scale_y_continuous(name="Mutations")
dev.off()

pdf("figures/MSH6_MAF_all.pdf", width=3.5, height=2)
ggplot(PASS3 ,aes(x=depth_pct, fill=type))+
  geom_histogram()+
  facet_wrap(~trt)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Alt frequency")+
  scale_y_continuous(name="Mutations")
dev.off()

fisher.test(table(somatic$genic, somatic$trt))
fisher.test(table(germline$genic, germline$trt))

pdf("Figures/MSH6_genic_somatic.pdf", width=0.8, height=2)
fisher.test(table(somatic$genic, somatic$trt))
summary<-somatic[,.(genic_pct=sum(genic)/sum(!genic), genic=sum(genic), GC=sum(substr(mut,1,1)=="C")/.N, N=.N, transversion=sum(transversion)/.N, CtoA=sum(CtoA)/.N, depth=mean(tumor_depth), tumor_alt_depth=mean(tumor_alt_depth), depth_pct=mean(tumor_alt_depth/tumor_depth)), by=.(trt)][order(genic)]
ggplot(summary, aes(x=trt, y=(genic_pct), fill=trt))+geom_bar(stat="identity", col="black")+theme_classic(base_size = 6)+
  scale_fill_manual(values=c("dodgerblue","orange"),name="Genotype", guide="none")+
  scale_y_continuous(name="genic/non-genic mutations ")+
  scale_x_discrete(name="Genotype", labels=c("wt/wt + wt/msh6","msh6/msh6"))+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

pdf("Figures/MSH6_loc_somatic.pdf", width=1.5, height=2)
summary<-dcast(somatic[,.(.N), by=.(trt, loc)], loc~trt)
summary$MSH6_pct<-summary$MSH6/summary$WT
chisq.test(summary[,2:3])
ggplot(summary, aes(x=loc, y=(MSH6_pct), fill=loc))+
  geom_bar(stat="identity", col="black")+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c("gray20","gray35","gray90"), guide="none")+
  scale_y_continuous(name="", position="right")+
  scale_x_discrete(name="Location", labels=c("Gene\nbodies","Intergenic","TE"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))
dev.off()

pdf("Figures/MSH6_loc_germline.pdf", width=2, height=2)
summary<-dcast(germline[,.(.N), by=.(trt, loc)], loc~trt)
summary$MSH6_pct<-summary$MSH6/summary$WT
chisq.test(summary[,2:3])
ggplot(summary, aes(x=loc, y=(MSH6_pct), fill=loc))+
  geom_bar(stat="identity", col="black")+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c("gray20","gray35","gray90"), guide="none")+
  scale_y_continuous(name="msh6/wt\nmutation\nratio ")+
  scale_x_discrete(name="Location", labels=c("Gene\nbodies","Intergenic","TE"))+
  theme(axis.title.y = element_text(angle=0, hjust=0.5, vjust=0.5))
dev.off()


pdf("Figures/MSH6_genic_somatic_SALK_089638.pdf", width=0.8, height=2)
msh6<-somatic[!grepl("5.+|2L", file)]
msh6$genotype<-ifelse(grepl("2.+|1G",msh6$file), "wt/wt","msh6/msh6")
msh6$genotype<-factor(msh6$genotype, levels=c("wt/wt","msh6/msh6"))
fisher.test(table(msh6$genic, msh6$genotype!="msh6/msh6"))
summary<-msh6[,.(genic_pct=sum(genic)/sum(!genic), genic=sum(genic), GC=sum(substr(mut,1,1)=="C")/.N, N=.N, transversion=sum(transversion)/.N, CtoA=sum(CtoA)/.N, depth=mean(tumor_depth), tumor_alt_depth=mean(tumor_alt_depth), depth_pct=mean(tumor_alt_depth/tumor_depth)), by=.( genotype)][order(genic)]
ggplot(summary, aes(x=genotype, y=(genic_pct), fill=genotype))+geom_bar(stat="identity", col="black")+theme_classic(base_size = 6)+
  scale_fill_manual(values=c("dodgerblue","orange"),name="Genotype", guide="none")+
  scale_y_continuous(name="genic/non-genic mutations ")+
  scale_x_discrete(name="SALK_089638")+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

pdf("Figures/MSH6_genic_somatic_SALK_037557.pdf", width=1, height=2)
msh6_2<-somatic[grepl("5.+|2L", file)]
msh6_2$genotype<-ifelse(msh6_2$file=="5I_S16","msh6/msh6","wt/wt")
msh6_2$genotype[msh6_2$file=="5J_S17"]<-"wt/msh6"
msh6_2$genotype[msh6_2$file=="2L_S14"]<-"wt/msh6"
msh6_2$genotype<-factor(msh6_2$genotype, levels=c("wt/wt","wt/msh6","msh6/msh6"))
fisher.test(table(msh6_2$genic, msh6_2$genotype!="msh6/msh6"))
summary<-msh6_2[,.(genic_pct=sum(genic)/sum(!genic), genic=sum(genic), GC=sum(substr(mut,1,1)=="C")/.N, N=.N, transversion=sum(transversion)/.N, CtoA=sum(CtoA)/.N, depth=mean(tumor_depth), tumor_alt_depth=mean(tumor_alt_depth), depth_pct=mean(tumor_alt_depth/tumor_depth)), by=.( genotype)][order(genic)]
ggplot(summary, aes(x=genotype, y=(genic_pct), fill=genotype))+geom_bar(stat="identity", col="black")+theme_classic(base_size = 6)+
  scale_fill_manual(values=c("dodgerblue","gray", "orange"),name="Genotype", guide="none")+
  scale_y_continuous(name="genic/non-genic mutations ")+
  scale_x_discrete(name="SALK_037557")+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

windows<-make_feature_windows(encode_data = genes, 
                              2, gene = T)
setkey(windows, chr, start, stop)

  setkey(variants, CHROM, POSITION, POSITION2)
  overlaps <- foverlaps(windows, somatic[trt=="MSH6"], type = "any")
  overlaps<-overlaps[!is.na(POS)]
  sum<-overlaps[,.(N=.N, length=mean(length)), by=.(region, pos)]
  P1<-ggplot(sum, aes(x=pos, y=N/(length)))+
    geom_bar(stat="identity", fill="dodgerblue", col="black")+
    scale_y_continuous(name="Variants/b.p.")+
    theme_classic(base_size = 6)+
    scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")
  
  plot(P1)


WT<-plot_peaks(genes, "Genes, SBS","Gene body", 2,somatic[trt=="WT"], gene=T)[[2]];WT$geno="wt/wt + wt/msh6";WT$pct<-prop.table(WT$mut/WT$length)
MSH6<-plot_peaks(genes, "Genes, SBS","Gene body", 2,somatic[trt=="MSH6"], gene=T)[[2]];MSH6$geno="msh6/msh6";MSH6$pct<-prop.table(MSH6$mut/WT$length)

all<-rbind(WT, MSH6)

pdf("Figures/MSH6_gene_body.pdf", width=3.4, height=2)
sums<-all[,.(mut=sum(mut), sum=sum(length), pct=pct),by=.(pos,region, geno)]
sums$geno<-factor(sums$geno, levels=c("wt/wt + wt/msh6","msh6/msh6"))
fwrite(sums, "data/mutation_rate_gene_bodies_MSH6.csv")
ggplot(sums, aes(x=pos, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_line()+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c("dodgerblue","orange"), guide="none")+
  scale_y_continuous(name="Proportion of mutations")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")
ggplot(sums, aes(x=pos, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_line()+
  facet_wrap(.~geno)+
  scale_fill_manual(values=c("dodgerblue","orange"), guide="none")+
  scale_y_continuous(name="Proportion of mutations")+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")
dev.off()


chr<-rbindlist(lapply(list.files("data/Arabidopsis Epigenome Profiles/", full.names = T, pattern="*200*"), function(f){
  chr<-fread(f)
  chr$CHROM<-as.character(chr$chr)
  chr$START<-chr$start
  chr$STOP<-chr$stop
  return(chr)
}))
setkey(chr, CHROM, start, stop)


pdf("Figures/MSH6_KO_H3K4me1.pdf", width=1.5, height=2)
chr$MSH6<-mutations_in_features(chr, somatic[ trt=="MSH6"])
chr$WT<-mutations_in_features(chr, somatic[ trt=="WT"])
summary<-chr[,.(MSH6=sum(MSH6), WT=sum(WT), ratio=sum(MSH6)/sum(WT), H3K4me1_mean=mean(`H3K4me1_mean`), length=.N*200), by=.(H3K4me1=as.numeric(Hmisc::cut2(`H3K4me1_mean`, g = 3)))]
ggplot(summary, aes(x=H3K4me1, y=ratio, fill=as.character(H3K4me1)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=c("gray90","gray35","gray20"), guide="none")+
  scale_y_continuous(name="")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="H3K4me1 ChIP-seq", breaks=1:3, labels=c("33 %ile","66 %ile","100%ile"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))
chisq.test(summary[,2:3])
dev.off()


pdf("Figures/MSH6_KO_H3K4me1_germline.pdf", width=2, height=2)
chr$MSH6<-mutations_in_features(chr, germline[ trt=="MSH6"])
chr$WT<-mutations_in_features(chr, germline[ trt=="WT"])
summary<-chr[,.(MSH6=sum(MSH6), WT=sum(WT), ratio=sum(MSH6)/sum(WT), H3K4me1_mean=mean(`H3K4me1_mean`), length=.N*200), by=.(H3K4me1=as.numeric(Hmisc::cut2(`H3K4me1_mean`, g = 3)))]
ggplot(summary, aes(x=H3K4me1, y=ratio, fill=as.character(H3K4me1)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=c("gray90","gray35","gray20"), guide="none")+
  scale_y_continuous(name="msh6/wt\nmutation\nratio ")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="H3K4me1 ChIP-seq", breaks=1:3, labels=c("33 %ile","66 %ile","100%ile"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))
chisq.test(summary[,2:3])
dev.off()



#Fixed germline mutations (how many mutations are fixed between WT and the SALK_lines?)

results<-fread("data/strelka2_results.csv")
PASS<-results[Mut_N>0 & CHROM %in% 1:5 & FILTER=="PASS"]
PASS$unique3<-paste(PASS$unique, PASS$file, sep = "_")
tab<-data.table(table(PASS$unique3))
PASS$N<-tab$N[match(PASS$unique3, tab$V1)]
PASS$geno<-substr(PASS$file,1,1)
PASS$trt<-ifelse(PASS$geno=="2", "WT","MSH6")
PASS$trt[PASS$file %in% c("1G_S6","5J_S17","5B_S15")]<-"WT"

PASS$tumor_ref_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
PASS$tumor_alt_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
PASS$tumor_depth<-PASS$tumor_ref_depth+PASS$tumor_alt_depth
PASS$depth_pct<-PASS$tumor_alt_depth/PASS$tumor_depth

PASS$chi<-apply(PASS,1,  function(x){
  chisq.test(c(as.numeric(x["tumor_alt_depth"]),as.numeric(x["tumor_depth"])-as.numeric(x["tumor_alt_depth"])))$p.value
})

PASS$type<-ifelse(PASS$depth_pct>0.5 & PASS$chi<0.05,"homozygous","somatic")
PASS$type[PASS$chi>(0.001)]<-"heterozygous"

unique4<-unique(PASS[type=="homozygous" & N+Mut_N==17,c("CHROM","POS","REF","ALT","unique","trt","genic","file","Mut_N","POSITION","N"), with=F])

PASS3<-unique4

PASS3$geno<-ifelse(substr(PASS3$file,1,1)=="1", "WT","SALK_037557")
PASS3$geno[substr(PASS3$file,1,1)=="5"]<-"SALK_089638"
PASS3$geno[substr(PASS3$file,1,2)=="2L"]<-"SALK_089638"
PASS3$geno[substr(PASS3$file,1,2)=="1G"]<-"SALK_037557"
PASS3$unique_geno<-paste(PASS3$unique, PASS3$geno)

counts<-data.table(table(PASS3$unique_geno, PASS3$geno))[N>0]
fixed<-counts[,.(pct=N/max(N),N=N, unique=V1), by=.(geno=V2)][pct==1]
PASS3_fixed<-PASS3[unique_geno %in% fixed[]$unique]
table(PASS3_fixed$unique)
d<-dist(t(as.matrix(table(PASS3_fixed$unique, PASS3_fixed$file))), method = "manhattan")
clust<-hclust(d)
plot(clust)


#PASS vs no PASS genic?

results<-fread("data/strelka2_results.csv")
results$context_SBS<-contexts(vars = results, fasta = genome)
results$geno<-substr(results$file,1,1)
results$trt<-ifelse(results$geno=="2", "WT","MSH6")
results$trt[results$file %in% c("1G_S6","5J_S17","5B_S15")]<-"WT"
contexts_plot(results$context_SBS)
singles<-results[Mut_N==1 & FILTER=="PASS"]
contexts_plot(results$context_SBS, full = T)
table(singles$unique)

mod<-glm(FILTER=="PASS"~genic+trt+context_SBS, singles, family="binomial")
summary(mod)
plot(mod)

