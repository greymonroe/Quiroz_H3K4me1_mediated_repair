library(data.table)
library(seqinr)
library(ggplot2)
library(polymorphology)
source("code/libraries_functions.R")


gff<-fread("data/A_thaliana_gff.txt")
genes<-gff[type=="gene"]
genes$START<-genes$start
genes$CHROM<-as.character(genes$chr)
genes$STOP<-genes$stop
mutations<-fread("data/A_thaliana_germline_mutations.txt")[!is.na(POS)]
mutations$CHROM<-as.character(mutations$CHROM)
mutations$POSITION<-mutations$POS
genes$mutations<-mutations_in_features(genes, mutations)

CDS<-fread("data/A_thaliana_genes.txt")
CDS$chr<-as.character(CDS$chr)
genes<-merge(CDS,genes)
setkey(genes, chr, start, stop)

# data from Niu et al. 2021, Nature Communications:
infiles<-c("data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph", 
           "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")
chipfiles<-c("data/GSM4668646_RDM15_FLAG_ChIP_rep1_unique_reads.bedGraph",
             "data/GSM4668648_RDM15_FLAG_ChIP_rep2_unique_reads.bedGraph", 
             "data/GSM4711900_RDM15_FLAG_rep1_uniq_more_reads.bedGraph",
             "data/GSM4711901_RDM15_FLAG_rep2_uniq_more_reads.bedGraph")

inputs<-lapply(infiles, function(i){ 
  chip_overlaps(i, genes)
})

chips1<-lapply(chipfiles[c(1,3)], function(i){ 
  chip_overlaps(bedfile = i, featureobject = genes)
})

chips2<-lapply(chipfiles[c(2,4)], function(i){ 
  chip_overlaps(i, genes)
})

in1_total<-chip_total(infiles[1])
in2_total<-chip_total(infiles[1])
chip1_total<-chip_total(chipfiles[1])+chip_total(chipfiles[3])
chip2_total<-chip_total(chipfiles[2])+chip_total(chipfiles[4])

genes$chip1<-rowSums(as.data.table(chips1), na.rm=F)
genes$input1<-rowSums(as.data.table(inputs)[,c(1)], na.rm=F)
genes$enrich1<-log2((1+genes$chip1)/chip1_total) - log2((1+genes$input1)/in1_total)

genes$chip2<-rowSums(as.data.table(chips2), na.rm=F)
genes$input2<-rowSums(as.data.table(inputs)[,c(2)], na.rm=F)
genes$enrich2<-log2((1+genes$chip2)/chip2_total) - log2((1+genes$input2)/in2_total)
genes$enrich<-rowMeans(genes[,c("enrich1","enrich2"), with=F], na.rm=T)

fwrite(genes, "data/A_thal_genes_PDS5_enrich.csv")

pdf("figures/PDS5C_chip_enrich_constraint.pdf", width=1.5, height=1.5)
means_cds<-genes[is.finite(NI),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich2, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(NI, g=20)))]
cor<-cor.test(genes[is.finite(NI)]$NI, genes[is.finite(NI)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp*5, y=enrich))+
  geom_line(linewidth=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="NI %ile")

means_cds<-genes[,.(mutations=sum(mutations, na.rm=T), length=sum(stop-start), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(enrich, g=20)))]
cor<-cor.test((genes$mutations)/(genes$stop-genes$start), genes$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=mutations/length))+
  geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="Mutations/b.p.")+
  #geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="PDS5C enrichment")

means_cds<-genes[is.finite(TajimasD),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich2, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(TajimasD, g=20)))]
cor<-cor.test(genes[is.finite(TajimasD)]$TajimasD, genes[is.finite(TajimasD)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp*5, y=enrich))+
  geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Tajimas D %ile")

means_cds<-genes[is.finite(tissue_breadth),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=(Hmisc::cut2(tissue_breadth, g=10)))]
cor<-cor.test(genes[is.finite(tissue_breadth)]$tissue_breadth, genes[is.finite(tissue_breadth)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=enrich))+
  geom_line(size=0.25, col="green4", group=1)+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_discrete(name="Tissue breadth")

means_cds<-genes[is.finite(DnDs),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(DnDs, g=100)))]
cor<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=enrich))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Dn/Ds %ile")

means_cds<-genes[is.finite(PnPs),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(PnPs, g=100)))]
cor<-cor.test(genes[is.finite(PnPs)]$PnPs, genes[is.finite(PnPs)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp*5, y=enrich))+
  geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Pn/Ps %ile")

means_cds<-genes[is.finite(H3K4me1),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(H3K4me1, g=100)))]
cor<-cor.test(genes$H3K4me1, genes$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=enrich))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="H3K4me1 %ile")

means_cds<-genes[is.finite(H3K4me1)&is.finite(DnDs),.(DnDs=mean(DnDs, na.rm=T), length=sum(`CDS length`), se=sd(DnDs, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(H3K4me1, g=100)))]
cor<-cor.test(genes[is.finite(H3K4me1)&is.finite(DnDs)]$H3K4me1, genes[is.finite(H3K4me1)&is.finite(DnDs)]$DnDs)
ggplot(means_cds[!is.na(grp)], aes(x=grp*5, y=DnDs))+
  geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="DnDs")+
  geom_errorbar(aes(ymin=DnDs-se, ymax=DnDs+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="H3K4me1 %ile")

means_cds<-genes[is.finite(DnDs),.(enrich=mean(H3K4me1, na.rm=T), length=sum(`CDS length`), se=sd(H3K4me1, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(DnDs, g=100)))]
cor<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$H3K4me1)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=enrich))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="H3K4me1")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Dn/Ds %ile")

dev.off()

summary(lm(enrich~H3K4me2+H3K4me1+H3K4me3+H3K27me1+H3K27ac+H3K36ac+H3K36me3+H3K9me1+H3K56ac+H3K9me2+H3K4me3, genes))

pdf("figures/PDS5C_chip_enrich_ESN.pdf", width=1.5, height=1.5)
  ESN<-fread("data/A_thaliana_essential_genes.txt")
  ESN[!ESN$gene %in% genes$gene]
  ESN_genes<-merge(ESN, genes, by="gene")
  means_esn<-ESN_genes[,.(enrich=mean(enrich, na.rm=T), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp)]
  means_esn$grp<-factor(means_esn$grp, levels=means_esn$grp[order(means_esn$enrich)])
  summary(lm(enrich~grp=="ESN",ESN_genes ))
  ggplot(means_esn, aes(x=grp, y=enrich))+
    geom_point(size=0.5, col="green4")+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="PDS5C\nenrichment")+
    geom_errorbar(aes(ymin=enrich-2*se, ymax=enrich+2*se), width=0, col="green4")+
    scale_x_discrete(name="Gene function", labels=c("Cellular","Environmental","Morphological","Essential"))+
    theme(axis.text.x = element_text(angle=90, hjust=1))
dev.off()


germs<-fread("data/A_thaliana_regions.txt")
germs$chr<-as.character(germs$chr)
germs$gene<-germs$info
setkey(germs, chr, start, stop)

pdf("figures/A_thal_features_mutation_germline.pdf", width=1.25, height=1.25)

type_means<-germs[,.(mut=sum(MA_snp), length=sum(length), pct=sum(MA_snp)/sum(length)), by=.(grp=type)]
type_means$grp<-c("Upstream","5' UTR","Coding","Intron","Downstream", "3' UTR")
type_means$grp<-factor(type_means$grp, levels=c("Upstream","5' UTR","Intron","Coding","3' UTR","Downstream"))
chisq.test(type_means[,2:3])

ggplot(type_means, aes(x=grp, y=pct/(107*25+400*8)))+
  geom_bar(stat="identity", size=0.5)+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="Mutations/bp")+
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())+
  scale_x_discrete()
dev.off()

inputs<-lapply(infiles, function(i){ 
  chip_overlaps(i, germs)
})

chips1<-lapply(chipfiles[c(1,3)], function(i){ 
  chip_overlaps(i, germs)
})

chips2<-lapply(chipfiles[c(2,4)], function(i){ 
  chip_overlaps(i, germs)
})

germs$chip1<-rowSums(as.data.table(chips1), na.rm=F)
germs$input1<-rowSums(as.data.table(inputs)[,c(1)], na.rm=F)
germs$enrich1<-log2((1+germs$chip1)/sum(germs$chip1, na.rm=T)) - log2((1+germs$input1)/sum(germs$input1, na.rm=T))

germs$chip2<-rowSums(as.data.table(chips2), na.rm=F)
germs$input2<-rowSums(as.data.table(inputs)[,c(2)], na.rm=F)
germs$enrich2<-log2((1+germs$chip2)/sum(germs$chip2, na.rm=T)) - log2((1+germs$input2)/sum(germs$input2, na.rm=T))
germs$enrich<-rowMeans(germs[,c("enrich1","enrich2"), with=F], na.rm=T)

means<-germs[,.(mut=sum(MA_snp), length=sum(length), pct=sum(MA_snp)/sum(length), enrich=mean(enrich, na.rm=T), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(enrich, g=20)))]
chisq.test(means[,2:3])

pdf("figures/PDS5C_chip_enrich_mutationrate.pdf", width=1.5, height=1.5)
ggplot(means, aes(x=grp*5, y=pct/(107*25+400*8)))+
  geom_point(size=0.5, col="green4")+
  theme_classic(base_size = 6)+
  geom_line(size=0.25, col="green4")+
  scale_y_continuous(name="A. thaliana\ngermline\nmutation rate")+
  scale_x_continuous(name="PDS5C enrichment\n%ile")
dev.off()

pdf("figures/PDS5C_chip_enrich_CDS.pdf", width=.75, height=1.25)

  type_means<-germs[,.(mut=sum(MA_snp), length=sum(length), pct=sum(MA_snp)/sum(length), enrich=mean(enrich, na.rm=T), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=type)]
  type_means$grp<-c("Upstream","5' UTR","Coding","Intron","Downstream", "3' UTR")
  type_means$grp<-factor(type_means$grp, levels=c("Upstream","5' UTR","Intron","Coding","3' UTR","Downstream"))
  TukeyHSD(aov(enrich~type, germs))
  ggplot(type_means, aes(x=grp, y=enrich))+
    #geom_point(size=1)+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="PDS5C  enrichment")+
    geom_errorbar(aes(ymin=enrich-2*se, ymax=enrich+2*se), width=0.5, col="green4")+
    theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank())+
    scale_x_discrete()
  
  cds_intron_means<-germs[type %in% c("CDS","intron"),.(mut=sum(MA_snp), length=sum(length), pct=sum(MA_snp)/sum(length), enrich=mean(enrich, na.rm=T), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=type)]
  ggplot(cds_intron_means, aes(x=grp, y=enrich))+
   # geom_point(size=1)+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="PDS5C enrich")+
    geom_errorbar(aes(ymin=enrich-2*se, ymax=enrich+2*se), width=0.1, col="green4")+
    theme(axis.text.x = element_text(angle=90, hjust=1))

  cds_intron_means<-germs[,.(mut=sum(MA_snp), length=sum(length), pct=sum(MA_snp)/sum(length), enrich=mean(enrich, na.rm=T), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=type=="CDS")]
  ggplot(cds_intron_means, aes(x=grp, y=enrich))+
    # geom_point(size=1)+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="PDS5C enrich")+
    geom_errorbar(aes(ymin=enrich-2*se, ymax=enrich+2*se), width=0.1, col="green4")+
    theme(axis.text.x = element_text(angle=90, hjust=1))
  
dev.off()

pdf("figures/PDS5C_chip_enrich_mutation.pdf", width=1.5, height=1.5)
genes$length<-genes$stop-genes$start
means_cds<-genes[is.finite(PnPs)&mutations>0,.(mutations=sum(mutations, na.rm=T), length=sum(stop-start), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(PnPs, g=10)))]
cor<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=mutations/length))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="Mutations/b.p.")+
  #geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Dn/Ds %ile")

means_cds<-genes[is.finite(DnDs)&mutations>0,.(mutations=sum(mutations, na.rm=T), length=sum(stop-start), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(H3K4me1, g=100)))]
cor<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$enrich)
ggplot(means_cds[!is.na(grp)], aes(x=grp, y=mutations/length))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="mutation rates")+
  #geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="H3K4me1")

means_cds<-genes[,.(mutations=sum(mutations, na.rm=T), length=sum(stop-start), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(enrich, g=100)))]
cor<-cor.test((genes$mutations)/genes$length, genes$enrich)
spearman<-cor.test((genes$mutations)/genes$length, genes$enrich, method = "spearman")

ggplot(means_cds[!is.na(grp)], aes(x=grp, y=log(mutations/length)))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="Mutations/b.p.")+
  #geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="PDS5C enrichment")

means_cds<-genes[is.finite(DnDs),.(enrich=mean(enrich, na.rm=T), length=sum(`CDS length`), se=sd(enrich, na.rm=T)/sqrt(.N)), by=.(grp=as.numeric(Hmisc::cut2(DnDs, g=100)))]
cor<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$enrich)
spear<-cor.test(genes[is.finite(DnDs)]$DnDs, genes[is.finite(DnDs)]$enrich, method="spearman")

ggplot(means_cds[!is.na(grp)], aes(x=grp, y=enrich))+
  #geom_line(size=0.25, col="green4")+
  geom_point(size=0.5)+
  theme_classic(base_size = 6)+
  ggtitle(paste("r =", round(cor$estimate, digits = 2)))+
  scale_y_continuous(name="PDS5C enrichment")+
  geom_errorbar(aes(ymin=enrich-se, ymax=enrich+se), width=0)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_x_continuous(name="Dn/Ds %ile")

dev.off()



