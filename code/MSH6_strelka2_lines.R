source("code/libraries_functions.R")
library(polymorphology2)
library(gridExtra)


mutations<-fread("data/Strelka2_mutations.csv")
mutations$H3K4me1_enriched<-mutations$enrich_H3K4me1>1.5
genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genes[,ID:=1:nrow(genes)]
genes$DIRECTION<-genes$direction
mutations$tissue_all<-features_in_sites(genes[tissue_breadth==54], mutations)$overlaps

genome<-read.fasta("data/TAIR10_chr_all.fa")
names(genome)<-gsub("Chr","",names(genome))

hp<-make_homopolymer(genome, 3)

mutationshp<-homopolymer_var_annotate(mutations, hp, 3, 1)
mutations$HP<-mutationshp$HP

enrichment_line<-line_data(mutations[type=="somatic"])
enrichment_line$trt<-factor(enrichment_line$trt, levels=c("WT","MSH6"))

pdf("figures/line_boxplots2.pdf", width=1.1, height=2)
enrichment_line<-line_data(mutations[type=="somatic"])
grid.arrange(
  line_boxplot(enrichment_line, "genic_pct","% genic", title="Somatic"),
  line_boxplot(enrichment_line, "H3K4me1_pct","% H3K4me1 enriched", title="Somatic"),
  ncol=2)
enrichment_line<-line_data(mutations)

enrichment_line<-line_data(mutations[type!="somatic"])
grid.arrange(
  line_boxplot(enrichment_line, "genic_pct","% genic", title="Germline"),
  line_boxplot(enrichment_line, "H3K4me1_pct","% H3K4me1 enriched", title="Germline"),
  ncol=2)
enrichment_line<-line_data(mutations)

grid.arrange(
  line_boxplot(enrichment_line, "genic_pct","% genic", title="All"),
  line_boxplot(enrichment_line, "H3K4me1_pct","% H3K4me1 enriched", title="All"),
  ncol=2)
dev.off()

gene_windows<-feature_windows(genes, breaks = 2, dist = 2000, directed = T, IDcol = "gene")
gene_windows$wt<-sites_in_features(gene_windows, mutations[trt=="WT" ], mode="counts")$counts
gene_windows$msh6<-sites_in_features(gene_windows, mutations[trt=="MSH6"], mode="counts")$counts

# Melt the data
melted <- melt(gene_windows, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")

# Sum and calculate pct
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]

dcast_sums<-dcast(sums, RELATIVEPOS+LENGTH~geno, value.var = "sum")

pdf("Figures/MSH6_gene_body_all.pdf", width=3.5, height=2)
ggplot(sums, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_line()+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c("dodgerblue","orange"), guide="none")+
  scale_y_continuous(name="Prop")+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"), name="Region relative to gene bodies")

dev.off()

pdf("figures/MSH6_MAF_all.pdf", width=2, height=2)
mutations$trt<-factor(mutations$trt, levels=c("WT","MSH6"))
ggplot(mutations ,aes(x=depth_pct, fill=trt, alpha=type))+
  geom_histogram(bins=30, col="black", linewidth=0.25)+
  scale_alpha_manual(values = c(0.75, 0.5,  1), guide="none")+
  facet_wrap(~trt)+
  scale_fill_manual(values=c("dodgerblue","orange"), guide="none")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="VAF", breaks=c(0, 0.5, 1), labels=c("0","0.5","1"))+
  scale_y_continuous(name="Mutations")
dev.off()




chr<-rbindlist(lapply(list.files("data/Arabidopsis Epigenome Profiles/", full.names = T, pattern="*200*"), function(f){
  chr<-fread(f)
  chr$CHROM<-as.character(chr$chr)
  chr$START<-chr$start
  chr$STOP<-chr$stop
  return(chr)
}))
setkey(chr, CHROM, start, stop)


pdf("Figures/MSH6_KO_H3K4me1_mutations.pdf", width=1, height=1.5)
chr$MSH6<-mutations_in_features(chr, mutations[ trt=="MSH6"])
chr$WT<-mutations_in_features(chr, mutations[ trt=="WT"])
summary<-chr[,.(MSH6=sum(MSH6), WT=sum(WT), ratio=sum(MSH6)/sum(WT), H3K4me1_mean=mean(`H3K4me1_mean`), length=sum(stop-start)), by=.(H3K4me1=as.numeric(Hmisc::cut2(`H3K4me1_mean`, g = 3)))]
summary$MSH6_pct<-summary$MSH6/summary$length
summary$WT_pct<-summary$WT/summary$length
ggplot(summary, aes(x=H3K4me1, y=ratio, fill=as.character(H3K4me1)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=c("gray90","palegreen","green4"), guide="none")+
  scale_y_continuous(name="")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="H3K4me1 ChIP-seq", breaks=1:3, labels=c("33\n%ile","66\n%ile","100\n%ile"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))
chisq.test(summary[,2:3])

ggplot(summary, aes(x=H3K4me1, y=WT_pct, fill=as.character(H3K4me1)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=c("gray90","palegreen","green4"), guide="none")+
  scale_y_continuous(name="")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="H3K4me1 ChIP-seq", breaks=1:3, labels=c("33 %ile","66 %ile","100%ile"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))

ggplot(summary, aes(x=H3K4me1, y=MSH6_pct, fill=as.character(H3K4me1)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=c("gray90","palegreen","green4"), guide="none")+
  scale_y_continuous(name="")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="H3K4me1 ChIP-seq", breaks=1:3, labels=c("33 %ile","66 %ile","100%ile"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))

chisq.test(summary[,2:3])

dev.off()


pdf("Figures/MSH6_loc_mutations.pdf", width=1, height=1.5)
summary<-dcast(mutations[,.(.N), by=.(trt, loc)], loc~trt)
summary$loc<-factor(summary$loc, levels=c("TE","intergenic","genic"))
summary$MSH6_pct<-summary$MSH6/summary$WT
chisq.test(summary[,2:3])
ggplot(summary, aes(x=loc, y=(MSH6_pct), fill=loc))+
  geom_bar(stat="identity", col="black")+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c("gray90","palegreen","green4"), guide="none")+
  scale_y_continuous(name="", position="left")+
  scale_x_discrete(name="Location", labels=c("TE", "Inter\ngenic", "Gene\nbodies"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))
dev.off()




