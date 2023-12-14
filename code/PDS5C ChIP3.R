library(polymorphology2)
gff<-fread("data/A_thaliana_gff.txt")
colnames(gff)[1:9] <- c("CHROM","SOURCE","TYPE","START","STOP","SCORE","DIRECTION","PHASE","INFO")
CDS<-gff[CHROM %in% 1:5 & TYPE=="gene",.(CHROM, TYPE, START, STOP)]

# Create non-CDS regions
non_CDS <- CDS[, .(TYPE="non_CDS", 
                   START = data.table::shift(STOP, type = "lag") + 1, 
                   STOP = START - 1), 
               by = CHROM]
# Remove rows where START is NA or greater than STOP (i.e., the first row for each chromosome)
non_CDS <- non_CDS[!is.na(START) & START <= STOP]

CDS_non_CDS<-rbind(non_CDS, CDS)
CDS_non_CDS[,ID := 1:.N]

CDS_non_CDS_enrich1<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
                                         inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")
CDS_non_CDS_enrich2<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
                                          inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")

CDS_non_CDS$H3K4me1<-(CDS_non_CDS_enrich2$enrich+CDS_non_CDS_enrich1$enrich)/2
CDS_non_CDS$H3K4me1_input<-(CDS_non_CDS_enrich2$input+CDS_non_CDS_enrich1$input)


CDS_non_CDS_PDS51<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668646_RDM15_FLAG_ChIP_rep1_unique_reads.bedGraph",
                                          inputfile = "data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph")


CDS_non_CDS_PDS51.2<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4711900_RDM15_FLAG_rep1_uniq_more_reads.bedGraph",
                                        inputfile = "data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph")

CDS_non_CDS_PDS52<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668648_RDM15_FLAG_ChIP_rep2_unique_reads.bedGraph",
                                        inputfile = "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")


CDS_non_CDS_PDS52.2<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4711901_RDM15_FLAG_rep2_uniq_more_reads.bedGraph",
                                          inputfile = "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")


CDS_non_CDS$PDS5C<-(CDS_non_CDS_PDS51$enrich+CDS_non_CDS_PDS51.2$enrich+CDS_non_CDS_PDS52$enrich+CDS_non_CDS_PDS52.2$enrich)/4

CDS_non_CDS$PDS5_input<-CDS_non_CDS_PDS51$input+CDS_non_CDS_PDS51.2$input+CDS_non_CDS_PDS52$input+CDS_non_CDS_PDS52.2$input

#fwrite(CDS_non_CDS, "data/CDS_non_CDS_PDS5_H3K4me1.csv")
CDS_non_CDS<-fread("data/CDS_non_CDS_PDS5_H3K4me1.csv")

CDS_non_CDS_ranks<-CDS_non_CDS[PDS5_input>0 & H3K4me1_input>0,.(CDS=sum(TYPE=="gene"), non_CDS=sum(TYPE=="non_CDS"), N=.N), by=.(H3K4me1=as.numeric(cut(rank(H3K4me1), 100)), PDS5C=as.numeric(cut(rank(PDS5C), 100)))]
CDS_non_CDS_ranks$pct_CDS<-CDS_non_CDS_ranks$CDS/CDS_non_CDS_ranks$N

p1 <- ggplot(CDS_non_CDS_ranks, aes(x=(H3K4me1), y=(PDS5C), fill=pct_CDS, alpha=log(N))) +
  geom_tile() +
  #geom_point(shape=16, size=0.1) +
  scale_fill_gradient(high="dodgerblue", low="orange")+
  scale_alpha(range = c(0,1))+
  theme_minimal(base_size = 8) +
  coord_fixed() +
  labs(
    x = "H3K4me1 enrichment (%ile)",
    y = "PDS5C enrichment (%ile)",
    fill = "Percent genic"
  ) +
  theme(
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 4),          # Adjust text size
    legend.title = element_text(size = 5),         # Adjust title size
    legend.key.size = unit(0.25, "cm"),
    legend.position="none",# Adjust key size
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p1

pdf("figures/PDS5_H3K4me1_gene_intergene.pdf", width=1.75, height=1.75)
p1
dev.off()

genes<-fread("data/A_thal_genes_PDS5_enrich.csv")

genes_ranks<-genes[,.(essential=sum(lethal=="Lethal"), non_essential=sum(lethal=="Non-lethal"), N=.N), by=.(H3K4me1=as.numeric(cut(rank(enrich_H3K4me1), 10)), PDS5C=as.numeric(cut(rank(enrich), 10)))]
genes_ranks$pct_essential<-genes_ranks$essential/genes_ranks$N
summary(genes_ranks$pct_essential)

p3 <- ggplot(genes_ranks, aes(x=(H3K4me1)*10, y=(PDS5C)*10, color=log(pct_essential), size=(N), alpha=log(N))) +
  #geom_tile() +
  geom_point(shape=16) +
  scale_color_gradient(high="dodgerblue", low="orange")+
  scale_size(range=c(0,3), guide="none")+
  scale_alpha(range = c(0,1),guide="none")+
  theme_minimal(base_size = 8) +
  coord_fixed() +
  labs(
    x = "H3K4me1 enrichment (%ile)",
    y = "PDS5C enrichment (%ile)",
    fill = "Percent genic"
  ) +
  theme(
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 4),          # Adjust text size
    legend.title = element_text(size = 5),         # Adjust title size
    legend.key.size = unit(0.25, "cm"),
    legend.position="none",# Adjust key size
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p3

pdf("figures/PDS5_H3K4me1_essential.pdf", width=1.75, height=1.75)
p3
dev.off()


# Save the plot to a PDF
pdf("figures/PDS5_H3K4me1_relation.pdf", width=3.5, height=2)
print(p)
print(p2)
dev.off()

