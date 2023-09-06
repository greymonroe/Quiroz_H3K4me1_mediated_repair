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


CDS_non_CDS_PDS51<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668646_RDM15_FLAG_ChIP_rep1_unique_reads.bedGraph",
                                          inputfile = "data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph")


CDS_non_CDS_PDS51.2<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4711900_RDM15_FLAG_rep1_uniq_more_reads.bedGraph",
                                        inputfile = "data/GSM4668645_RDM15_FLAG_Input_rep1_unique_reads.bedGraph")

CDS_non_CDS_PDS52<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4668648_RDM15_FLAG_ChIP_rep2_unique_reads.bedGraph",
                                        inputfile = "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")


CDS_non_CDS_PDS52.2<-features_chip_enrich(CDS_non_CDS, chipfile = "data/GSM4711901_RDM15_FLAG_rep2_uniq_more_reads.bedGraph",
                                          inputfile = "data/GSM4668647_RDM15_FLAG_Input_rep2_unique_reads.bedGraph")


CDS_non_CDS$PDS5C<-(CDS_non_CDS_PDS51$enrich+CDS_non_CDS_PDS51.2$enrich+CDS_non_CDS_PDS52$enrich+CDS_non_CDS_PDS52.2$enrich)/4

fwrite("data/PDS5_H3K4me1/")



p2 <- ggplot(CDS_non_CDS, aes(x=rank(H3K4me1), y=rank(PDS5C))) +
  geom_density2d_filled(aes(fill = after_stat(level))) +
  facet_grid(~TYPE) +
  theme_minimal(base_size = 6) +
  coord_fixed() +
  labs(
    x = "H3K4me1 enrichment",
    y = "PDS5C enrichment",
    fill = "Density"
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 4),          # Adjust text size
    legend.title = element_text(size = 5),         # Adjust title size
    legend.key.size = unit(0.25, "cm"),
    legend.position="none",# Adjust key size
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


genes<-fread("data/A_thal_genes_PDS5_enrich.csv")


p <- ggplot(genes, aes(x=rank(enrich_H3K4me1), y=rank(enrich))) +
  geom_density2d_filled(aes(fill = after_stat(level))) +
  facet_grid(~lethal, labeller = as_labeller(c(`Lethal` = "Essential", `Non-lethal` = "Non-Essential"))) +
  theme_minimal(base_size = 6) +
  coord_fixed() +
  labs(
    x = "H3K4me1 enrichment",
    y = "PDS5C enrichment",
    fill = "Density"
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 4),          # Adjust text size
    legend.title = element_text(size = 5),         # Adjust title size
    legend.key.size = unit(0.25, "cm"),
    legend.position="none",# Adjust key size
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the plot to a PDF
pdf("figures/PDS5_H3K4me1_relation.pdf", width=3.5, height=2)
print(p)
print(p2)
dev.off()

