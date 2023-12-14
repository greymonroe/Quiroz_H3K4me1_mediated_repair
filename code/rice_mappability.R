
mapp<-fread("data/rice_genome/genmap/Osativa_204_softmasked.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")

gene_annotations_all_mapp<-features_in_features(gene_annotations_all,mapp, mode = "mean", value = "MAPPABILITY" )

gene_annotations_all$MAPPABILITY<-gene_annotations_all_mapp$mean[match(gene_annotations_all$ID, gene_annotations_all_mapp$ID)]

gene_windows_mapp<-features_in_features(gene_windows,mapp, mode = "mean", value = "MAPPABILITY" )

gene_windows$MAPPABILITY<-gene_windows_mapp$mean[match(gene_windows_mapp$ID, gene_windows_mapp$ID)]

plot_feature_windows(gene_windows, variable="MAPPABILITY", mode="mean")

overlaps <- foverlaps(features, features2)

mapp$LENGTH<-mapp$STOP-mapp$START
sum(mapp[MAPPABILITY==1]$LENGTH)
sum(mapp$LENGTH)
mapp$ID<-1:nrow(mapp)

mutations <- ns_s[, .(CHROM = Chromosome, POS = Position)]

counts<-sites_in_features(mapp, mutations, mode="counts")
mapp$mutations<-counts$counts[match(mapp$ID, counts$ID)]

mapp_sum<-mapp[,.(mutations=sum(mutations), length=sum(LENGTH), rate=sum(mutations)/sum(LENGTH), pct=sum(LENGTH)/374246721, N=.N, MAPPABILITY=mean(MAPPABILITY)), by=.(mapp=cut2(MAPPABILITY, g=10))]

library(ggrepel)

pdf("figures/rice_mappability_scatter.pdf", width=7.5, height=1.5)

p1<-plot_feature_windows(gene_windows, variable="MAPPABILITY", mode="mean")$plot+scale_y_continuous(name="MAPPABILITY")+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "START", "STOP", "+5kb"))

p2<-ggplot(mapp_sum, aes(x=MAPPABILITY, y=rate))+
  geom_point(shape=21, fill="red")+
  scale_y_continuous(name="Mutations/bp")+
  theme_classic(base_size = 6)+
  scale_x_log10()+
  geom_text_repel(aes(label=paste(round(pct*100, 2), "%")), nudge_y = 0.000005, size=2)+
  theme(plot.background = element_blank(),
        panel.background = element_blank())

mapp_summ_genes<-gene_annotations_all[is.finite(enrich),.(mutations=sum(mutations), length=sum(length), enrich=mean(enrich), sum(mutations)/sum(length), N=.N, MAPPABILITY=mean(MAPPABILITY)), by=.(enrichment=cut2(enrich, g=10))]

p3<-ggplot(mapp_summ_genes, aes(x=enrich, y=MAPPABILITY))+
  geom_point(shape=21, fill="red")+
  scale_y_continuous(name="MAPPABILITY")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="log(H3K4me1/H3)")

mapp_summ_genes<-gene_annotations_all[is.finite(pnps) & is.finite(enrich),.(mutations=sum(mutations), length=sum(length), enrich=mean(enrich), sum(mutations)/sum(length), N=.N, MAPPABILITY=mean(MAPPABILITY)), by=.(PnPs=cut2(pnps, g=4))]

p4<-ggplot(mapp_summ_genes, aes(x=PnPs, y=MAPPABILITY))+
  geom_point(shape=21, fill="red")+
  scale_y_continuous(name="MAPPABILITY")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Pn/Ps 3000 Genomes")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

grid.arrange(p2, p1, p3, p4, ncol=4, widths=c(2,2, 2, 1.1))
dev.off()




gene_annotations_all_HP<-features_in_features(gene_annotations_all,homopolymers, mode = "counts" )

gene_annotations_all$HP<-gene_annotations_all_HP$counts[match(gene_annotations_all$ID, gene_annotations_all_HP$ID)]

gene_windows_HP<-features_in_features(gene_windows,homopolymers, mode = "counts" )
gene_windows$HP<-gene_windows_HP$counts[match(gene_windows$ID, gene_windows_HP$ID)]

gene_windows_muts_noHP<-sites_in_features(gene_windows, ns_s_filt[HP==0], mode="counts")
gene_windows$mutations_noHP<-gene_windows_muts_noHP$counts[match(gene_windows$ID, gene_windows_muts_noHP$ID)]

gene_windows_muts<-sites_in_features(gene_windows, ns_s_filt, mode="counts")
gene_windows$mutations<-gene_windows_muts$counts[match(gene_windows$ID, gene_windows_muts$ID)]
gene_windows$pctHP<-(gene_windows$mutations-gene_windows$mutations_noHP)/gene_windows$mutations*100

maxpos <- max(gene_windows$RELATIVEPOS)


pdf("figures/rice_HP_scatter.pdf", width=7.5, height=1.5)

gene_windows_filt<-gene_windows[model %in% genes_filt$model & MAPPABILITY==1]
p1<-plot_feature_windows(gene_windows_filt, variable="HP", mode="percent")$plot+scale_y_continuous(name="HP sequences/bp")+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "START", "STOP", "+5kb"))

gene_windows_filt<-gene_windows[model %in% genes_filt$model & MAPPABILITY==1]
p2<-plot_feature_windows(gene_windows_filt, variable="mutations_noHP", mode="percent")$plot+geom_line(col = "lightblue")+scale_y_continuous(name="Mutations/bp")+
  geom_line(data=plot_feature_windows(gene_windows_filt, variable="mutations", mode="percent")$summary, col="black")+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "START", "STOP", "+5kb"))

p3<-plot_feature_windows(gene_windows_filt, variable="pctHP", mode="mean")$plot+scale_y_continuous(name="% Homopolymer mutation")+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "START", "STOP", "+5kb"))

mapp_summ_genes<-genes_filt[is.finite(enrich),.(mutations=sum(mutations), length=sum(length), enrich=mean(enrich), sum(mutations)/sum(length), N=.N, HP=sum(HP), MAPPABILITY=mean(MAPPABILITY)), by=.(enrichment=cut2(enrich, g=10))]

p4<-ggplot(mapp_summ_genes, aes(x=enrich, y=HP/length))+
  geom_point(shape=21, fill="red")+
  scale_y_continuous(name="HP sequences/bp")+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="log(H3K4me1/H3)")

mapp_summ_genes<-genes_filt[is.finite(pnps) & is.finite(enrich),.(mutations=sum(mutations), length=sum(length), enrich=mean(enrich), sum(mutations)/sum(length), HP=sum(HP), N=.N, MAPPABILITY=mean(MAPPABILITY)), by=.(PnPs=cut2(pnps, g=4))]

p5<-ggplot(mapp_summ_genes, aes(x=PnPs, y=HP/length))+
  geom_point(shape=21, fill="red")+
  scale_y_continuous(name="HP sequences/bp")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Pn/Ps")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

grid.arrange(p4, p5, p1, p3,p2, ncol=5, widths=c(2,1.1,2, 2, 2))

dev.off()
