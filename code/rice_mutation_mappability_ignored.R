
ns_s_filt<-ns_s[N<2]

# Adding specific mutation types counts
gene_annotations_all[, syn_mut := sites_in_features(gene_annotations_all, ns_s_filt[MutationType == "synonymous"], mode = "counts")$counts]
gene_annotations_all[, ns_mut := sites_in_features(gene_annotations_all, ns_s_filt[MutationType == "non-synonymous"], mode = "counts")$counts]
gene_annotations_all[, stop_mt := sites_in_features(gene_annotations_all, ns_s_filt[`StopCodon?` != ""], mode = "counts")$counts]

genes_filt <- gene_annotations_all[is.finite(enrich) & is.finite(pnps) & is_expressed == "Y"]

gene_windows_muts<-sites_in_features(gene_windows, ns_s_filt, mode="counts")
gene_windows$mutations<-gene_windows_muts$counts[match(gene_windows$ID, gene_windows_muts$ID)]

# PnPs lines --------------------------------------------------------------

pdf("figures/Lu_gene_lines_PnPs_no_MAPP_filt.pdf", height=1.75, width=2.3)

gene_sums<-genes_filt[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), stop_mt=sum(stop_mt), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  Pn=sum(Pn+stop_gained), Ps=sum(Ps), length=sum(LENGTH), meanpnps=mean(log(pnps))), by=.(PnPs=cut2(pnps, g=4))]

chisq.test(gene_sums[,.(ns_mut,syn_muts)])

ggplot(gene_sums, aes(x=PnPs, y=(ns_mut)/syn_muts, group=1))+
  geom_point(col="red")+
  geom_line(col="red")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(gene_sums$PnPs)))+
  scale_y_continuous(name="de novo MA N/S")+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45, hjust=1,vjust=1), plot.background = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(color="red"), axis.title.y = element_text(color="red"))

ggplot(gene_sums, aes(x=PnPs, y=enrich, group=1))+
  #geom_bar(stat="identity", fill="dodgerblue4", col="black", width=0.5)+
  geom_line(col="green3")+
  geom_point(shape=21, fill="green3")+
  theme_classic(base_size = 6)+
  geom_errorbar(aes(ymin=enrich-2*se_enrich, ymax=enrich+2*se_enrich), width=0, col="green3", alpha=0.5)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(gene_sums$PnPs)))+
  scale_y_continuous(name="log(H3K4me1/H3)")+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45, hjust=1,vjust=1), plot.background = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(color="green4"), axis.title.y = element_text(color="green4"))

genes_filt$PnPs<-as.numeric(cut2(genes_filt$pnps, g=4))
gene_sums<-genes_filt[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), Pn=sum(Pn), Ps=sum(Ps), meanpnps=mean(log(pnps)), CDS_mutation=sum(CDS_mutation)), by=.(PnPs=PnPs)][order(PnPs)]


gene_sums$silent<-(gene_sums$mut)/(gene_sums$length)
chisq.test(gene_sums[,.(mut, length)])

CI<-rbindlist(lapply(1:max(genes_filt$PnPs), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=genes_filt[PnPs==i,], statistic=bootstrap_stat, R=100)
  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  return(dt)
}))

# Add the lower and upper CI as new columns in your gene_sums dataframe
gene_sums$lower_ci <- CI[,2]
gene_sums$upper_ci <- CI[,1]

ggplot(gene_sums, aes(x=factor(PnPs), y=silent, group=1))+
  geom_point()+
  geom_line(col="black")+
  theme_classic(base_size = 6)+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0, alpha=0.5)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(factor(gene_sums$PnPs))))+
  scale_y_continuous(name="Mutations/bp", position = "right", breaks=c(c(8,8.5,9,9.5)*10^-5), labels = real_sci_format(c(8,8.5,9,9.5)*10^-5))+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45, hjust=1,vjust=1), plot.background = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(color="black"), axis.title.y = element_text(color="black"))

gene_sums$silent<-(gene_sums$mut-gene_sums$CDS_mutation)/(gene_sums$length-gene_sums$CDS)
chisq.test(gene_sums[,.(mut-CDS_mutation, length-CDS)])


ggplot(gene_sums, aes(x=factor(PnPs), y=silent, group=1))+
  geom_point()+
  geom_line(col="black")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(factor(gene_sums$PnPs))))+
  scale_y_continuous(name="UTR + Intron mutations/bp", position = "left")+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45, hjust=1,vjust=1), plot.background = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(color="black"), axis.title.y = element_text(color="black"))

chisq.test(gene_sums[,.(syn_muts, CDS)])

ggplot(gene_sums, aes(x=factor(PnPs), y=syn_muts/CDS, group=1))+
  geom_point()+
  geom_line(col="black")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(factor(gene_sums$PnPs))))+
  scale_y_continuous(name="Syn Mutations/CDS bp", position = "left")+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45, hjust=1,vjust=1), plot.background = element_blank(), panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(color="black"), axis.title.y = element_text(color="black"))
dev.off()


# H3K4me1 enrich scatter --------------------------------------------------


pdf("figures/Lu_genic_enrich_scatter_no_MAPP_filt.pdf", height=1.75, width=2)


genes_filt$enrichment<-as.numeric(cut2(genes_filt$enrich, g=10))
gene_sums<-genes_filt[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), CDS_mutation=sum(CDS_mutation), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), meanpnps=mean(log(pnps))), by=.(enrichment)][order(enrichment)]

gene_sums$silent<-(gene_sums$mut)/(gene_sums$length)
chisq.test(gene_sums[,.(mut, length)])
cor.test(gene_sums$silent, gene_sums$enrich)

# Perform the bootstrapping
CI<-rbindlist(lapply(1:max(genes_filt$enrichment), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=genes_filt[enrichment==i,], statistic=bootstrap_stat, R=100)
  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  return(dt)
}))

# Add the lower and upper CI as new columns in your gene_sums dataframe
gene_sums$lower_ci <- CI[,2]
gene_sums$upper_ci <- CI[,1]

ggplot(gene_sums, aes(x=enrich, y=silent, fill=as.numeric(enrichment)*10))+
  geom_smooth(method="lm", col="gray")+
  scale_fill_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green2", 
    midpoint = 50,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), alpha=0.2, width=0) +
  scale_x_continuous(name="log(H3K4me1/H3)")+
  theme_classic(base_size = 6)+
  theme(
    legend.position = c(.85, .75),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
    panel.background = element_blank(), 
    plot.background = element_blank()
  )+
  geom_point(shape=21, col="black")+
  scale_y_continuous(name="Mutations/bp")

chisq.test(gene_sums[,.(syn_muts, CDS)])
ggplot(gene_sums, aes(x=enrich, y=(syn_muts/CDS), fill=as.numeric(enrichment)*10))+
  geom_smooth(method="lm", col="gray")+
  scale_fill_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green2", 
    midpoint = 50,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  scale_x_continuous(name="log(H3K4me1/H3)")+
  theme_classic(base_size = 6)+
  theme(
    legend.position = c(.85, .75),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
    panel.background = element_blank(), 
    plot.background = element_blank()
  )+
  geom_point(shape=21, col="black")+
  scale_y_continuous(name="Syn Mutations/CDS bp")

gene_sums$silent<-(gene_sums$mut-gene_sums$CDS_mutation)/(gene_sums$length-gene_sums$CDS)
chisq.test(gene_sums[,.(mut-CDS_mutation, length-CDS)])

ggplot(gene_sums, aes(x=enrich, y=silent, fill=as.numeric(enrichment)*10))+
  geom_smooth(method="lm", col="gray")+
  scale_fill_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green2", 
    midpoint = 50,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  scale_x_continuous(name="log(H3K4me1/H3)")+
  theme_classic(base_size = 6)+
  theme(
    legend.position = c(.85, .75),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
    panel.background = element_blank(), 
    plot.background = element_blank()
  )+
  geom_point(shape=21, col="black")+
  scale_y_continuous(name="UTR + Intron mutations/bp")

chisq.test(gene_sums[,.(ns_mut, syn_muts)])

ggplot(gene_sums, aes(x=enrich, y=(ns_mut/syn_muts), fill=as.numeric(enrichment)*10))+
  geom_smooth(method="lm", col="gray")+
  scale_x_continuous(name="log(H3K4me1/H3)")+
  theme_classic(base_size = 6)+
  theme(
    legend.position = "none",
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
    panel.background = element_blank(), 
    plot.background = element_blank()
  )+
  geom_point(shape=21, col="black", fill="red")+
  scale_y_continuous(name="Non-Syn/Syn Mutations")


dev.off()


# Gene windows enrich -----------------------------------------------------


pdf("figures/Lu_gene_window_enrich_no_MAPP_filt.pdf", height=1.7, width=2.25)

gene_windows_filt<-gene_windows[model %in% gene_annotations_all[is.finite(pnps) & is.finite(enrich) & is_expressed=="Y"]$model]

gene_annotations_all$enrichment<-as.numeric(cut2(gene_annotations_all$enrich, g=4))

gene_windows_filt$enrichment<-gene_annotations_all$enrichment[match(gene_windows_filt$model, gene_annotations_all$model)]

maxpos <- max(gene_windows_filt$RELATIVEPOS)

enrich_summary <- gene_windows_filt[enrichment %in% c(1, 4), .(enrich = log(sum(H3K4me1)/sum(H3)), muts=sum(mutations), mut=sum(mutations)/sum(LENGTH),length=sum(LENGTH)), by = .(REGION, RELATIVEPOS, enrichment)]


ggplot(enrich_summary, aes(x = RELATIVEPOS, y = mut, col=as.numeric(enrichment), group=as.numeric(enrichment))) +
  #geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed") +
  theme_classic(base_size = 6) +
  geom_line()+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="Mutations/bp")+
  geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed", col="gray") +
  theme(panel.grid = element_blank(), legend.position = "none",panel.background = element_blank(), 
        plot.background = element_blank())+
  #geom_hline(yintercept = 0)+
  scale_color_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green3", 
    midpoint = 2.5,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "0%", "100%", "+5kb"), name="Gene body")


ggplot(enrich_summary, aes(x = RELATIVEPOS, y = enrich, col=as.numeric(enrichment), group=as.numeric(enrichment))) +
  geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed") +
  theme_classic(base_size = 6) +
  geom_line()+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="log(H3K4me1/H3)")+
  geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed", col="gray") +
  theme(panel.grid = element_blank(), legend.position = "none",panel.background = element_blank(), 
        plot.background = element_blank())+
  #geom_hline(yintercept = 0)+
  scale_color_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green3", 
    midpoint = 2.5,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "0%", "100%", "+5kb"), name="Gene body")


ggplot(enrich_summary, aes(x=enrich, y=mut, col=as.numeric(enrichment),shape=REGION, group=factor(enrichment)))+
  scale_color_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green3", 
    midpoint = 2.5,
    guide = "none"
  )+
  geom_smooth(method = "lm", col="gray")+
  geom_point()+
  theme_classic(base_size = 6)+
  theme(legend.position = "none",panel.background = element_blank(), 
        plot.background = element_blank())+
  scale_x_continuous(name="log(H3K4me1/H3)")+
  scale_y_continuous(name="Mutations/bp")

dev.off()

# gene windows enrich simple ----------------------------------------------



pdf("figures/Lu_gene_window_enrich_simple_no_MAPP_filt.pdf", height=1.7, width=2.25)

gene_windows_filt<-gene_windows[model %in% genes_filt$model ]

enrich_summary_region <- gene_windows_filt[, .(enrich = log(sum(H3K4me1)/sum(H3)), muts=sum(mutations), mut=sum(mutations)/sum(LENGTH),length=sum(LENGTH)), by = .(REGION=="gene body")]

maxpos <- max(gene_windows$RELATIVEPOS)
enrich_summary <- gene_windows_filt[, .(enrich = log(sum(H3K4me1)/sum(H3)), muts=sum(mutations), mut=sum(mutations)/sum(LENGTH),length=sum(LENGTH)), by = .(REGION, RELATIVEPOS)]

ggplot(enrich_summary, aes(x = RELATIVEPOS, y = mut)) +
  geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed", col="gray") +
  theme_classic(base_size = 6) +
  geom_line()+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="Mutations/bp", position = "right")+
  theme(panel.grid = element_blank(), panel.background = element_blank(), plot.background = element_blank(), axis.line.y = element_blank())+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "0%", "100%", "+5kb"), name="Gene body")

ggplot(enrich_summary, aes(x = RELATIVEPOS, y = enrich)) +
  theme_classic(base_size = 6) +
  geom_line(col="green3")+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="log(H3K4me1/H3)")+
  theme(panel.grid = element_blank(), panel.background = element_blank(), plot.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(colour = "green4"), axis.title.y = element_text(color="green4"))+
  scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-5kb", "0%", "100%", "+5kb"), name="Gene body")

dev.off()





# H3K4me1 histogram -------------------------------------------------------

pdf("figures/Lu_genic_histogram_enrich_no_MAPP_filt.pdf", height=1.75, width=2.5)


# Calculate the 0.5th and 99.5th percentiles of the enrich values
lower_bound <- quantile(genes_filt$enrich, probs = 0.005, na.rm = TRUE)
upper_bound <- quantile(genes_filt$enrich, probs = 0.995, na.rm = TRUE)

# Subset the genes_filt data frame to keep only the rows where enrich values are within the 99th percentile range
genes_filt_99 <- genes_filt[genes_filt$enrich >= lower_bound & genes_filt$enrich <= upper_bound, ]

ggplot(genes_filt_99, aes(x=enrich, group=as.numeric(cut2(enrich,g=10))*10, fill=as.numeric(cut2(enrich,g=10))*10))+
  geom_histogram(bins=50, col="black", linewidth=0.25 )+
  scale_fill_gradient2(
    name = "%ile", 
    low = "orange2", 
    mid = "gray", 
    high = "green2", 
    midpoint = 50,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  scale_x_continuous(name="log(H3K4me1/H3)")+
  theme_classic(base_size = 6)+
  theme(
    legend.position = c(.3, .7),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
    panel.background = element_blank(), 
    plot.background = element_blank()
  )+
  scale_y_continuous(name="# genes")

dev.off()
