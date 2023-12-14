

# load  --------------------------------------------------------------

source("code/libraries_functions.R")
library(polymorphology2)
source("code/load_rice_data.R")

# gene windows regression ------------------------------------------------

model<-glm(mut~H3K4me1_dep + H3K36me3_dep + H3K27me3_dep + H3K56ac_dep + H3K4me3_dep+H3_dep+Input_dep, windows500_filt, family="binomial")
summary(model)

# model<-glm(mut~H3K4me1 + H3K36me3 + H3K27me3 + H3K56ac + H3K4me3+H3, windows200_filt, family="binomial")
# summary(model)

# Extract coefficients, z-values, and p-values
coefs <- summary(model)$coefficients
coefs <- as.data.frame(coefs)
coefs$Predictor <- rownames(coefs)

# Filter out unwanted predictors and remove "_dep" from names
coefs <- coefs[!coefs$Predictor %in% c("(Intercept)", "H3_dep", "Input_dep"), ]
coefs$Predictor <- gsub("_dep", "", coefs$Predictor)

# Rank by z-value
coefs <- coefs[order(coefs$`z value`), ]
coefs$Predictor <- factor(coefs$Predictor, levels = coefs$Predictor)

# Add multivariate model type to coefs
coefs$model <- "multivariate"

# List of predictors for univariate models
predictors <- c("H3K4me1_dep", "H3K36me3_dep", "H3K27me3_dep", "H3K56ac_dep", "H3K4me3_dep")

# Extract z-scores and p-values from each univariate model and add to coefs
for (predictor in predictors) {
  formula <- as.formula(paste("mut ~", predictor, "+ H3_dep + Input_dep"))
  model_uni <- glm(formula, data = windows500_filt, family = "binomial")
  coefs_uni <- as.data.frame(summary(model_uni)$coefficients)
  coefs_uni$Predictor <- gsub("_dep", "", rownames(coefs_uni))
  coefs_uni <- coefs_uni[coefs_uni$Predictor == gsub("_dep", "", predictor), , drop = FALSE]
  coefs_uni$model <- "univariate"
  coefs <- rbind(coefs, coefs_uni)
}

coefs<-data.table(coefs)

pdf("figures/Lu_2019_mutation_model_windows2.pdf", width=2.5, height=1.5)
p <- ggplot(coefs, aes(y = Predictor, x = `z value`, fill = -log10(`Pr(>|z|)`), color = -log10(`Pr(>|z|)`), group=model, linetype=model)) +
  geom_segment(data=coefs[model=="multivariate"], aes(y=as.numeric(Predictor)-.1, yend=as.numeric(Predictor)-.1, x=0, xend=`z value`))+
  geom_segment(data=coefs[model=="univariate"], aes(y=as.numeric(Predictor)+.1, yend=as.numeric(Predictor)+.1, x=0, xend=`z value`), linetype="dashed")+
  geom_point(data=coefs[model=="multivariate"], aes(y=as.numeric(Predictor)-.1), shape=21, col="black")+
  geom_point(data=coefs[model=="univariate"], aes(y=as.numeric(Predictor)+.1), shape=21, col="black")+
  labs(y = "Histone modification", x = "z-score", fill = "-log10(p)")+
  scale_fill_gradient(low = "gray", high = "green4", 
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_color_gradient(low = "gray", high = "green4", 
                       guide="none")+
  theme_classic(base_size = 6) +
  scale_x_continuous(limits=c(-max(abs(coefs$`z value`)), max(abs(coefs$`z value`))))+
  theme(legend.key.size = unit(0.25, "cm"), 
        legend.position = c(0.8,0.5),
        legend.background = element_blank(),
        panel.grid = element_blank(),panel.background = element_blank(), plot.background = element_blank())+
  scale_y_continuous(labels=unique(coefs$Predictor))+
  geom_vline(xintercept = 0)
  p

dev.off()


# CDS Regression ---------------------------------------------------------
# 
# 
# gff<-read.GFF("data/all.clean.gff")
# gff$name=substr(gff$INFO, 4,17)
# gff$model=gsub("(ID=)|:.+|;.+","", gff$INFO)
# cds<-gff[TYPE=="CDS"]
# cds[, ID := 1:nrow(cds)]
# cds<-cds[CHROM %in% paste0("Chr",1:12)]
# cds$LENGTH<-cds$STOP-cds$START
# 
# 
# 
# # List of bedGraph files and their corresponding names
# marks <- list(
#   H3K36me3 = "data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph",
#   H2A.Z = "data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph",
#   H3K56ac = "data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph",
#   H3K4me3 = "data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph",
#   Input = "data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph",
#   H3 = "data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph",
#   H3K27me3 = "data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph",
#   H3K4me1 = "data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph"
# )
# 
# # Read in the bedGraph files and compute the values
# for (mark in names(marks)) {
#   message(mark)
#   bedGraph_data <- read.bedGraph(marks[[mark]])
#   mark_data <- features_in_features(cds, bedGraph_data, mode="sumxlength", value = "DEPTH")
#   cds[[mark]]<-mark_data$sumxlength[match(cds$ID, mark_data$ID)]
# }
# 
# mutations <- ns_s[MutationType == "synonymous"][, .(CHROM = Chromosome, POS = Position)]
# counts<-sites_in_features(cds, mutations, mode="counts")
# cds$mutations<-counts$counts[match(cds$ID, counts$ID)]
# cds$LENGTH<-cds$STOP-cds$START
# 
# cds_unique <- unique(cds[, -c("INFO","ID","model"), with = FALSE])
# 
# cds_summary <- cds_unique[, .(
#   mutation_rate = sum(mutations) / sum(LENGTH),
#   H3K36me3 = sum(H3K36me3) / sum(LENGTH),
#   H2A.Z = sum(H2A.Z) / sum(LENGTH),
#   H3K56ac = sum(H3K56ac) / sum(LENGTH),
#   H3K4me3 = sum(H3K4me3) / sum(LENGTH),
#   Input = sum(Input) / sum(LENGTH),
#   H3 = sum(H3) / sum(LENGTH),
#   H3K27me3 = sum(H3K27me3) / sum(LENGTH),
#   H3K4me1 = sum(H3K4me1) / sum(LENGTH),
#   LENGTH=sum(LENGTH),
#   mutations=sum(mutations)
# ), by = .(gene = name)]
# 
# model<-lm(mutation_rate~H3K4me1 + H3K36me3 + H3K27me3 + H3K56ac + H3K4me3+H3+H2A.Z+Input, cds_summary)
# summary(model)
# 
# # Extract coefficients, t-values, and p-values
# coefs <- summary(model)$coefficients
# coefs <- as.data.frame(coefs)
# coefs$Predictor <- rownames(coefs)
# 
# # Filter out unwanted predictors and remove "_dep" from names
# coefs <- coefs[!coefs$Predictor %in% c("(Intercept)", "H3", "Input"), ]
# coefs$Predictor <- gsub("_dep", "", coefs$Predictor)
# 
# # Rank by z-value
# coefs <- coefs[order(coefs$`t value`), ]
# coefs$Predictor <- factor(coefs$Predictor, levels = coefs$Predictor)
# coefs$model <- "multivariate"
# 
# # List of predictors for univariate models
# predictors <- c("H3K4me1", "H3K36me3", "H3K27me3", "H3K56ac", "H3K4me3")
# 
# # Extract t-scores and p-values from each univariate model and add to coefs
# for (predictor in predictors) {
#   formula <- as.formula(paste("mutation_rate ~", predictor, "+ H3 + Input"))
#   model_uni <- lm(formula, data = cds_summary)
#   coefs_uni <- as.data.frame(summary(model_uni)$coefficients)
#   coefs_uni$Predictor <- gsub("_dep", "", rownames(coefs_uni))
#   coefs_uni <- coefs_uni[coefs_uni$Predictor == gsub("_dep", "", predictor), , drop = FALSE]
#   coefs_uni$model <- "univariate"
#   coefs <- rbind(coefs, coefs_uni)
# }
# 
# # Create the plot
# pdf("figures/Lu_2019_mutation_model_synonymous_CDS.pdf", width=3.5, height=1.5)
# p <- ggplot(coefs, aes(x = Predictor, y = `t value`, fill = -log10(`Pr(>|t|)`))) +
#   geom_bar(stat = "identity", position = "dodge", col = "black") +
#   scale_fill_gradient(low = "white", high = "green4") +
#   labs(x = "Predictor", y = "t-score", fill = "-log10(p)") +
#   theme_classic(base_size = 6) +
#   theme(legend.key.size = unit(0.25, "cm"), 
#         axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
#   facet_wrap(~ model, ncol = 2)
# 
# print(p)
# dev.off()
# 


# gene load and prepare data ----------------------------------------------

#genes_filt<-gene_annotations_all[is.finite(enrich) & is.finite(pnps) & MAPPABILITY==1 & is_expressed=="Y"]
#genes_filt$pnps_cut<-cut2(genes_filt$pnps, g=5)




# PnPs lines --------------------------------------------------------------

pdf("figures/Lu_gene_lines_PnPs.pdf", height=1.75, width=2.3)

gene_sums<-genes_filt[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), stop_mt=sum(stop_mt), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  Pn=sum(Pn+stop_gained), Ps=sum(Ps), length=sum(LENGTH), meanpnps=mean(log(pnps))), by=.(PnPs=cut2(pnps, g=4))]

chisq.test(gene_sums[,.(ns_mut,syn_muts)])

ggplot(gene_sums, aes(x=PnPs, y=(ns_mut)/syn_muts, group=1))+
  geom_point(col="red")+
  geom_line(col="red")+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="3000 genomes natural Pn/Ps", limits = rev(levels(gene_sums$PnPs)))+
  scale_y_continuous(name="Non-Syn/Syn Mutations")+
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


pdf("figures/Lu_genic_enrich_scatter.pdf", height=1.75, width=2)


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


pdf("figures/Lu_gene_window_enrich.pdf", height=1.7, width=2.25)

gene_windows_filt<-gene_windows[ !is.na(MAPPABILITY) & model %in% gene_annotations_all[is.finite(pnps) & is.finite(enrich) & is_expressed=="Y"]$model & MAPPABILITY==1]

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



pdf("figures/Lu_gene_window_enrich_simple.pdf", height=1.7, width=2.25)

gene_windows_filt<-gene_windows[model %in% genes_filt$model & MAPPABILITY==1]

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

pdf("figures/Lu_genic_histogram_enrich.pdf", height=1.75, width=2.5)


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

