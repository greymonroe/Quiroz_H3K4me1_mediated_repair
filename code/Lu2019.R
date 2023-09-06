source("code/libraries_functions.R")
library(polymorphology2)

chip_overlaps_Lu<-function(bedfile, featureobject){
  cat("\nreading ");cat(bedfile)
  in1<-fread(bedfile)
  colnames(in1)<-c("chr","start","stop","depth")
  in1$length<-as.numeric(in1$stop-in1$start)
  in1$depth<-as.numeric(in1$depth)
  setkey(in1, chr, start, stop)
  out<-rbindlist(lapply(unique(in1$chr), function(c) {
    cat(c)
    CDS_input_overlap<-foverlaps(featureobject[chr==c], in1[chr==c],type="any")
    CDS_input<-CDS_input_overlap[,.(len=sum(length,na.rm=T), dep=sum(depth,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, window_ID)]
    CDS_input$input<-CDS_input$len*CDS_input$dep
    rm("CDS_input_overlap")
    input<-CDS_input$input
    return(CDS_input)
  }))
  return(out)
}

gene_annotations_all<-fread("data/all.locus_brief_info.7.0")

gene_windows<-make_gene_windows2(data = gene_annotations_all, window=200)

H3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph", gene_windows)
H3K27me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph", gene_windows)
H3K4me1<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph", gene_windows)
H3K36me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph", gene_windows)
H2A.Z<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph", gene_windows)
H3K56ac<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph", gene_windows)
H3K4me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph", gene_windows)
Input<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph", gene_windows)

H3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph")
H3K27me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph")
H3K4me1_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph")
H3K36me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph")
H2A.Z_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph")
H3K56ac_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph")
H3K4me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph")
Input_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph")


nuclear$H3<-log2((1+H3$input)/H3_total) - log2((1+Input$input)/Input_total)
nuclear$H3K27me3<-log2((1+H3K27me3$input)/H3K27me3_total) - log2((1+Input$input)/Input_total)
nuclear$H3K4me1<-log2((1+H3K4me1$input)/H3K4me1_total) - log2((1+Input$input)/Input_total)
nuclear$H3K36me3<-log2((1+H3K36me3$input)/H3K36me3_total) - log2((1+Input$input)/Input_total)
nuclear$H2A.Z<-log2((1+H2A.Z$input)/H2A.Z_total) - log2((1+Input$input)/Input_total)
nuclear$H3K56ac<-log2((1+H3K56ac$input)/H3K56ac_total) - log2((1+Input$input)/Input_total)
nuclear$H3K4me3<-log2((1+H3K4me3$input)/H3K4me3_total) - log2((1+Input$input)/Input_total)
nuclear$H3K27me3<-log2((1+H3K27me3$input)/H3K27me3_total) - log2((1+Input$input)/Input_total)

nuclear<-gene_windows[chr %in% paste0("Chr",1:12)]
nuclear$H3K4me1_dep<-H3K4me1$input
nuclear$H3K36me3_dep<-H3K36me3$input
nuclear$H3_dep<-H3$input
nuclear$H3K27me3_dep<-H3K27me3$input
nuclear$H2A.Z_dep<-H2A.Z$input
nuclear$H3K56ac_dep<-H3K56ac$input
nuclear$H3K4me3_dep<-H3K4me3$input
nuclear$Input_dep<-Input$input



nuclear$mut<-add_vars_to_gene_windows(nuclear, ns_s[])


model<-glm(mut~H3K4me1_dep + H3K36me3_dep + H3K27me3_dep + H3K56ac_dep + H3K4me3_dep+H3_dep+H2A.Z_dep+Input_dep, nuclear, family="binomial")
summary(model)

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
predictors <- c("H3K4me1_dep", "H3K36me3_dep", "H3K27me3_dep", "H3K56ac_dep", "H3K4me3_dep", "H2A.Z_dep")

# Extract z-scores and p-values from each univariate model and add to coefs
for (predictor in predictors) {
  formula <- as.formula(paste("mut ~", predictor, "+ H3_dep + Input_dep"))
  model_uni <- glm(formula, data = nuclear, family = "binomial")
  coefs_uni <- as.data.frame(summary(model_uni)$coefficients)
  coefs_uni$Predictor <- gsub("_dep", "", rownames(coefs_uni))
  coefs_uni <- coefs_uni[coefs_uni$Predictor == gsub("_dep", "", predictor), , drop = FALSE]
  coefs_uni$model <- "univariate"
  coefs <- rbind(coefs, coefs_uni)
}

# Create the plot
pdf("figures/Lu_2019_mutation_model_windows.pdf", width=3.5, height=1.5)
p <- ggplot(coefs, aes(x = Predictor, y = `z value`, fill = -log10(`Pr(>|z|)`))) +
  geom_bar(stat = "identity", position = "dodge", col = "black") +
  scale_fill_gradient(low = "white", high = "green4") +
  labs(x = "Predictor", y = "z-score", fill = "-log10(p)") +
  theme_classic(base_size = 6) +
  theme(legend.key.size = unit(0.25, "cm"), 
        axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  facet_wrap(~ model, ncol = 2)

print(p)
dev.off()


gff<-read.GFF("data/all.clean.gff")
gff$name=substr(gff$INFO, 4,17)
gff$model=gsub("(ID=)|:.+|;.+","", gff$INFO)
cds<-gff[TYPE=="CDS"]
cds[, ID := 1:nrow(cds)]
cds<-cds[CHROM %in% paste0("Chr",1:12)]


# List of bedGraph files and their corresponding names
marks <- list(
  H3K36me3 = "data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph",
  H2A.Z = "data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph",
  H3K56ac = "data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph",
  H3K4me3 = "data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph",
  Input = "data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph",
  H3 = "data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph",
  H3K27me3 = "data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph",
  H3K4me1 = "data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph"
)

# Read in the bedGraph files and compute the values
for (mark in names(marks)) {
  message(mark)
  bedGraph_data <- read.bedGraph(marks[[mark]])
  mark_data <- features_in_features(cds, bedGraph_data, mode="sumxlength", value = "DEPTH")
  cds[[mark]]<-mark_data$sumxlength[match(cds$ID, mark_data$ID)]
}

mutations <- ns_s[MutationType == "synonymous"][, .(CHROM = Chromosome, POS = Position)]
counts<-sites_in_features(cds, mutations, mode="counts")
s$mutations<-counts$counts[match(cds$ID, counts$ID)]
cds$LENGTH<-cds$STOP-cds$START

cds_unique <- unique(cds[, -c("INFO","ID","model"), with = FALSE])

cds_summary <- cds_unique[, .(
  mutation_rate = sum(mutations) / sum(LENGTH),
  H3K36me3 = sum(H3K36me3) / sum(LENGTH),
  H2A.Z = sum(H2A.Z) / sum(LENGTH),
  H3K56ac = sum(H3K56ac) / sum(LENGTH),
  H3K4me3 = sum(H3K4me3) / sum(LENGTH),
  Input = sum(Input) / sum(LENGTH),
  H3 = sum(H3) / sum(LENGTH),
  H3K27me3 = sum(H3K27me3) / sum(LENGTH),
  H3K4me1 = sum(H3K4me1) / sum(LENGTH),
  LENGTH=sum(LENGTH),
  mutations=sum(mutations)
), by = .(gene = name)]

model<-lm(mutation_rate~H3K4me1 + H3K36me3 + H3K27me3 + H3K56ac + H3K4me3+H3+H2A.Z+Input, cds_summary)
summary(model)

# Extract coefficients, z-values, and p-values
coefs <- summary(model)$coefficients
coefs <- as.data.frame(coefs)
coefs$Predictor <- rownames(coefs)

# Filter out unwanted predictors and remove "_dep" from names
coefs <- coefs[!coefs$Predictor %in% c("(Intercept)", "H3", "Input"), ]
coefs$Predictor <- gsub("_dep", "", coefs$Predictor)

# Rank by z-value
coefs <- coefs[order(coefs$`t value`), ]
coefs$Predictor <- factor(coefs$Predictor, levels = coefs$Predictor)
coefs$model <- "multivariate"

# List of predictors for univariate models
predictors <- c("H3K4me1", "H3K36me3", "H3K27me3", "H3K56ac", "H3K4me3", "H2A.Z")

# Extract z-scores and p-values from each univariate model and add to coefs
for (predictor in predictors) {
  formula <- as.formula(paste("mutation_rate ~", predictor, "+ H3 + Input"))
  model_uni <- lm(formula, data = cds_summary)
  coefs_uni <- as.data.frame(summary(model_uni)$coefficients)
  coefs_uni$Predictor <- gsub("_dep", "", rownames(coefs_uni))
  coefs_uni <- coefs_uni[coefs_uni$Predictor == gsub("_dep", "", predictor), , drop = FALSE]
  coefs_uni$model <- "univariate"
  coefs <- rbind(coefs, coefs_uni)
}

# Create the plot
pdf("figures/Lu_2019_mutation_model_synonymous_CDS.pdf", width=3.5, height=1.5)
p <- ggplot(coefs, aes(x = Predictor, y = `t value`, fill = -log10(`Pr(>|t|)`))) +
  geom_bar(stat = "identity", position = "dodge", col = "black") +
  scale_fill_gradient(low = "white", high = "green4") +
  labs(x = "Predictor", y = "t-score", fill = "-log10(p)") +
  theme_classic(base_size = 6) +
  theme(legend.key.size = unit(0.25, "cm"), 
        axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  facet_wrap(~ model, ncol = 2)

print(p)
dev.off()


