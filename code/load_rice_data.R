
source("code/libraries_functions.R")
genome<-read.fasta("data/Osativa_204_softmasked.fa.gz")
load("data/homopolymer_rice.Rda")
homopolymers<-rbindlist(lapply(names(homopolymer_rice), function(x){
  tmp<-homopolymer_rice[[x]]
  data.table(tmp[,.(CHROM=x, START=start, STOP=end, LENGTH=length, BP=var)])
}))


gff<-fread("data/all.clean.gff", fill=T)
gff$name=substr(gff$V9, 4,17)
gff$model=gsub("(ID=)|:.+|;.+","", gff$V9)

#load rice data

LoF<-data.table(read.xlsx("data/Genes-Mutated-by-Loss-of-Function-Mutations-1504lines.xlsx", sheet=1, startRow = 4))
gene_annotations_all<-fread("data/all.locus_brief_info.7.0")
gene_annotations_all$length<-gene_annotations_all$stop-gene_annotations_all$start
gene_annotations_all$window_ID<-1:nrow(gene_annotations_all)
gene_annotations_all$ID<-1:nrow(gene_annotations_all)

gff_genes<-gff[model %in% gene_annotations_all$model]
cds<-gff_genes[V3=="CDS"]

gene_annotations_all$CDS<-sapply(gene_annotations_all$model, function(x){
  proteins<-cds[model == x]
  if(nrow(proteins)>0){
    cds_bp<-(unique(unlist(apply(unique(proteins[,-9]), 1, function(x){
      as.numeric(x["V4"]):as.numeric(x["V5"])
    }))))
    return(length(cds_bp))
  } else return(NA)
})

gene_annotations_all$length<-gene_annotations_all$stop-gene_annotations_all$start
gene_annotations_all$pct_CDS<-gene_annotations_all$CDS/gene_annotations_all$length
setkey(gene_annotations_all, chr, start, stop)
TE<-gene_annotations_all[is_TE=="Y" & is_representative=="Y"]
gene_annotations_basic<-gene_annotations_all[chr %in% paste0("Chr",1:12) & is_TE=="N" & is_representative=="Y"]
gene_annotations_basic$direction<-gene_annotations_basic$ori
gene_annotations_basic$type="gene"

snps<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=1, startRow = 4))
snps$CHROM<-snps$Chromosome
snps$POS<-snps$Position
snps$unique<-paste0(snps$CHROM, snps$POS, snps$Single.base.substitution)
snps$chr<-snps$CHROM
snps$start<-snps$POS
snps$stop<-snps$POS
snps<-snps[!is.na(start)]
snps$REF<-substr(snps$Single.base.substitution,1,1)
snps$ALT<-substr(snps$Single.base.substitution,3,3)
snps$context<-contexts(snps, genome)
context_table<-data.table(table(context=snps$context))
context_table$context_only<-substr(context_table$context, 1, 3)
context_table$mut<-paste(substr(context_table$context, 2,2),substr(context_table$context, 6,6), sep=">")
snps$context10<-long_context(snps, 10)
snps$minus10<-substr(snps$context10, 1, 10)
snps$plus10<-substr(snps$context10, 12, 21)
snps$ID<-1:nrow(snps)
snps_homopolymer_neighbor<-polymorphology2::homopolymer_var_annotate(snps, homopolymers, size = 3, dist = 1)
snps$HP<-snps_homopolymer_neighbor$HP[match(snps$ID, snps_homopolymer_neighbor$ID)]
table(snps$HP)

count<-data.table(table(u=snps$unique))
snps$N<-count$N[match(snps$unique, count$u)]
table(snps[N<2]$HP)


snps[, c("CHROM", "START", "STOP", "ID") := .(CHROM, POS, POS, ID)]

# Calculating mappability for each snp
mapp<-fread("data/rice_genome/genmap/Osativa_204_softmasked.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp<-mapp[CHROM %in% paste0("Chr",1:12)]
snps_mapp <- features_in_features(snps, mapp, mode = "mean", value = "MAPPABILITY")
snps[, MAPPABILITY := snps_mapp$mean[match(ID, snps_mapp$ID)]]


setkey(snps, chr, start, stop)


Athal<-fread("data/A_thal_germ_SBS.txt")
Athal_mut<-fread("data/A_thaliana_germline_mutations.txt")

del<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=2, startRow = 2))
del$CHROM<-del$Chromosome
del$POS<-del$Start.position
ins<-data.table(read.xlsx("data/Mutations-Identified-1504lines.xlsx", sheet=3, startRow = 2))
ins$CHROM<-ins$Chromosome
ins$POS<-ins$Position
indels<-rbind(ins[`Size.(bp)`<10], del[`Size.(bp)`<10], fill=T)
indels$unique<-paste0(indels$CHROM, indels$POS)
indels$start<-indels$POS
indels$stop<-indels$POS
indels$chr<-indels$CHROM
indels$unique<-paste0(indels$CHROM, indels$POS)
indels<-indels[!is.na(start)]
indels$ID<-1:nrow(indels)
setkey(indels, chr, start, stop)

#http://glab.hzau.edu.cn/RiceENCODE/download/Peaks/Nipponbare/Histone/

H3K4me1<-read_encode("data/H3K4me1_seedlings.bed")
H3K4me1_panicles<-read_encode("data/H3K4me1_panicles.bed")
H3K4me1_leaves<-read_encode("data/H3K4me1_mature_leaves.bed")
H3K9me1<-read_encode("data/SRR3213598_H3K9me1_seedlings.bed")
H3K36me3<-read_encode("data/SRR094791_H3K36me3_seedlings.bed")
H3K4me3<-read_encode("data/H3K4me3_seedlings.bed")
H3K9me2<-read_encode("data/H3K9me2_seedlings.bed")
H3K27ac<-read_encode("data/H3K27ac_seedlings.bed")
PII<-read_encode("data/PII_seedlings.bed")
H3K27me3<-read_encode("data/H3K27me3_seedlings.bed")
H3K4ac<-read_encode("data/SRR3213599_H3K4ac_seedlings.bed")
H3K12ac<-read_encode("data/SRR6510886_H4K12ac_seedlings.bed")
H3K9ac<-read_encode("data/SRR6795643_H3K9ac_seedlings.bed")

allmods<-rbindlist(list(H3K4me1,
H3K9me1,
H3K36me3,
H3K4me3,
H3K9me2,
H3K27ac,
PII,
H3K27me3,
H3K4ac,
H3K12ac,
H3K9ac))

allmodssum<-allmods[,.(N=.N, length=mean(length)), by=V9]

# 3000 rice genomes

missense<-fread("data/rice_missense.txt");colnames(missense)<-c("CHROM","POS")
missense$CHROM<-paste0('Chr',missense$CHROM)

synonymous<-fread("data/rice_synonymous.txt");colnames(synonymous)<-c("CHROM","POS")
synonymous$CHROM<-paste0('Chr',synonymous$CHROM)

stop_gained<-fread("data/rice_stop_gained.txt");colnames(stop_gained)<-c("CHROM","POS")
stop_gained$CHROM<-paste0('Chr',stop_gained$CHROM)



# ns_s --------------------------------------------------------------------


ns_s<-fread("data/Mutations_Final_StopCodonAnnotated_IDAdded.csv")
ns_s$CHROM<-ns_s$Chromosome
ns_s$POS<-ns_s$Position
ns_s$unique<-paste0(ns_s$CHROM, ns_s$POS, ns_s$`Single base substitution`)
ns_s$START<-ns_s$POS
ns_s$STOP<-ns_s$POS
ns_s$REF<-ns_s$Old
ns_s$ALT<-ns_s$New

count<-data.table(table(u=ns_s$unique))
ns_s$N<-count$N[match(ns_s$unique, count$u)]
ns_s$chr<-ns_s$CHROM
ns_s$start<-ns_s$POS
ns_s$Mutant.ID<-ns_s$`Mutant ID`
ns_s$stop<-ns_s$POS
ns_s$gene<-gsub("\\..+", "", ns_s$ID)
ns_s$ID<-1:nrow(ns_s)
ns_s_homopolymer_neighbor<-polymorphology2::homopolymer_var_annotate(ns_s, homopolymers, size = 3, dist = 1)
ns_s$HP<-snps_homopolymer_neighbor$HP[match(ns_s$ID, snps_homopolymer_neighbor$ID)]
ns<-ns_s[MutationType=="non-synonymous"]
s<-ns_s[MutationType=="synonymous"]
ns_s$context<-contexts(ns_s, genome)

setkey(ns_s, chr, start, stop)
setkey(s, chr, start, stop)
setkey(ns, chr, start, stop)


# Calculating mappability for each gene
mapp<-fread("data/rice_genome/genmap/Osativa_204_softmasked.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp<-mapp[CHROM %in% paste0("Chr",1:12)]
ns_s_mapp <- features_in_features(ns_s, mapp, mode = "mean", value = "MAPPABILITY")
ns_s[, MAPPABILITY := ns_s_mapp$mean[match(ID, ns_s_mapp$ID)]]
ns_s$mut<-paste0(substr(ns_s$context, 2,2), ">", substr(ns_s$context, 6,6))
table(ns_s$mut)

ns_s_filt<-ns_s[N<2 & MAPPABILITY==1]
fwrite(ns_s_filt[,.(CHROM, POS, REF, ALT, LINE=Mutant.ID, MAPPABILITY, CONTEXT=context, HP, CDS=WithinCDS, EFFECT=MutationType)], "tables/S3_rice_mutations.csv")


table(ns_s_filt$HP)

pdf("figures/kitaake_SBS.pdf", width=6, height=1.5)
plot_tricontexts(ns_s_filt$context)
dev.off()

# genes filtered ----------------------------------------------------------

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


# Renaming columns for consistency
gene_annotations_all[, c("CHROM", "START", "STOP", "ID") := .(chr, start, stop, .I)]
gene_annotations_all<-gene_annotations_all[CHROM %in% paste0("Chr",1:12)]
# Calculating mappability for each gene
gene_annotations_all_mapp <- features_in_features(gene_annotations_all, mapp, mode = "mean", value = "MAPPABILITY")
gene_annotations_all[, MAPPABILITY := gene_annotations_all_mapp$mean[match(ID, gene_annotations_all_mapp$ID)]]

# Processing bedGraph data for specific marks and adding to gene annotations
for (mark in c("H3K4me1", "H3")) {
  message("Processing mark: ", mark)
  bedGraph_data <- read.bedGraph(marks[[mark]])
  mark_data <- features_in_features(gene_annotations_all, bedGraph_data, mode = "sumxlength", value = "DEPTH")
  gene_annotations_all[[mark]] <- mark_data$sumxlength[match(gene_annotations_all$ID, mark_data$ID)]
}

# Adding mutation counts to gene annotations
mutations <- ns_s_filt[N < 2, .(CHROM = Chromosome, POS = Position)]
counts <- sites_in_features(gene_annotations_all, mutations, mode = "counts")
gene_annotations_all[, mutations := counts$counts[match(ID, counts$ID)]]

# Adding specific mutation types counts
gene_annotations_all[, syn_mut := sites_in_features(gene_annotations_all, ns_s_filt[ MutationType == "synonymous"], mode = "counts")$counts]
gene_annotations_all[, ns_mut := sites_in_features(gene_annotations_all, ns_s_filt[ MutationType == "non-synonymous"], mode = "counts")$counts]
gene_annotations_all[, stop_mt := sites_in_features(gene_annotations_all, ns_s_filt[ `StopCodon?` != ""], mode = "counts")$counts]

# Adding missense and synonymous mutation counts

gene_annotations_all[, `:=`(Pn = sites_in_features(gene_annotations_all, missense, mode = "counts")$counts,
                            Ps = sites_in_features(gene_annotations_all, synonymous, mode = "counts")$counts,
                            stop_gained = sites_in_features(gene_annotations_all, stop_gained, mode = "counts")$counts)]

# Calculating pnps and gene length
gene_annotations_all[, `:=`(pnps = (Pn+stop_gained) / Ps, LENGTH = STOP - START)]

# Calculating enrichment
gene_annotations_all[, enrich := log(H3K4me1 / H3)]

# Processing GFF data
gff <- read.GFF("data/all.clean.gff")
gff[, `:=`(name = substr(INFO, 4, 17), model = gsub("(ID=)|:.+|;.+","", INFO))]
cds <- gff[TYPE == "CDS", .(CHROM, START, STOP, TYPE, model, ID = .I)]

# Filtering CDS by chromosome and calculating length
cds <- cds[CHROM %in% paste0("Chr", 1:12)]
cds[, LENGTH := STOP - START + 1]

# Adding mutation counts to CDS
cds_mutation <- sites_in_features(cds, ns_s_filt[N < 2, .(CHROM = Chromosome, POS = Position)], mode = "counts")
cds[, mutation := cds_mutation$counts[match(ID, cds_mutation$ID)]]

# Summarizing mutations and lengths by model
cds_sum <- cds[, .(CDS_mutation = sum(mutation), CDS_length = sum(LENGTH)), by = model]
gene_annotations_all[, `:=`(CDS_mutation = cds_sum$CDS_mutation[match(model, cds_sum$model)],
                            CDS = cds_sum$CDS_length[match(model, cds_sum$model)])]

# Processing introns and UTRs
models <- gff[TYPE %in% c("mRNA", "exon"), .(CHROM, START, STOP, TYPE, model)]
introns <- models[TYPE != "mRNA", find_introns(.SD), by = model]
introns[, CHROM := gene_annotations_all$CHROM[match(model, gene_annotations_all$model)]]

UTR <- gff[grepl("UTR", TYPE), .(CHROM, START, STOP, TYPE, model)]
nCDS <- rbind(introns[TYPE == "intron"], UTR)
nCDS[, `:=`(LENGTH = STOP - START + 1, ID = .I)]

# Adding mutation counts to nCDS
nCDS_mutation <- sites_in_features(nCDS, ns_s_filt[N < 2, .(CHROM = Chromosome, POS = Position)], mode = "counts")
nCDS[, mutation := nCDS_mutation$counts[match(ID, nCDS_mutation$ID)]]

# Summarizing mutations and lengths by model for nCDS
nCDS_sum <- nCDS[, .(nCDS_mutation = sum(mutation), nCDS_length = sum(LENGTH)), by = model]
gene_annotations_all[, `:=`(nCDS_mutation = nCDS_sum$nCDS_mutation[match(model, nCDS_sum$model)],
                            nCDS = nCDS_sum$nCDS_length[match(model, nCDS_sum$model)])]
gene_annotations_all[, `:=`(nCDS_mutation = fifelse(is.na(nCDS_mutation), 0, nCDS_mutation),
                            nCDS = fifelse(is.na(nCDS), 0, nCDS))]

# Filtering genes
genes_filt <- gene_annotations_all[is.finite(enrich) & is.finite(pnps) & MAPPABILITY == 1 & is_expressed == "Y"]

# gene windows ------------------------------------------------------------

gene_annotations_all$DIRECTION<-gene_annotations_all$ori
gene_windows<-polymorphology2::feature_windows(gene_annotations_all, breaks = 10, dist=5000, directed = T, IDcol = "model")

for (mark in c("H3K4me1","H3")) {
  message(mark)
  bedGraph_data <- read.bedGraph(marks[[mark]])
  mark_data <- features_in_features(gene_windows, bedGraph_data, mode="sumxlength", value = "DEPTH")
  gene_windows[[mark]]<-mark_data$sumxlength[match(gene_windows$ID, mark_data$ID)]
}


gene_windows_mapp<-features_in_features(gene_windows,mapp, mode = "mean", value = "MAPPABILITY" )

gene_windows$MAPPABILITY<-gene_windows_mapp$mean[match(gene_windows_mapp$ID, gene_windows_mapp$ID)]

gene_windows$enrich<-log(gene_windows$H3K4me1/gene_windows$H3)
gene_windows_muts<-sites_in_features(gene_windows, ns_s_filt, mode="counts")
gene_windows$mutations<-gene_windows_muts$counts[match(gene_windows$ID, gene_windows_muts$ID)]



# windows_for_regression_analysis -----------------------------------------



windowsize=500
windows500<-rbindlist(apply(gene_annotations_all, 1, function(r){
  message(r["model"])
  CHROM=r["CHROM"]
  START=as.numeric(r["START"])-5000
  STOP=as.numeric(r["STOP"])+5000
  starts=seq(from=START, to=STOP-windowsize, by=windowsize)
  stops=seq(from=START+windowsize, to=STOP, by=windowsize)
  return(data.table(CHROM, START=starts, STOP=stops))
}))

#windows500<-fread("data/rice_500bp_windows_epigenome_lu2019.csv")
windows500$ID<-1:nrow(windows500)


for (mark in names(marks)) {
  message(mark)
  bedGraph_data <- read.bedGraph(marks[[mark]])
  mark_data <- features_in_features(windows500, bedGraph_data, mode="sumxlength", value = "DEPTH")
  windows500[[paste0(mark,"_dep")]]<-mark_data$sumxlength[match(windows500$ID, mark_data$ID)]
}

# H3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph", windows500)
# H3K27me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph", windows500)
# H3K4me1<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph", windows500)
# H3K36me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph", windows500)
# H2A.Z<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph", windows500)
# H3K56ac<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph", windows500)
# H3K4me3<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph", windows500)
# Input<-chip_overlaps_Lu("data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph", windows500)

H3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674680_ChIP_Rice_7days_leaf_H3.bedGraph")
H3K27me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674681_ChIP_Rice_7days_leaf_H3K27me3.bedGraph")
H3K4me1_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph")
H3K36me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674682_ChIP_Rice_7days_leaf_H3K36me3.bedGraph")
H2A.Z_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674686_ChIP_Rice_7days_leaf_H2A.Z.bedGraph")
H3K56ac_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674683_ChIP_Rice_7days_leaf_H3K56ac.bedGraph")
H3K4me3_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674684_ChIP_Rice_7days_leaf_H3K4me3.bedGraph")
Input_total<-chip_total("data/Lu2019/GSE128434_RAW/GSM3674687_ChIP_Rice_7days_leaf_Input.bedGraph")

# windows500$H3K4me1_dep<-H3K4me1$input
# windows500$H3K36me3_dep<-H3K36me3$input
# windows500$H3_dep<-H3$input
# windows500$H3K27me3_dep<-H3K27me3$input
# windows500$H2A.Z_dep<-H2A.Z$input
# windows500$H3K56ac_dep<-H3K56ac$input
# windows500$H3K4me3_dep<-H3K4me3$input
# windows500$Input_dep<-Input$input

windows500$H3<-log2((1+windows500$H3_dep)/H3_total) - log2((1+windows500$Input_dep)/Input_total)
windows500$H3K27me3<-log2((1+windows500$H3K27me3_dep)/H3K27me3_total) - log2((1+windows500$H3_dep)/H3_total)
windows500$H3K4me1<-log2((1+windows500$H3K4me1_dep)/H3K4me1_total) - log2((1+windows500$H3_dep)/H3_total)
windows500$H3K36me3<-log2((1+windows500$H3K36me3_dep)/H3K36me3_total) - log2((1+windows500$H3_dep)/H3_total)
windows500$H2A.Z<-log2((1+windows500$H2A.Z_dep)/H2A.Z_total) - log2((1+windows500$H3_dep)/H3_total)
windows500$H3K56ac<-log2((1+windows500$H3K56ac_dep)/H3K56ac_total) - log2((1+windows500$H3_dep)/H3_total)
windows500$H3K4me3<-log2((1+windows500$H3K4me3_dep)/H3K4me3_total) - log2((1+windows500$H3_dep)/H3_total)

windows500_muts<-sites_in_features(windows500,  ns_s_filt, mode="counts")
windows500$mut_counts<-windows500_muts$counts[match(windows500$ID, windows500_muts$ID)]
windows500$mutations<-windows500$mut_counts
windows500$mut<-windows500$mutations>0

windows500_mapp<-features_in_features(windows500,mapp, mode = "mean", value = "MAPPABILITY" )
windows500$MAPPABILITY<-windows500_mapp$mean[match(windows500$ID, windows500_mapp$ID)]

windows500_filt<-windows500[!is.na(MAPPABILITY) & MAPPABILITY==1]
fwrite(windows500_filt, "data/rice_500bp_windows_epigenome_lu2019.csv")


