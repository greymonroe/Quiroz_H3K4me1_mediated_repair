library(seqinr)
library(openxlsx)
library(data.table)
library(polymorphology) #devtools::install_github("greymonroe/polymorphology)
library(rtracklayer) #BiocManager::install("rtracklayer")
library(MASS)
library(boot)
library(polymorphology2) #devtools::install_github("greymonroe/polymorphology2)



nonsyn_syn_exp3 <- function(muts, genome, CDS, reps) {
  
  bp_freq<-data.table(table(REF=genome[[1]]))[!REF %in% c("N","n")]
  bp_freq$REF<-toupper(bp_freq$REF)
  
  # Create an empty data frame to store the mutations
  mutations <- data.frame(REF = character(0), ALT = character(0))
  
  # Loop through each nucleotide and find its possible mutations
  for (ref in nucleotides) {
    alts <- nucleotides[nucleotides != ref]
    temp_df <- data.frame(REF = ref, ALT = alts)
    mutations <- rbind(mutations, temp_df)
  }
  
  mutations<-data.table((table(REF=muts$REF, ALT=muts$ALT)))[REF!=ALT]
  mutations<-merge(mutations, bp_freq, by="REF")
  mutations$PROB=prop.table(mutations$N.x/mutations$N.y)
  mutations$PROB2=prop.table(mutations$N.x)
  unique_refs <- unique(mutations$REF)
  ref_probs <- sapply(unique_refs, function(ref) sum(mutations[REF == ref, PROB]))
  ref_probsdt<-data.table(REF=names(ref_probs), PROB=ref_probs)
  
  CDS_dt_all<-rbindlist(lapply(sample(1:length(CDS), reps), function(prot){
    random_seq<-toupper(CDS[[prot]])
    CDS_dt<-data.table(POS=1:length(random_seq), REF=random_seq)
    CDS_dt<-merge(CDS_dt, mutations[,.(REF, ALT, PROB)], by="REF",  allow.cartesian=TRUE)[order(POS),.(MUT=.I, POS, REF, ALT, PROB=prop.table(PROB))]
    mutation<-CDS_dt[sample(CDS_dt$MUT, round(length(random_seq)*.02), prob = CDS_dt$PROB)]
    
    # 6. Translate both sequences
    original_protein <- seqinr::translate(random_seq)
    
    # 5. Introduce the mutation
    
    mutation$effects<-apply(mutation, 1, function(r){
      mutated_sequence <- random_seq
      pos<-as.numeric(r["POS"])
      mutated_base<-r["ALT"]
      mutated_sequence[pos] <- mutated_base
      mutated_protein <- seqinr::translate(mutated_sequence)
      non_syn <- paste0(original_protein, collapse="") != paste0(mutated_protein, collapse="")
      non_syn<-ifelse(non_syn, "Non-Synonymous","Synonymous")
    })
    return(mutation)
  }))
  
}

ratio_calc<-function(test, bootN, N, s){
  ratios<-sapply(1:bootN, function(i){
    Ns=sum(test[sample(1:nrow(test),N)]$effects=="Non-Synonymous")-N*s
    S=sum(test[sample(1:nrow(test),N)]$effects=="Synonymous")
    ratio<-Ns/S
  })
}

# Function to find introns for each mRNA model
find_introns <- function(dt) {
  # Filter exons and sort by start position
  exons <- dt[TYPE == "exon", .(START = as.integer(START), STOP = as.integer(STOP))]
  exons <- exons[order(START)]
  exons[, TYPE := "exon"]  # Add TYPE column to exons
  
  # Calculate introns
  if (nrow(exons) > 1) {
    introns <- data.table(
      START = as.integer(exons$STOP[-nrow(exons)] + 1),
      STOP = as.integer(exons$START[-1] - 1),
      TYPE = "intron"
    )
    # Combine exons and introns and sort
    return(rbindlist(list(exons, introns), use.names = TRUE)[order(START)])
  } else {
    return(exons)
  }
}



silent_rate <- function(data) {
  data<-data[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), CDS_mutation=sum(CDS_mutation), meanpnps=mean(log(pnps)))]
  
  silent<-(data$mut) / (data$length )
  
  return(silent)
}

# Define the statistic function for bootstrapping
bootstrap_stat <- function(data, indices) {
  # Draw bootstrap samples
  sample_data <- data[indices, ]
  # Calculate the silent rate for the bootstrap sample
  silent <- silent_rate(sample_data)
  return(mean(silent))
}

bootstrap_mutation_silent<-function(data, indices){
  sample_data <- data.table(data[indices, ])
  
  genes_filt_99$syn_mut<-sites_in_features(genes_filt_99, sample_data[MutationType=="synonymous"], mode = "counts")$counts
  
  genes_filt_99$ns_mut<-sites_in_features(genes_filt_99, sample_data[MutationType=="non-synonymous"], mode = "counts")$counts
  
  gene_sums<-genes_filt_99[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), meanpnps=mean(log(pnps))), by=.(group)][order(group)]
  gene_sums$silent<-gene_sums$syn_muts/gene_sums$CDS+(gene_sums$mut-gene_sums$syn_muts-gene_sums$ns_mut)/(gene_sums$length-gene_sums$CDS)
  
  return(gene_sums$silent)
}


make_feature_windows<-function(encode_data, deciles=10, gene=F){

  windows<-rbindlist(apply(encode_data, 1, function(x) {
    
    chr=x["chr"]
    body_starts=round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), length.out=deciles+1)[-(deciles+1)])
    body_stops<-round(seq(as.numeric(x["start"]), as.numeric(x["stop"]), length.out=deciles+1)[-1])
    upstream_starts<-seq(as.numeric(x["start"])-2000, as.numeric(x["start"]), length.out=deciles+1)[-(deciles+1)]
    upstream_stops<-seq(as.numeric(x["start"])-2000, as.numeric(x["start"]), length.out=deciles+1)[-1]
    downstream_starts<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+2000, length.out=deciles+1)[-(deciles+1)]
    downstream_stops<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+2000, length.out=deciles+1)[-1]

    out<-data.table(chr=x["chr"],
                    start=c(upstream_starts, body_starts, downstream_starts),
                    stop=c(upstream_stops, body_stops, downstream_stops),
                    region=c(rep("upstream", length(upstream_starts)),rep("gene body", length(body_starts)),rep("downstream", length(downstream_starts))),
                    ID=as.numeric(x["ID"]))

    out$pos<-1:nrow(out)
    out$length<-out$stop-out$start
    
    if(gene == T){
      direction=x["direction"]
      if(direction=="-"){
        out$pos<-rev(out$pos)
        out$region<-rev(out$region)
      }
    }
    
    return(out)
    
  }))
  
  setkey(windows, chr, start, stop)
  return(windows)
}

read_encode<-function(file){
  encode_data<-fread(file)
  encode_data$V1<-gsub("chr","Chr", encode_data$V1)
  encode_data$chr<-gsub("chr","Chr", encode_data$V1)
  encode_data$ID<-1:nrow(encode_data)
  encode_data$start<-encode_data$V2
  encode_data$stop<-encode_data$V3
  encode_data$length<-encode_data$stop-encode_data$start+1
  
  setkey(encode_data, chr, start, stop)
  
  return(encode_data)
}

make_gene_windows<-function(data, window=150){
  deciles<-3000/window
  windows<-rbindlist(apply(data, 1, function(x) {
    
    chr=x["chr"]
    body_starts=seq(as.numeric(x["start"]), as.numeric(x["stop"]), by=3000/deciles);body_starts<-body_starts[-length(body_starts)]
    body_stops<-seq(as.numeric(x["start"]), as.numeric(x["stop"]),  by=3000/deciles)[-1]
    upstream_starts<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-(deciles+1)]
    upstream_stops<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-1]
    downstream_starts<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-(deciles+1)]
    downstream_stops<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-1]
    
    out<-data.table(chr=x["chr"], 
                    start=c(upstream_starts, body_starts, downstream_starts),
                    stop=c(upstream_stops, body_stops, downstream_stops),
                    region=c(rep("upstream", length(upstream_starts)),rep("gene body", length(body_starts)),rep("downstream", length(downstream_starts))),
                    ID=x["ID"], 
                    gene=x["locus"])
    out$pos<-1:nrow(out)
    out$length<-out$stop-out$start
    return(out)
    
  }))
  windows$window_ID<-1:nrow(windows)
  setkey(windows, chr, start, stop)
  gene_windows<-windows
  gene_windows$H3K4me1<-encode_overlap(H3K4me1, gene_windows)
  gene_windows$H3K9me1<-encode_overlap(H3K9me1, gene_windows)
  gene_windows$H3K4me3<-encode_overlap(H3K4me3, gene_windows)
  gene_windows$H3K36me3<-encode_overlap(H3K36me3, gene_windows)
  gene_windows$H3K9me2<-encode_overlap(H3K9me2, gene_windows)
  gene_windows$H3K27me3<-encode_overlap(H3K27me3, gene_windows)
  gene_windows$H3K27ac<-encode_overlap(H3K27ac, gene_windows)
  gene_windows$PII<-encode_overlap(PII, gene_windows)
  gene_windows$H3K4ac<-encode_overlap(H3K4ac, gene_windows)
  gene_windows$H3K12ac<-encode_overlap(H3K12ac, gene_windows)
  gene_windows$H3K9ac<-encode_overlap(H3K9ac, gene_windows)
  
  return(gene_windows)
}

make_genome_windows<-function(genome, window, overlap=T){
  lengths<-lapply(genome, length)
  names(lengths)<-names(genome)
  chrs<-rbindlist(lapply(names(genome)[1:12], function(c){
    starts<-seq(1, lengths[[c]], by=window);starts<-starts[-length(starts)]
    stops<-c(starts[-1],lengths[[c]])
    return(data.table(chr=c, start=starts, stop=stops))
  }))
  gene_windows<-chrs
  gene_windows$window_ID<-1:nrow(gene_windows)
  setkey(gene_windows, chr, start, stop)
  if(overlap==T){
  gene_windows$H3K4me1<-encode_overlap(H3K4me1, gene_windows)
  gene_windows$H3K9me1<-encode_overlap(H3K9me1, gene_windows)
  gene_windows$H3K4me3<-encode_overlap(H3K4me3, gene_windows)
  gene_windows$H3K36me3<-encode_overlap(H3K36me3, gene_windows)
  gene_windows$H3K9me2<-encode_overlap(H3K9me2, gene_windows)
  gene_windows$H3K27me3<-encode_overlap(H3K27me3, gene_windows)
  gene_windows$H3K27ac<-encode_overlap(H3K27ac, gene_windows)
  gene_windows$PII<-encode_overlap(PII, gene_windows)
  gene_windows$H3K4ac<-encode_overlap(H3K4ac, gene_windows)
  gene_windows$H3K12ac<-encode_overlap(H3K12ac, gene_windows)
  gene_windows$H3K9ac<-encode_overlap(H3K9ac, gene_windows)
  gene_windows$snp<-add_vars_to_gene_windows(gene_windows, snps)
  gene_windows$snp_noCT<-add_vars_to_gene_windows(gene_windows, snps[!Single.base.substitution %in% c("C>T","G>A")])
  gene_windows$snp_homoz<-add_vars_to_gene_windows(gene_windows, snps[Genotype=="Homozygous"])
  gene_windows$indel<-add_vars_to_gene_windows(gene_windows, indels)
  gene_windows$chr<-factor(gene_windows$chr, levels=paste0("Chr",1:12))
  } else {
    gene_windows$H3K4me1<-encode_hits(H3K4me1, gene_windows)
    gene_windows$H3K9me1<-encode_hits(H3K9me1, gene_windows)
    gene_windows$H3K4me3<-encode_hits(H3K4me3, gene_windows)
    gene_windows$H3K36me3<-encode_hits(H3K36me3, gene_windows)
    gene_windows$H3K9me2<-encode_hits(H3K9me2, gene_windows)
    gene_windows$H3K27me3<-encode_hits(H3K27me3, gene_windows)
    gene_windows$H3K27ac<-encode_hits(H3K27ac, gene_windows)
    gene_windows$PII<-encode_hits(PII, gene_windows)
    gene_windows$H3K4ac<-encode_hits(H3K4ac, gene_windows)
    gene_windows$H3K12ac<-encode_hits(H3K12ac, gene_windows)
    gene_windows$H3K9ac<-encode_hits(H3K9ac, gene_windows)
    gene_windows$snp<-add_vars_hits_to_gene_windows(gene_windows, snps)
    gene_windows$snp_noCT<-add_vars_hits_to_gene_windows(gene_windows, snps[!Single.base.substitution %in% c("C>T","G>A")])
    gene_windows$snp_homoz<-add_vars_hits_to_gene_windows(gene_windows, snps[Genotype=="Homozygous"])
    gene_windows$indel<-add_vars_hits_to_gene_windows(gene_windows, indels)
    gene_windows$chr<-factor(gene_windows$chr, levels=paste0("Chr",1:12))
    
  }
  setkey(gene_windows, chr, start, stop)
  
  return(gene_windows)
}

encode_overlap<-function(encode_data, gene_windows){
  encode_overlap<-foverlaps(gene_windows, encode_data,type="any")
  marked<-encode_overlap[,.(marked=ifelse(sum(!is.na(ID))==0, "unmarked","marked")), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  return(as.factor(marked$marked))
}

encode_hits<-function(encode_data, gene_windows, out="marked"){
  encode_overlap<-foverlaps(gene_windows, encode_data,type="any")
  marked<-encode_overlap[,.(marked=sum(!is.na(ID)), length=sum(length,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  if(out=="length"){
  return(as.numeric(marked$length))
  } else  return(as.numeric(marked$marked))
}

plot_peaks<-function(encode_data, ggtitle, xtitle, deciles, var_data, gene=F, marky=F){
  
  windows<-make_feature_windows(encode_data =  encode_data, deciles, gene=gene)
  windows_vars<-foverlaps(windows, var_data,type="any")
  #windows_vars[ID==1 & region=="upstream"]
  windows_vars[ID==1 & region=="gene body"]
  windows_means<-windows_vars[,.(mut=sum(!is.na(start)), N=.N, length=mean(length, na.rm=T)),by=.(pos, region, i.ID)]
  windows_means<-windows_means[,.(pct=sum(mut)/sum(length), mut=sum(mut), length=sum(length)), by=.(pos, region)]
  if(marky==T){windows_means$pct<-windows_means$mut}
  plot<-ggplot(windows_means, aes(x=pos, y=pct, col=region=="gene body", group=region))+
    geom_vline(xintercept = c(deciles, deciles*2), linetype="dashed", size=0.25)+
    geom_line()+
    scale_color_manual(values=c("gray75","green3"), guide="none")+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="Mutations/bp")+
    ggtitle(ggtitle)+
    scale_x_continuous(breaks=c(0,max(windows_means$pos)/3, max(windows_means$pos)/3*2, max(windows_means$pos)), labels=c("-2kb","0%","100%","+2kb"), name=xtitle)
return(list(plot, windows_means))
  }

add_vars_to_gene_windows<-function(gene_windows, var_object){
  vars_overlap<-foverlaps(gene_windows, var_object,type="any")
  vars<-vars_overlap[,.(mutations=sum(!is.na(Mutant.ID))), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  vars$mutated<-vars$mutations>0
  return(vars$mutated)
}

add_vars_hits_to_gene_windows<-function(gene_windows, var_object){
  vars_overlap<-foverlaps(gene_windows, var_object,type="any")
  vars<-vars_overlap[,.(mutations=sum(!is.na(Mutant.ID))), by=.(chr, start=i.start, stop=i.stop, window_ID)]
  vars$mutated<-vars$mutations>0
  return(vars$mutations)
}


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

plot_model<-function(model_sum, ggtitle){
  ggplot(model_sum[predictor!="(Intercept)"], aes(x=predictor, y=y, fill=log10(`P`)))+
    geom_bar(stat="identity", col="black",size=0.25)+
    theme_classic(base_size = 6)+
    scale_fill_gradientn(colors=c("dodgerblue","white"), name="-log10(P)")+
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.key.size = unit(.3,"line"), legend.key=element_rect(color="black"))+
    scale_x_discrete(name="Predictor")+
    ggtitle(ggtitle)
}

log_model<-function(gene_windows, variable, aic=F){
  form<-formula(paste0("as.numeric(",variable,")~H3K4me1+H3K9me1+H3K4me3+H3K27me3+H3K9me2+H3K27ac+H3K36me3+PII+H3K4ac+H3K12ac+H3K9ac"))
  model<-glm(form, gene_windows, family="binomial")
  if(aic==T){
    model<-stepAIC(model, direction="both", trace = F)
  }
  model_sum<-summary(model)
  model_sum<-data.table(model_sum$coefficients)
  model_sum$predictor<-gsub("unmarked","",row.names(summary(model)$coefficients))
  model_sum$Estimate<--model_sum$Estimate
  model_sum$P<-model_sum$`Pr(>|z|)`
  model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(-model_sum$`z value`)])
  model_sum$y<--model_sum$`z value`
  return(model_sum)
}

log_model_single<-function(gene_windows, variable){
  
  preds<-c("H3K4me1","H3K9me1","H3K4me3","H3K27me3","H3K9me2","H3K27ac","H3K36me3","PII","H3K4ac","H3K12ac","H3K9ac")
  models<-lapply(preds, function(pred){
    cat(pred)
    form<-formula(paste0("as.numeric(",variable,")~",pred))
    model<-summary(glm(form, gene_windows, family="binomial"))
    model_sum<-data.table(model$coefficients)
    model_sum$predictor<-gsub("unmarked","",row.names(model$coefficients))
    model_sum$Estimate<--model_sum$Estimate
    model_sum$P<-model_sum$`Pr(>|z|)`
    model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(-model_sum$`z value`)])
    model_sum$y<--model_sum$`z value`
    return(model_sum)
  })
  models<-rbindlist(models)[predictor!="(Intercept)"]
  models$predictor<-factor(models$predictor, levels=models$predictor[order(-models$`z value`)])
  return(models)
}

lm_model<-function(gene_windows, variable, aic=F){
  form<-formula(paste0("as.numeric(",variable,")~H3K4me1+H3K4me3+H3K9me1+H3K27me3+H3K9me2+H3K27ac+H3K36me3+PII+H3K4ac+H3K12ac+H3K9ac"))
  model<-lm(form, gene_windows)
  if(aic==T){model<-MASS::stepAIC(model, direction="both", trace = F)}
  model_sum<-data.table(summary(model)$coefficients)
  model_sum$predictor<-gsub("unmarked","",row.names(summary(model)$coefficients))
  model_sum$Estimate<-model_sum$Estimate
  model_sum$P<-model_sum$`Pr(>|t|)`
  model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(model_sum$`t value`)])
  model_sum$y<-model_sum$`t value`
  return(model_sum)
}

poisson_model<-function(gene_windows, variable, aic=F){
  form<-formula(paste0("as.numeric(",variable,")~H3K4me1+H3K4me3+H3K9me1+H3K27me3+H3K9me2+H3K27ac+H3K36me3+PII+H3K4ac+H3K12ac+H3K9ac"))
  model<-glm(form, gene_windows, family="poisson")
  if(aic==T){model<-MASS::stepAIC(model, direction="both", trace = F)}
  model_sum<-data.table(summary(model)$coefficients)
  model_sum$predictor<-gsub("unmarked","",row.names(summary(model)$coefficients))
  model_sum$Estimate<-model_sum$Estimate
  model_sum$P<-model_sum$`Pr(>|z|)`
  model_sum$predictor<-factor(model_sum$predictor, levels=model_sum$predictor[order(model_sum$`z value`)])
  model_sum$y<-model_sum$`z value`
  return(model_sum)
}

chip_total<-function(bedfile){
  in1<-fread(bedfile)
  total<-sum(in1$V4*(in1$V3-in1$V2))
}

chip_overlaps<-function(bedfile, featureobject){
  cat("\nreading ");cat(bedfile)
  in1<-fread(bedfile)
  colnames(in1)<-c("chr","start","stop","depth")
  in1$length<-as.numeric(in1$stop-in1$start)
  in1$depth<-as.numeric(in1$depth)
  setkey(in1, chr, start, stop)
  out<-unlist(lapply(1:5, function(c) {
    cat(" chr ");cat(c)
    CDS_input_overlap<-foverlaps(featureobject[chr==c], in1[chr==c],type="any")
    CDS_input<-CDS_input_overlap[,.(len=sum(length,na.rm=T), dep=sum(depth,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, gene)]
    CDS_input$input<-CDS_input$len*CDS_input$dep
    rm("CDS_input_overlap")
    input<-CDS_input$input
    return(input)
  }))
  return(out)
}

chip_overlaps_window<-function(bedfile, featureobject){
  cat("\nreading ");cat(bedfile)
  in1<-fread(bedfile)
  featureobject$ID<-1:nrow(featureobject)
  colnames(in1)<-c("chr","start","stop","depth")
  in1$length<-as.numeric(in1$stop-in1$start)
  in1$depth<-as.numeric(in1$depth)
  setkey(in1, chr, start, stop)
  out<-unlist(lapply(1:5, function(c) {
    cat(" chr ");cat(c)
    CDS_input_overlap<-foverlaps(featureobject[chr==c], in1[chr==c],type="any")
    CDS_input<-CDS_input_overlap[,.(len=sum(length,na.rm=T), dep=sum(depth,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, ID)]
    CDS_input$input<-CDS_input$len*CDS_input$dep
    rm("CDS_input_overlap")
    input<-CDS_input$input
    return(input)
  }))
  return(out)
}

peaks_randomized<-function(featureobject){
  rand<-rbindlist(apply(rbindlist(list(featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject,featureobject)), 1, function(x){
    chr<-x["chr"]
    start<-sample(1:length(genome[[chr]]), 1)
    stop<-start+as.numeric(x["length"])
    return(data.table(chr=chr, start=start, stop=stop))
  }))
  rand$length<-rand$stop-rand$start
  rand$ID<-1:nrow(rand)
  setkey(rand, chr, start, stop)
}

plot_bars_rice<-function(sumstable, yvar, xlab, ylab, ggtitle){
  ggplot(gene_annotations_all_sums, aes_string(x="grp", y=yvar))+
    geom_bar(stat="identity", fill="dodgerblue4", col="black", width=0.5)+
    theme_classic(base_size = 6)+
    scale_x_discrete(name=xlab)+
    scale_y_continuous(name=ylab)+
    ggtitle(ggtitle)
  
}


bw_read<-function(bwfile){
  bw_dt<-data.table(data.frame(import(con=bwfile)))
  colnames(bw_dt)[1]<-"chr"
  bw_dt$chr<-as.numeric(gsub("chr","",bw_dt$chr))
  colnames(bw_dt)[3]<-"stop"
  colnames(bw_dt)[4]<-"length"
  colnames(bw_dt)[6]<-"depth"
  setkey(bw_dt, chr, start, stop)
}

bw_overlaps<-function(bw_object, window_object){
  out<-unlist(lapply(unique(window_object$chr), function(c) {
    cat(" chr ");cat(c)
    window_input_overlap<-foverlaps(window_object[chr==c], bw_object[chr==c],type="any")
    window_input<-window_input_overlap[,.(len=sum(length,na.rm=T), dep=sum(depth,na.rm=T)), by=.(chr, start=i.start, stop=i.stop, window_id)]
    window_input$input<-window_input$len*window_input$dep
    rm("window_input_overlap")
    input<-window_input$input
    return(input)
  }))
  return(out)
}

convert_to_format <- function(num) {
  # Convert the number to scientific notation
  sci_num <- format(num, scientific = TRUE, digits = 1)
  
  # Separate the coefficient and exponent parts
  parts <- strsplit(sci_num, split = "e")[[1]]
  
  # Format the string as needed using bquote
  formatted_num <- bquote("p=" ~ .(as.numeric(parts[1])) %*% 10^.(as.numeric(parts[2])))
  formatted_num <- bquote(italic(p) == .(as.numeric(parts[1])) %*% 10^.(as.numeric(parts[2])))
  
  return(formatted_num)
}

line_data<-function(mutations){
  enrichment_line<-mutations[,
                             .(H3K4me1_mean=median(enrich_H3K4me1), 
                               H3K4me1=sum(enrich_H3K4me1>1)/sum(enrich_H3K4me1<1), 
                               H3K4me1_pct=sum(enrich_H3K4me1>1)/.N*100,
                               HP_pct=sum(HP)/.N*100,
                               H3K4me1_1=sum(enrich1_H3K4me1>0)/.N, 
                               H3K4me1_2=sum(enrich2_H3K4me1>0)/.N,
                             all_tissue=sum(tissue_all)/sum(genic)*100,
                               #genic=sum(genic)/sum(!genic), 
                               genic_pct=sum(genic)/.N*100,
                               ox=sum(mut %in% c("C>A","T>G"))/sum(!mut %in% c("C>A","T>G")),
                               N=.N), 
                             by=.(file, trt)]
  return(enrichment_line)
}


line_barplot<-function(enrichment_line, response, yaxis, title){
  
  tmp<-enrichment_line[,c("trt",response), with=F]
  colnames(tmp)[2]<-"response"
  tmp<-tmp[is.finite(response)]
  tmp$trt<-factor(tmp$trt, levels=c("WT","MSH6"))
  
  
  test<-t.test(tmp$response~tmp$trt)
  p<-test$p.value
  means<-tmp[,.(mean=mean(response), se=sd(response)/sqrt(.N)),by=.(trt)]
  plot<-ggplot(means, aes(x=trt, y=mean,  fill=trt))+
    theme_classic(base_size = 6)+
    geom_bar(stat="identity", col="black", width=0.75)+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0)+
    scale_fill_manual(values = c("dodgerblue","orange2"), guide="none")+
    scale_y_continuous(name=yaxis)+
    scale_x_discrete(name="Genotype", labels=c("WT","msh6"))+
    theme(plot.title = element_text(hjust = 0.5, size=6))+
    ggtitle(label=title, subtitle = convert_to_format(p))
  plot
}

line_boxplot<-function(enrichment_line, response, yaxis, title){
  
  tmp<-enrichment_line[,c("trt",response), with=F]
  colnames(tmp)[2]<-"response"
  tmp$trt<-factor(tmp$trt, levels=c("WT","MSH6"))
  
  
  test<-t.test(tmp$response~tmp$trt)
  p<-test$p.value
  plot<-ggplot(tmp, aes(x=trt, y=response, col=trt))+
    theme_classic(base_size = 6)+
    geom_boxplot(fill=NA, outlier.size = 0.5)+
    scale_color_manual(values = c("dodgerblue","orange2"), guide="none")+
    scale_y_continuous(name=yaxis)+
    scale_x_discrete(name="Genotype", labels=c("WT","msh6"))+
    theme(plot.title = element_text(hjust = 0.5, size=6), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))+
    ggtitle(label=title, subtitle = convert_to_format(p))
  plot
}


# Define the statistic function for bootstrapping
bootstrap_stat_wt <- function(data, indices, region) {
  # Draw bootstrap samples
  sample_data <- data[indices, ]
  # Calculate the silent rate for the bootstrap sample
  
  gene_windows$wt<-sites_in_features(gene_windows, sample_data, mode="counts")$counts
  
  # Sum and calculate pct
  wt<-sum(gene_windows[RELATIVEPOS==region]$wt)/sum(gene_windows[RELATIVEPOS==region]$LENGTH)
  
  return(wt)
}


boot_H3K4me1<-function(data, indices){
  
  sample_data <- data[indices, ]
  
  chr_filt$MSH6<-mutations_in_features(chr_filt, sample_data[ trt=="MSH6"])
  chr_filt$WT<-mutations_in_features(chr_filt, sample_data[ trt=="WT"])
  summary<-chr_filt[,.(MSH6=sum(MSH6), WT=sum(WT), ratio=sum(MSH6)/sum(WT), H3K4me1_enrich=mean(`H3K4me1_enrich`), length=sum(stop-start)), by=.(group)][order(group)]
  return(summary$ratio)
  
}


make_gene_windows2<-function(data, window=150){
  deciles<-3000/window
  windows<-rbindlist(apply(data, 1, function(x) {
    
    chr=x["chr"]
    body_starts=seq(as.numeric(x["start"]), as.numeric(x["stop"]), by=3000/deciles);body_starts<-body_starts[-length(body_starts)]
    body_stops<-seq(as.numeric(x["start"]), as.numeric(x["stop"]),  by=3000/deciles)[-1]
    upstream_starts<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-(deciles+1)]
    upstream_stops<-seq(as.numeric(x["start"])-3000, as.numeric(x["start"]), length.out=deciles+1)[-1]
    downstream_starts<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-(deciles+1)]
    downstream_stops<-seq(as.numeric(x["stop"]), as.numeric(x["stop"])+3000, length.out=deciles+1)[-1]
    
    out<-data.table(chr=x["chr"], 
                    start=c(upstream_starts, body_starts, downstream_starts),
                    stop=c(upstream_stops, body_stops, downstream_stops),
                    region=c(rep("upstream", length(upstream_starts)),rep("gene body", length(body_starts)),rep("downstream", length(downstream_starts))),
                    ID=x["ID"], 
                    gene=x["locus"])
    out$pos<-1:nrow(out)
    out$length<-out$stop-out$start
    return(out)
    
  }))
  windows$window_ID<-1:nrow(windows)
  setkey(windows, chr, start, stop)
  gene_windows<-windows
  return(gene_windows)
}


