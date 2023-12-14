library(boot)
# Since your silent rate calculation is a bit complex, we'll define a custom function
silent_rate <- function(data) {
  data<-data[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), CDS_mutation=sum(CDS_mutation), meanpnps=mean(log(pnps)))]
  
  
  silent<-data$syn_muts/data$CDS + (data$mut - data$syn_muts - data$ns_mut) / (data$length - data$CDS)
  silent<-(data$mut - data$CDS_mutation) / (data$length - data$CDS)
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

# # Perform the bootstrapping
# CI<-rbindlist(lapply(1:max(genes_filt_99$enrichment), function(i){
#   set.seed(123) # For reproducibility
#   boot_results <- boot(data=genes_filt_99[enrichment==i,], statistic=bootstrap_stat, R=100)
#   
#   # Calculate the 95% confidence intervals
#   boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
#   dt<-data.table(t(boot_ci))
#   colnames(dt)<-c("upper","lower")
#   return(dt)
# }))
# 
# # Add the lower and upper CI as new columns in your gene_sums dataframe
# gene_sums$lower_ci <- CI[,2]
# gene_sums$upper_ci <- CI[,1]
# 

bootstrap_mutation_silent<-function(data, indices){
  sample_data <- data.table(data[indices, ])
  
  genes_filt_99$syn_mut<-sites_in_features(genes_filt_99, sample_data[MutationType=="synonymous"], mode = "counts")$counts
  
  genes_filt_99$ns_mut<-sites_in_features(genes_filt_99, sample_data[MutationType=="non-synonymous"], mode = "counts")$counts
  
  gene_sums<-genes_filt_99[,.(mut=sum(mutations),syn_muts=sum(syn_mut),ns_mut=sum(ns_mut), H3K4me1=sum(H3K4me1), enrich=mean(enrich, na.rm=T), CDS=sum(CDS), N=.N, se_enrich=sd(enrich, na.rm=T)/sqrt(.N),  length=sum(LENGTH), meanpnps=mean(log(pnps))), by=.(group)][order(group)]
  gene_sums$silent<-gene_sums$syn_muts/gene_sums$CDS+(gene_sums$mut-gene_sums$syn_muts-gene_sums$ns_mut)/(gene_sums$length-gene_sums$CDS)
  
  return(gene_sums$silent)
}



