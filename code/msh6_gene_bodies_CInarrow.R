library(polymorphology2)


# Windows genes -----------------------------------------------------------
table(mutations$type, mutations$trt)

mutations<-fread("data/Strelka2_mutations.csv")
genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genes[,ID:=1:nrow(genes)]
genes$DIRECTION<-genes$direction

gene_windows<-feature_windows(genes[enrich_H3K4me1>0.8322], breaks = 1, dist = 2000, directed = T, IDcol = "gene")

gene_windows$wt<-sites_in_features(gene_windows, mutations[trt=="WT"], mode="counts")$counts
gene_windows$msh6<-sites_in_features(gene_windows, mutations[trt=="MSH6"], mode="counts")$counts

body_vs<-gene_windows[,.(wt=sum(wt), msh6=sum(msh6)), by=.(REGION=="gene body")]
chisq.test(body_vs[,2:3])

# Melt the data
melted <- melt(gene_windows, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


ggplot(sums, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"))


mapp<-fread("data/genmap_TAIR/TAIR10_chr_all.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp[,CHROM:=gsub("Chr","", CHROM)]



gene_windows_mapp<-features_in_features(gene_windows,mapp, mode = "mean", value = "MAPPABILITY" )

gene_windows$MAPPABILITY<-gene_windows_mapp$mean[match(gene_windows_mapp$ID, gene_windows_mapp$ID)]

gene_windows2<-gene_windows[MAPPABILITY==1]

# Melt the data
melted <- melt(gene_windows2, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")

# Sum and calculate pct
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


dcast_sums<-dcast(sums, RELATIVEPOS+LENGTH~geno, value.var = "sum")
chisq.test(dcast_sums[,3:4])
ggplot(dcast_sums, aes(x=RELATIVEPOS, y=msh6/wt))+
  geom_bar(stat="identity")

# Perform the bootstrapping
CI_wt<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="WT",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"wt"
  dt$mean<-mean(boot_results$t[,1])
  return(dt)
}))


# Perform the bootstrapping
CI_msh6<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="MSH6",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"msh6"
  dt$mean<-mean(boot_results$t[,1])
  return(dt)
}))

CI<-rbind(CI_msh6, CI_wt)

CI_merge<-merge(sums, CI, by=c("RELATIVEPOS","geno"))

CI_merge$geno<-factor(CI_merge$geno, levels=c("wt","msh6"), labels=c("Wildtype","msh6"))

pdf("figures/MSH6_genebodies_all_CI_highH3K4me1_narrow.pdf", width=1.5, height=2)
ggplot(CI_merge, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_bar(stat="identity", position="dodge", col="black")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3), labels=c("-2kb","Gene bodies","+2kb"))

dev.off()


library(polymorphology2)



# low ---------------------------------------------------------------------


# Windows genes -----------------------------------------------------------
table(mutations$type, mutations$trt)

mutations<-fread("data/Strelka2_mutations.csv")
genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genes[,ID:=1:nrow(genes)]
genes$DIRECTION<-genes$direction

gene_windows<-feature_windows(genes[enrich_H3K4me1<0.8322], breaks = 1, dist = 2000, directed = T, IDcol = "gene")

gene_windows$wt<-sites_in_features(gene_windows, mutations[trt=="WT"], mode="counts")$counts
gene_windows$msh6<-sites_in_features(gene_windows, mutations[trt=="MSH6"], mode="counts")$counts

body_vs<-gene_windows[,.(wt=sum(wt), msh6=sum(msh6)), by=.(REGION=="gene body")]
chisq.test(body_vs[,2:3])

# Melt the data
melted <- melt(gene_windows, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


ggplot(sums, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"))


mapp<-fread("data/genmap_TAIR/TAIR10_chr_all.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp[,CHROM:=gsub("Chr","", CHROM)]

# 
# enrich1<-features_chip_enrich(gene_windows, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
#                                           inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")
# enrich2<-features_chip_enrich(gene_windows, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
#                                           inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")
# 
# gene_windows$H3K4me1<-(enrich1$enrich+enrich2$enrich)/2

gene_windows_mapp<-features_in_features(gene_windows,mapp, mode = "mean", value = "MAPPABILITY" )

gene_windows$MAPPABILITY<-gene_windows_mapp$mean[match(gene_windows_mapp$ID, gene_windows_mapp$ID)]

gene_windows2<-gene_windows[MAPPABILITY==1]

# Melt the data
melted <- melt(gene_windows2, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")

# Sum and calculate pct
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


dcast_sums<-dcast(sums, RELATIVEPOS+LENGTH~geno, value.var = "sum")
chisq.test(dcast_sums[,3:4])
ggplot(dcast_sums, aes(x=RELATIVEPOS, y=msh6/wt))+
  geom_bar(stat="identity")

# Perform the bootstrapping
CI_wt<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="WT",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"wt"
  return(dt)
}))


# Perform the bootstrapping
CI_msh6<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="MSH6",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"msh6"
  return(dt)
}))

CI<-rbind(CI_msh6, CI_wt)

CI_merge<-merge(sums, CI, by=c("RELATIVEPOS","geno"))

CI_merge$geno<-factor(CI_merge$geno, levels=c("wt","msh6"), labels=c("Wildtype","msh6"))

pdf("figures/MSH6_genebodies_all_CI_lowH3K4me1_narrow.pdf", width=1.5, height=2)
ggplot(CI_merge, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_bar(stat="identity", position="dodge", col="black")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3), labels=c("-2kb","Gene bodies","+2kb"))

dev.off()




# all ---------------------------------------------------------------------

gene_windows<-feature_windows(genes, breaks = 1, dist = 2000, directed = T, IDcol = "gene")

gene_windows$wt<-sites_in_features(gene_windows, mutations[trt=="WT"], mode="counts")$counts
gene_windows$msh6<-sites_in_features(gene_windows, mutations[trt=="MSH6"], mode="counts")$counts

body_vs<-gene_windows[,.(wt=sum(wt), msh6=sum(msh6)), by=.(REGION=="gene body")]
chisq.test(body_vs[,2:3])

# Melt the data
melted <- melt(gene_windows, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


ggplot(sums, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"))


mapp<-fread("data/genmap_TAIR/TAIR10_chr_all.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp[,CHROM:=gsub("Chr","", CHROM)]



gene_windows_mapp<-features_in_features(gene_windows,mapp, mode = "mean", value = "MAPPABILITY" )

gene_windows$MAPPABILITY<-gene_windows_mapp$mean[match(gene_windows_mapp$ID, gene_windows_mapp$ID)]

gene_windows2<-gene_windows[MAPPABILITY==1]

# Melt the data
melted <- melt(gene_windows2, id.vars = c("RELATIVEPOS", "LENGTH"), measure.vars = c("wt", "msh6"), variable.name = "geno", value.name = "value")

# Sum and calculate pct
sums <- melted[, .(sum=sum(value), pct = sum(value) / sum(LENGTH), LENGTH=sum(LENGTH)), by = .(RELATIVEPOS, geno)]


dcast_sums<-dcast(sums, RELATIVEPOS+LENGTH~geno, value.var = "sum")
chisq.test(dcast_sums[,3:4])
ggplot(dcast_sums, aes(x=RELATIVEPOS, y=msh6/wt))+
  geom_bar(stat="identity")

# Perform the bootstrapping
CI_wt<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="WT",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"wt"
  return(dt)
}))


# Perform the bootstrapping
CI_msh6<-rbindlist(lapply(1:max(gene_windows2$RELATIVEPOS), function(i){
  set.seed(123) # For reproducibility
  boot_results <- boot(data=mutations[mutations$trt=="MSH6",], 
                       statistic=function(data, indices) {
                         bootstrap_stat_wt(data, indices, region=i)
                       }, R=100)  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$RELATIVEPOS<-i
  dt$geno<-"msh6"
  return(dt)
}))

CI<-rbind(CI_msh6, CI_wt)

CI_merge<-merge(sums, CI, by=c("RELATIVEPOS","geno"))

CI_merge$geno<-factor(CI_merge$geno, levels=c("wt","msh6"), labels=c("Wildtype","msh6"))

pdf("figures/MSH6_genebodies_all_CI_allH3K4me1_narrow.pdf", width=1.5, height=2)
ggplot(CI_merge, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_bar(stat="identity", position="dodge", col="black")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3), labels=c("-2kb","Gene bodies","+2kb"))

dev.off()

