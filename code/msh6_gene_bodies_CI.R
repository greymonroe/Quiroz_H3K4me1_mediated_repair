library(polymorphology2)


# Windows genes -----------------------------------------------------------
table(mutations$type, mutations$trt)

mutations<-fread("data/Strelka2_mutations.csv")
genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
genes[,ID:=1:nrow(genes)]
genes$DIRECTION<-genes$direction

gene_windows<-feature_windows(genes, breaks = 2, dist = 2000, directed = T, IDcol = "gene")

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

pdf("figures/MSH6_genebodies_all_CI_H3K4me1.pdf", width=2.5, height=2)
ggplot(CI_merge, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_bar(stat="identity", position="dodge", col="black")+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_continuous(name="Mutations/bp", breaks=c(5e-7,1e-6, 2e-6, 4e-6, 6e-6), labels=real_sci_format(c(5e-7,1e-6, 2e-6, 4e-6, 6e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"))

ggplot(CI_merge, aes(x=RELATIVEPOS, y=pct, fill=geno))+
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", size=0.25)+
  geom_point(shape=21)+
  geom_line(aes(group=1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.1)+
  facet_wrap(.~geno, scales="free")+
  scale_fill_manual(values=c(Wildtype="dodgerblue",msh6="orange"), guide="none")+
  scale_y_log10(name="Mutations/bp", breaks=c(1e-6, 3e-6), labels=real_sci_format(c(1e-6, 3e-6)))+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank())+
  scale_x_continuous(breaks=c(1,2,3.5,5,6), labels=c("-2kb","-1kb","Gene bodies","+1kb","+2kb"))

dev.off()


chr<-rbindlist(lapply(list.files("data/Arabidopsis Epigenome Profiles/", full.names = T, pattern="*200*"), function(f){
  chr<-fread(f)
  chr$CHROM<-as.character(chr$chr)
  chr$START<-chr$start
  chr$STOP<-chr$stop
  return(chr)
}))
setkey(chr, CHROM, start, stop)
chr$ID<-1:nrow(chr)

mapp<-fread("data/genmap_TAIR/TAIR10_chr_all.genmap.bedgraph")
colnames(mapp)<-c("CHROM","START","STOP","MAPPABILITY")
mapp[,CHROM:=gsub("Chr","", CHROM)]

chr_mapp<-features_in_features(chr,mapp, mode = "mean", value = "MAPPABILITY" )

chr$MAPPABILITY<-chr_mapp$mean[match(chr_mapp$ID, chr_mapp$ID)]

chr$MSH6<-mutations_in_features(chr, mutations[ trt=="MSH6"])
chr$WT<-mutations_in_features(chr, mutations[ trt=="WT"])

enrich1<-features_chip_enrich(chr, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
                                          inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")
enrich2<-features_chip_enrich(chr, chipfile = "data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
                                          inputfile = "data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")

chr$H3K4me1_enrich<-(enrich1$enrich+enrich2$enrich)/2
chr$input<-(enrich1$input+enrich2$input)
chr$chip<-(enrich1$chip+enrich2$chip)


chr_filt <- chr[MAPPABILITY==1, ]

chr_filt$group=as.numeric(cut2(chr_filt$H3K4me1_enrich, g=10))
summary<-chr_filt[,.(MSH6=sum(MSH6), WT=sum(WT), ratio=sum(MSH6)/sum(WT), H3K4me1_enrich=mean(`H3K4me1_enrich`), length=sum(stop-start)), by=.(group)][order(group)]
summary$MSH6_pct<-summary$MSH6/summary$length
summary$WT_pct<-summary$WT/summary$length
summary
chisq.test(summary[,2:3])

pdf("figures/msh6_H3K4me1_mut_diff.pdf", width=2, height=2)
ggplot(summary, aes(x=H3K4me1_enrich, y=log(ratio), fill=group*10))+
  geom_smooth(method = "lm", alpha=0.25, col="gray") +
  geom_point(shape=21, col="black")+
  scale_fill_gradient2(
    name = "%ile", 
    low = "gray90", 
    mid = "palegreen", 
    high = "green4", 
    midpoint = median(chr_filt$group)*10,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  )+
  theme_classic(base_size = 6)+
  theme(
    legend.position = c(.85, .25),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6),
    legend.background = element_blank()
  )+
 
  scale_y_continuous(name="log(MSH6:WT)")+
  scale_x_continuous(name="H3K4me1 enrichment")
  

long_summary <- melt(summary, id.vars = c("group", "length","H3K4me1_enrich"), measure.vars = c("MSH6", "WT"))

ggplot(long_summary, aes(x = H3K4me1_enrich, y = (value/length), fill=variable, color = variable, group = variable)) +
  geom_smooth(method = "lm", alpha=0.25) +
  geom_point(shape = 21, col="black") +
  scale_color_manual(values = c("WT" = "dodgerblue", "MSH6" = "orange"), name="") +
  scale_fill_manual(values = c("WT" = "dodgerblue", "MSH6" = "orange"), name="") +
  theme_classic(base_size = 6) +
  scale_y_log10(name="Mutations/bp", breaks=c(5e-6, 1e-6, 2e-6, 5e-7), labels=real_sci_format(c(5e-6, 1e-6, 2e-6, 5e-7)))+
  scale_x_continuous(name="H3K4me1 enrichment")+
  theme(
    legend.position = c(.25, .25),
    legend.key.size = unit(0.3, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6),
    legend.background = element_blank()
  )
dev.off()

pdf("figures/Athal_H3K4me1windows.pdf", width=1, height=0.6)

ggplot(chr_filt, aes(H3K4me1_enrich, fill=group*10, group=group))+
  geom_histogram(col="black", size=0.25)+
  scale_fill_gradient2(
    name = "%ile", 
    low = "gray90", 
    mid = "palegreen", 
    high = "green4", 
    midpoint = median(chr_filt$group)*10,
    guide = "none"
  )+
  scale_x_continuous(name="H3K4me1")+
  theme_classic(base_size = 6)+
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = c(.25, .75),
    legend.key.size = unit(0.2, 'cm'), # Adjust size as needed
    legend.text = element_text(size = 6), # Adjust text size as needed
    legend.key = element_rect(colour = "black", linewidth = 0.5),
  )+
  scale_y_continuous(name="N", breaks=c(0, 30000, 60000), labels=c("0","30K","60K"))

dev.off()




boot_results <- boot(data=mutations, statistic=boot_H3K4me1, R=100)


CI<-rbindlist(lapply(1:max(chr$group), function(i){
  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, index=i, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$group=i
  return(dt)
}))

summary$lower_ci <- CI[,2]
summary$upper_ci <- CI[,1]



boot_loc<-function(data, indices){
  sample_data <- data[indices, ]
  
  table(loc=sample_data$loc, trt=sample_data$trt)
  summary<-data.table(table(loc=sample_data$loc, trt=sample_data$trt))
  summary <- dcast(summary, loc ~ trt, value.var = "N")
  summary$ratio<-summary$MSH6/summary$WT
  return((summary$ratio))
}


summary<-data.table(table(loc=mutations$loc, trt=mutations$trt))
summary <- dcast(summary, loc ~ trt, value.var = "N")
summary$ratio<-summary$MSH6/summary$WT
chisq.test(summary[,c(2:3)])

boot_results <- boot(data=mutations, statistic=boot_loc, R=1000)
boot_ci <- boot.ci(boot_results,type="perc")$percent[4:5]
dt<-data.table(t(boot_ci))
colnames(dt)<-c("upper","lower")


CI<-rbindlist(lapply(1:3, function(i){
  
  # Calculate the 95% confidence intervals
  boot_ci <- boot.ci(boot_results, index=i, type="perc")$percent[4:5]
  dt<-data.table(t(boot_ci))
  colnames(dt)<-c("upper","lower")
  dt$group=i
  return(dt)
}))

summary$lower_ci <- CI[,2]
summary$upper_ci <- CI[,1]

pdf("Figures/MSH6_loc_mutations_CI.pdf", width=1, height=2)
summary$loc<-factor(summary$loc, levels=c("TE","intergenic","genic"))
ggplot(summary, aes(x=loc, y=log(ratio), fill=as.character(loc)))+
  geom_bar(stat="identity", position="dodge", col="black")+
  scale_fill_manual(values=rev(c("gray90","palegreen","green4")), guide="none")+
  scale_y_continuous(name="log(msh6:WT mutations)")+
  geom_errorbar(aes(ymin=log(lower_ci), ymax=log(upper_ci)), width=0.1)+
  theme_classic(base_size = 6)+
  scale_x_discrete(name="Location", labels=c("TE", "Inter\ngenic", "Gene\nbodies"))+
  theme(axis.title.y = element_text(angle=90, hjust=0.5, vjust=0.5))

dev.off()



pdf("figures/SBS_mutations.pdf", width=6, height=1.25)

plot_tricontexts(mutations[ trt=="WT"]$context_SBS)
plot_tricontexts(mutations[type=="somatic" & trt=="WT"]$context_SBS)$plot+ggtitle("Wildtype somatic")
plot_tricontexts(mutations[type!="somatic" & trt=="WT"]$context_SBS)$plot+ggtitle("Wildtype germline")

plot_tricontexts(mutations[trt=="MSH6"]$context_SBS)
plot_tricontexts(mutations[type=="somatic" & trt=="MSH6"]$context_SBS)$plot+ggtitle("msh6 somatic")
plot_tricontexts(mutations[type!="somatic" & trt=="MSH6"]$context_SBS)$plot+ggtitle("msh6 germline")

dev.off()



pdf("figures/SBS_mutations_simple.pdf", width=1, height=1.25)

plot_tricontexts(mutations[ trt=="WT"]$context_SBS, full=F)
plot_tricontexts(mutations[type=="somatic" & trt=="WT"]$context_SBS, full=F)$plot+ggtitle("WT\nsomatic")
plot_tricontexts(mutations[type!="somatic" & trt=="WT"]$context_SBS, full=F)$plot+ggtitle("WT\ngermline")
plot_tricontexts(mutations[trt=="MSH6"]$context_SBS, full=F)
plot_tricontexts(mutations[type=="somatic" & trt=="MSH6"]$context_SBS, full=F)$plot+ggtitle("msh6\nsomatic")
plot_tricontexts(mutations[type!="somatic" & trt=="MSH6"]$context_SBS, full=F)$plot+ggtitle("msh6\ngermline")

dev.off()

pdf("figures/SBS_mutations_simple_lines_somatic.pdf", width=7, height=4)

meta<-data.table(table(file=mutations$file, genotype=mutations$genotype, line=mutations$line))[N>0][order(genotype, decreasing = T)]
line_contexts<-lapply(meta$file, function(f){
  muts<-mutations[file==f & type=="somatic"]
  genotype<-unique(muts$genotype)
  line<-unique(muts$line)
  plot_tricontexts(muts$context_SBS, full=F)$plot+ggtitle(label = f, subtitle = paste(genotype,line,sep="\n"))+
    scale_y_continuous(limits=c(0,30), name="N")
  
})

library(gridExtra)
grobs <- lapply(line_contexts, ggplotGrob)
do.call("grid.arrange", c(line_contexts, ncol = 7))
dev.off()

pdf("figures/SBS_mutations_simple_lines.pdf", width=7, height=4)

meta<-data.table(table(file=mutations$file, genotype=mutations$genotype, line=mutations$line))[N>0][order(genotype, decreasing = T)]
line_contexts<-lapply(meta$file, function(f){
  muts<-mutations[file==f]
  genotype<-unique(muts$genotype)
  line<-unique(muts$line)
  plot_tricontexts(muts$context_SBS, full=F)$plot+ggtitle(label = f, subtitle = paste(genotype,line,sep="\n"))+
    scale_y_continuous(limits=c(0,50), name="N")
  
})

library(gridExtra)
grobs <- lapply(line_contexts, ggplotGrob)
do.call("grid.arrange", c(line_contexts, ncol = 7))
dev.off()

