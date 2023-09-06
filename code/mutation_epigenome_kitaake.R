source("code/libraries_functions.R")
source("code/load_rice_data.R")

pdf("figures/genome.pdf", width=5.7, height=1.5)
    genome_windows<-make_genome_windows(genome[1:12], 3000000, overlap=F)
    genome_windows_scale<-data.table(scale(genome_windows[,-c(1:4)]))
    genome_windows_melt<-melt(cbind(genome_windows[,c(1,2)], genome_windows_scale), id.vars = c("chr", "start"))
    starts<-genome_windows[start==1]
    ggplot(genome_windows_melt, aes(x=start/1e6, group=chr, y=variable, fill=value))+
      geom_tile()+
      scale_fill_gradientn(colors=c("orange","white","blue","black"))+
      theme_bw(base_size = 6)+
      facet_grid(~chr, scales = "free", space = "free_x")+
      theme(legend.key.size = unit(.5,"line"), legend.key=element_rect(color="black"))+
      scale_x_continuous(name="Mb")
dev.off()

pdf("figures/genome_window_regression.pdf", width=1.7, height=1.5)

  genome_windows_log<-make_genome_windows(genome[1:12], 100, overlap=T)
  model_sum<-log_model(genome_windows_log, "snp")
  plot_model(model_sum, "SBS, 100bp windows, log")
  model_sum<-log_model(genome_windows_log, "snp_noCT")
  plot_model(model_sum, "non-C>T SBS, 100bp windows, log")
  model_sum<-log_model(genome_windows_log, "snp_homoz")
  plot_model(model_sum, "Homozygous SBS, 100bp windows, log")
  model_sum<-log_model(genome_windows_log, "indel")
  plot_model(model_sum, "InDel, 100bp windows, log")
  
  model_sum<-log_model(genome_windows_log, "snp", aic=T)
  plot_model(model_sum, "SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(genome_windows_log, "snp_noCT", aic=T)
  plot_model(model_sum, "non-C>T SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(genome_windows_log, "snp_homoz", aic=T)
  plot_model(model_sum, "Homozygous SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(genome_windows_log, "indel", aic=T)
  plot_model(model_sum, "InDel, 100bp windows, \nlog aic")
  
  model_sum<-log_model_single(genome_windows_log, "snp")
  plot_model(model_sum, "SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(genome_windows_log, "snp_noCT")
  plot_model(model_sum, "non-C>T SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(genome_windows_log, "snp_homoz")
  plot_model(model_sum, "Homozygous SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(genome_windows_log, "indel")
  plot_model(model_sum, "InDel, 100bp windows, \nlog single")
  
  genome_windows<-make_genome_windows(genome[1:12], 100, overlap=F)
  model_sum<-lm_model(genome_windows, "snp", aic=F)
  plot_model(model_sum, "SBS, 100bp windows, lm")
  model_sum<-lm_model(genome_windows, "snp_noCT, lm", aic=F)
  plot_model(model_sum, "non-C>T SBS, 100bp windows, lm")
  model_sum<-lm_model(genome_windows, "snp_homoz", aic=F)
  plot_model(model_sum, "non-C>T SBS, 100bp windows, lm")

dev.off()

gene_windows<-make_gene_windows(data = gene_annotations_basic, window=100)
gene_windows$snp<-add_vars_to_gene_windows(gene_windows, snps)
gene_windows$snp_noCT<-add_vars_to_gene_windows(gene_windows, snps[!Single.base.substitution %in% c("C>T","G>A")])
gene_windows$snp_homoz<-add_vars_to_gene_windows(gene_windows, snps[Genotype=="Homozygous"])
gene_windows$indel<-add_vars_to_gene_windows(gene_windows, indels)
gene_windows$snp_count<-add_vars_hits_to_gene_windows(gene_windows, snps)
gene_windows$snp_homopolymer_free<-add_vars_to_gene_windows(gene_windows, snps[homopolymer_neighbor!=T])

pdf("figures/gene_windows_regression.pdf", width=1.8, height=1.5)
 ##### all variable pred
  model_sum<-log_model(gene_windows, "snp")
  plot_model(model_sum, "SBS")
  model_sum<-log_model(gene_windows, "snp_homopolymer_free")
  plot_model(model_sum, "SBS homopolymers removed")
  model_sum<-log_model(gene_windows, "snp_homoz")
  plot_model(model_sum, "Homozygous SBS")
  model_sum<-log_model(gene_windows, "snp_noCT")
  plot_model(model_sum, "non-C>T SBS")
  model_sum<-log_model(gene_windows[gene %in% LoF$Gene.ID], "snp_noCT")
  plot_model(model_sum, "non-C>T SBS, LoF only")
  model_sum<-log_model(gene_windows[gene %in% LoF$Gene.ID], "snp")
  plot_model(model_sum, "SBS, LoF only")
  model_sum<-log_model(gene_windows, "indel")
  plot_model(model_sum, "indel")
  ##### aic variable pred
  model_sum<-log_model(gene_windows, "snp", aic=T)
  plot_model(model_sum, "SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(gene_windows, "snp_noCT", aic=T)
  plot_model(model_sum, "non-C>T SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(gene_windows, "snp_homoz", aic=T)
  plot_model(model_sum, "Homozygous SBS, 100bp windows, \nlog aic")
  model_sum<-log_model(gene_windows[gene %in% LoF$Gene.ID], "snp", aic=T)
  plot_model(model_sum, "SBS, LoF only, \nlog aic")
  ##### single variable pred
  model_sum<-log_model_single(gene_windows, "snp")
  plot_model(model_sum, "SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(gene_windows, "snp_noCT")
  plot_model(model_sum, "non-C>T SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(gene_windows, "snp_homoz")
  plot_model(model_sum, "Homozygous SBS, 100bp windows, \nlog single")
  model_sum<-log_model_single(gene_windows[gene %in% LoF$Gene.ID], "snp")
  plot_model(model_sum, "SBS, LoF only, \nlog single")
  ##### Additional models
  model_sum<-lm_model(gene_windows, "snp_count", aic=T)
  plot_model(model_sum, "SBS, 100bp windows, \nlm aic")
  model_sum<-lm_model(gene_windows, "snp_count", aic=F)
  plot_model(model_sum, "SBS, 100bp windows, \nlm")
  model_sum<-poisson_model(gene_windows, "snp_count", aic=T)
  plot_model(model_sum, "SBS, 100bp windows, \npoisson aic")
  model_sum<-poisson_model(gene_windows, "snp_count", aic=F)
  plot_model(model_sum, "SBS, 100bp windows, \npoisson")
dev.off()

pdf("figures/gene_windows_regression_all_SBS.pdf", width=1.8, height=1.5)

model_sum<-lm_model(gene_windows, "snp_count", aic=T)
plot_model(model_sum, "SBS, 100bp windows, \nlm aic")
model_sum<-lm_model(gene_windows, "snp_count", aic=F)
plot_model(model_sum, "SBS, 100bp windows, \nlm")
model_sum<-poisson_model(gene_windows, "snp_count", aic=T)
plot_model(model_sum, "SBS, 100bp windows, \npoisson aic")
model_sum<-poisson_model(gene_windows, "snp_count", aic=F)
plot_model(model_sum, "SBS, 100bp windows, \npoisson")
## around peaks

pdf("figures/peaks_mutations.pdf", width=1.7, height=1.5)
  out_H3K4me1_genes<-plot_peaks(gene_annotations_basic, "H3K4me1","Gene body", 30,H3K4me1, gene=T, marky=T)
  print(out_H3K4me1_genes[[1]])
  out_H3K4me1_genes_random<-plot_peaks(gene_annotations_basic, "Random","Gene body", 30,peaks_randomized(H3K4me1), gene=T, marky=T)
  
  rand<-out_H3K4me1_genes_random[[2]]
  rand$pct<-rand$pct/10
  ggplot(out_H3K4me1_genes[[2]], aes(x=pos, y=pct, col=region=="gene body", group=region))+
    geom_line(data=rand, col="black", group=1)+
    geom_vline(xintercept = c(30, 30*2), linetype="dashed", size=0.25)+
    geom_line()+
    scale_color_manual(values=c("gray75","green3"), guide="none")+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="Peaks")+
    ggtitle("H3K4me1")+
    scale_x_continuous(breaks=c(0,max(out_H3K4me1_genes[[2]]$pos)/3, max(out_H3K4me1_genes[[2]]$pos)/3*2, max(out_H3K4me1_genes[[2]]$pos)), labels=c("-2kb","0%","100%","+2kb"), name="Gene body")
  
  out_H3K4me1<-plot_peaks(H3K4me1, "H3K4me1, SBS","Peak region", 30, ns_s)
  print(out_H3K4me1[[1]])
  out_H3K4me1_random<-plot_peaks(peaks_randomized(H3K4me1), "Random, SBS","Region", 30, ns_s)

  ggplot(out_H3K4me1[[2]], aes(x=pos, y=pct, col=region=="gene body", group=region))+
    geom_line(data=out_H3K4me1_random[[2]], col="black", group=1)+
    geom_vline(xintercept = c(30, 30*2), linetype="dashed", size=0.25)+
    geom_line()+
    scale_color_manual(values=c("gray75","green3"), guide="none")+
    theme_classic(base_size = 6)+
    scale_y_continuous(name="Mutations/bp")+
    ggtitle("SBS")+
    scale_x_continuous(breaks=c(0,max(out_H3K4me1[[2]]$pos)/3, max(out_H3K4me1[[2]]$pos)/3*2, max(out_H3K4me1[[2]]$pos)), labels=c("-2kb","0%","100%","+2kb"), name="Peak region")

    out_genebody<-plot_peaks(gene_annotations_basic, "Genes, SBS","Gene body", 30,ns_s, gene=T)
    print(out_genebody[[1]])
    
    plot_peaks(gene_annotations_basic, "Genes, homozygous SBS","Gene body", 30,ns_s[Genotype=="Homozygous"], gene=T)
    plot_peaks(gene_annotations_basic, "Genes, non-C>T SBS","Gene body", 30,ns_s[!`Single base substitution` %in% c("C>T","G>A")], gene=T)
    plot_peaks(gene_annotations_basic[locus %in% LoF$Gene.ID], "LoF Only Genes, SBS","Gene body", 30,ns_s, gene=T)
dev.off()

pdf("figures/peaks_mutations_extra.pdf", width=1.7, height=1.5)

  plot_peaks(H3K4ac, "H3K4ac, SBS","Peak region", 30, ns_s)
  plot_peaks(H3K9me2, "H3K9me2, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K36me3, "H3K36me3, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K4me3, "H3K4me3, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K9me1, "H3K9me1, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K12ac, "H3K12ac, SBS","Peak region", 30,ns_s)
  plot_peaks(PII, "PII, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K27me3, "H3K27me3, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K4me1, "H3K4me1, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K27ac, "H3K27ac, SBS","Peak region", 30,ns_s)
  plot_peaks(H3K9ac, "H3K9ac, SBS","Peak region", 30,ns_s)
  
  plot_peaks(gene_annotations_basic, "PII","Gene body", 30,PII, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K9me2","Gene body", 30,H3K9me2, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K4me3","Gene body", 30,H3K4me3, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K9ac","Gene body", 30,H3K9ac, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K36me3","Gene body", 30,H3K36me3, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K12ac","Gene body", 30,H3K12ac, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K27me3","Gene body", 30,H3K27me3, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K9me1","Gene body", 30,H3K9me1, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K27ac","Gene body", 30,H3K27ac, gene=T, marky=T)
  plot_peaks(gene_annotations_basic, "H3K4me1","Gene body", 30,H3K4me1, gene=T, marky=T)
  
dev.off()


# SBS ---------------------------------------------------------

pdf("figures/kitaake_SBS.pdf", width=6, height=1.5)
  ggplot(context_table, aes(x=context_only, y=(N), fill=mut))+
    geom_bar(stat="identity", width=0.5)+
    facet_grid(~mut, scales = "free")+
    theme_classic(base_size = 6)+
    scale_x_discrete(name="Context")+
    ggtitle("Rice")+
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(x = 0.3, units = "line"))+
    scale_fill_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")
  ggplot(Athal, aes(x=context_only, y=(N), fill=mut))+
    geom_bar(stat="identity", width=0.5)+
    facet_grid(~mut, scales = "free")+
    ggtitle("A. thaliana")+
    theme_classic(base_size = 6)+
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(x = 0.3, units = "line"))+
    scale_fill_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")
dev.off()

pdf("figures/kitaake_Athal_SBS_scatter.pdf", width=1.5, height=1.5)
  SBS_dt<-data.table(A_thal=Athal$N, rice=context_table$N, mut=context_table$mut)
  ggplot(SBS_dt, aes(x=log(A_thal),y=log(rice), col=mut))+
    geom_point(size=0.5)+
    scale_x_continuous(name="log(A. thaliana SBS)")+
    scale_y_continuous(name="log(Kitaake rice SBS)")+
    scale_color_manual(values=c("cyan3","black","red4","gray","green3","pink3"), guide="none")+
    theme_classic(base_size = 6)
dev.off()

pdf("figures/tss_tts_all.pdf", width=3, height=1.5)
  out<-tss_tts.variants(gff = gene_annotations_basic,vcf=snps)
  tss_tts.variants.plot(out, window=200)+ggtitle("SBS")
  out<-tss_tts.variants(gff = gene_annotations_basic[locus %in% LoF$Gene.ID],vcf=snps)
  tss_tts.variants.plot(out, window=200)+ggtitle("SBS LoF only")
  out<-tss_tts.variants(gff = gene_annotations_basic,vcf=snps[Genotype=="Homozygous"])
  tss_tts.variants.plot(out, window=200)+ggtitle("SBS Homozygous only")
  out<-tss_tts.variants(gff = gene_annotations_basic,vcf=snps[!Single.base.substitution %in% c("C>T","G>A")])
  tss_tts.variants.plot(out, window=200)+ggtitle("SBS non-C>T only")
  out<-tss_tts.variants(gff = gene_annotations_basic,vcf=indels)
  tss_tts.variants.plot(out, window=200)+ggtitle("InDels")
dev.off()


# NS/S -------------------------------------------------------
          
      mut_table<-data.table(table(ref=s$Old, alt=s$New))[N>0]
      write.table(mut_table, "data/mut_table.txt")
      null<-Null_ns_s(mutations = "data/mut_table.txt", composition = "data/Base_pair_frequency_of_CDS.csv")
      
      stop=sum(ns_s$`StopCodon?`=="Stop_codon")
      Ns=sum(ns_s$MutationType=="non-synonymous")-stop
      S=sum(ns_s$MutationType=="synonymous")
      
      ns_s_genes<-ns_s[ns_s$gene %in% gene_annotations_basic$locus & `StopCodon?`!="Stop_codon"]
      Ns_gene=sum(ns_s_genes$MutationType=="non-synonymous")
      S_gene=sum(ns_s_genes$MutationType=="synonymous")
      
      ns_s_TE<-ns_s[ns_s$gene %in% TE$locus & `StopCodon?`!="Stop_codon"]
      Ns_TE=sum(ns_s_TE$MutationType=="non-synonymous")
      S_TE=sum(ns_s_TE$MutationType=="synonymous")
      
      gene_annotations_all$syn<-add_vars_hits_to_gene_windows(gene_annotations_all, synonymous)
      gene_annotations_all$missense<-add_vars_hits_to_gene_windows(gene_annotations_all, missense)
      gene_annotations_all$ns<-add_vars_hits_to_gene_windows(gene_annotations_all, ns_s[MutationType=="non-synonymous" &`StopCodon?`!="Stop_codon"])
      gene_annotations_all$s<-add_vars_hits_to_gene_windows(gene_annotations_all, ns_s[MutationType=="synonymous"])
      gene_annotations_all$pnps<-gene_annotations_all$missense/gene_annotations_all$syn
      gene_annotations_all_means<-gene_annotations_all[,.(syn=sum(syn), missense=sum(missense), ns_s=sum(missense)/sum(syn)), by=.(is_TE)]
      
      chisq.test(matrix(c(Ns_TE, S_TE, Ns_gene, S_gene), nrow=2))
      chisq.test(c(Ns_gene, S_gene), p = prop.table(c(null[[1]], 1)))
      chisq.test(matrix(c(Ns, S, nrow(missense), nrow(synonymous)), nrow=2))
      
      dt<-data.table(ratio=c(null[[1]], Ns/S, Ns_TE/S_TE, Ns_gene/S_gene,  nrow(missense)/nrow(synonymous), gene_annotations_all_means$ns_s[1], gene_annotations_all_means$ns_s[2]), 
                     src=c("Random expectation","MA KitaakeX all","MA KitaakeX \nprotein coding genes", "MA KitaakeX \ntransposable elements", "3,010 accessions all","3,010 accessions \nprotein coding genes", "3,010 accessions \ntransposable elements"),
                     Data=c("Neutral simulation","Mutation accumulation (MA)","Mutation accumulation (MA)","Mutation accumulation (MA)", "Natural populations","Natural populations","Natural populations"))
      
      pdf("figures/Non-synonymous.pdf", width=3.5, height=2)
      ggplot(dt, aes(x=src, y=ratio, fill=Data))+
        geom_bar(stat="identity", position='dodge', col="black", width=0.5)+
        theme_classic(base_size = 6)+
        scale_y_continuous(name="Non-syn/Syn")+
        scale_x_discrete(name="")+
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.4, units="line"))
      dev.off()
      
      table(ns_s$gene=="")
      gene_annotations_all$CHROM<-gene_annotations_all$chr
      gene_annotations_all$START<-gene_annotations_all$start
      gene_annotations_all$STOP<-gene_annotations_all$stop
      ns_s$CHROM<-ns_s$chr
      ns_s$POSITION<-ns_s$POS
      ns_s$genic<-features_overlap_mutation(gene_annotations_all[is_TE=="N" & is_representative=="Y"], ns_s)
      ns_s$TE<-features_overlap_mutation(gene_annotations_all[is_TE=="Y" & is_representative=="Y"], ns_s)
      ns_s$loc<-ifelse(ns_s$genic, "genic","intergenic")
      ns_s$loc[ns_s$TE]<-"TE"
     
      gene_annotations_all[is_representative=="Y",.(sum(length)), by=is_TE]
      table(ns_s$genic)
      
      pdf("figures/mutation_gene_intergene_TE.pdf", width=1.5, height=1.5)
      sum<-data.table(table(loc=ns_s$loc))
      sum$length<-c(110152630, 374471240-110152630-54427496, 54427496)
      sum$pct<-sum$N/sum$length
      chisq.test(sum[,2:3])
      ggplot(sum, aes(x=loc, y=pct))+
        geom_bar(stat="identity", col="black", fill="dodgerblue4")+
        scale_x_discrete(name="Location")+
        scale_y_continuous(name="Mutations/bp")+
        theme_classic(base_size = 6)
      sum<-data.table(table(loc=ns_s[!`Single base substitution` %in% c("C>T","G>A")]$loc))
      sum$length<-c(110152630, 374471240-110152630-54427496, 54427496)
      sum$pct<-sum$N/sum$length
      chisq.test(sum[,2:3])
      ggplot(sum, aes(x=loc, y=pct))+
        geom_bar(stat="identity", col="black", fill="dodgerblue4")+
        scale_x_discrete(name="Location")+
        scale_y_continuous(name="Mutations (non C>T)/bp")+
        theme_classic(base_size = 6)
      dev.off()
      
      
      
# rice constraints
gene_annotations_all$mut<-add_vars_hits_to_gene_windows(gene_annotations_all, ns_s)
gene_annotations_all$mut_nonCT<-add_vars_hits_to_gene_windows(gene_annotations_all,ns_s[!`Single base substitution` %in% c("C>T","G>A")])
gene_annotations_all$H3K4me1<-encode_hits(H3K4me1, gene_annotations_all, out="length")
gene_annotations_all$H3K36me3<-encode_hits(H3K36me3, gene_annotations_all, out="length")
gene_annotations_all$coding<-add_vars_hits_to_gene_windows(gene_annotations_all, ns_s[WithinCDS=="Yes"])
gene_annotations_all$noncoding<-add_vars_hits_to_gene_windows(gene_annotations_all, ns_s[WithinCDS=="No"])

pdf("figures/rice_mutation_constraints.pdf", width=1.5, height=1.5)

# genes with or without H3K4me1
gene_annotations_all_sums<-gene_annotations_all[is_representative=="Y" & chr %in% paste0("Chr",1:12) & is.finite(pnps),.(mut=sum(mut),  length=sum(length), ns=sum(ns), s=sum(s), ns_s=sum(ns)/sum(s), pct=sum(mut)/sum(length), H3K4me1=sum(H3K4me1),nonH3K4me1=sum(H3K4me1==0),H3K4me1_pct=sum(H3K4me1)/sum(length), N=.N, mut_nonCT=sum(mut_nonCT), pct_nonCT=sum(mut_nonCT)/sum(length), mut_syn=sum(s), pct_mut_syn=sum(s)/sum(CDS),  mut_noncoding=sum(noncoding), pct_mut_noncoding=sum(noncoding)/sum(length-CDS), CDS=sum(CDS), noncoding=sum(length-CDS)), by=.(grp=H3K4me1>0)]

chisq.test(gene_annotations_all_sums[,c("mut_syn","CDS"), with=F])
(gene_annotations_all_sums$pct_mut_syn[2]-gene_annotations_all_sums$pct_mut_syn[1])/gene_annotations_all_sums$pct_mut_syn[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_mut_syn", xlab="H3K4me1",ylab="Synonymous Mutations/CDS bp",ggtitle = "p = 2e-4")

chisq.test(gene_annotations_all_sums[,2:3])
(gene_annotations_all_sums$pct[2]-gene_annotations_all_sums$pct[1])/gene_annotations_all_sums$pct[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct", xlab="H3K4me1",ylab="Mutations/bp",ggtitle = "p < 2e-16")

chisq.test(gene_annotations_all_sums[,4:5])
plot_bars_rice(gene_annotations_all_sums, yvar="ns_s", xlab="H3K4me1",ylab="de novo MA N/S",ggtitle = "p = 0.99")

chisq.test(gene_annotations_all_sums[,c(8,3)])
plot_bars_rice(gene_annotations_all_sums, yvar="H3K4me1_pct", xlab="H3K4me1",ylab="H3K4me1 peaks (% of genes)",ggtitle = "p < 2e-16")


gene_annotations_all_sums<-gene_annotations_all[is_representative=="Y" & chr %in% paste0("Chr",1:12) & is.finite(pnps),.(mut=sum(mut),  length=sum(length), ns=sum(ns), s=sum(s), ns_s=sum(ns)/sum(s), pct=sum(mut)/sum(length), H3K4me1=sum(H3K4me1),nonH3K4me1=sum(H3K4me1==0),H3K4me1_pct=sum(H3K4me1)/sum(length), N=.N, mut_nonCT=sum(mut_nonCT), pct_nonCT=sum(mut_nonCT)/sum(length), mut_syn=sum(s), pct_mut_syn=sum(s)/sum(CDS),  mut_noncoding=sum(noncoding), pct_mut_noncoding=sum(noncoding)/sum(length-CDS), CDS=sum(CDS), noncoding=sum(length-CDS)), by=.(grp=is_expressed)]
#all
chisq.test(gene_annotations_all_sums[,2:3])
(gene_annotations_all_sums$pct[2]-gene_annotations_all_sums$pct[1])/gene_annotations_all_sums$pct[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct", xlab="Expressed",ylab="Mutations/bp",ggtitle = "p = 1e-15")
#nonCT
chisq.test(gene_annotations_all_sums[,c("mut_nonCT","length"), with=F])
(gene_annotations_all_sums$pct_nonCT[2]-gene_annotations_all_sums$pct_nonCT[1])/gene_annotations_all_sums$pct_nonCT[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_nonCT", xlab="Expressed",ylab="Mutations/bp",ggtitle = "p = 1e-6")
#synonly
chisq.test(gene_annotations_all_sums[,c("mut_syn","CDS"), with=F])
(gene_annotations_all_sums$pct_mut_syn[2]-gene_annotations_all_sums$pct_mut_syn[1])/gene_annotations_all_sums$pct_mut_syn[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_mut_syn", xlab="Expressed",ylab="Synonymous Mutations/CDS bp",ggtitle = "p = 1e-7")
#non-coding
chisq.test(gene_annotations_all_sums[,c("mut_noncoding","noncoding"), with=F])
(gene_annotations_all_sums$pct_mut_noncoding[2]-gene_annotations_all_sums$pct_mut_noncoding[1])/gene_annotations_all_sums$pct_mut_noncoding[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_mut_noncoding", xlab="PnPs",ylab="non-coding Mutations/bp",ggtitle = "p = 1e-12")

#N/S
chisq.test(gene_annotations_all_sums[,4:5])
plot_bars_rice(gene_annotations_all_sums, yvar="ns_s", xlab="Expressed",ylab="de novo MA N/S",ggtitle = "p = 0.98")

#H3K4me1
chisq.test(gene_annotations_all_sums[,c(8,3)])
plot_bars_rice(gene_annotations_all_sums, yvar="H3K4me1_pct", xlab="Expressed",ylab="H3K4me1 peaks (% of genes)",ggtitle = "p < 2e-16")

#PnPs
gene_annotations_all_sums<-gene_annotations_all[ is_representative=="Y" & chr %in% paste0("Chr",1:12) & is.finite(pnps),.(mut=sum(mut), length=sum(length), ns=sum(ns), s=sum(s), ns_s=sum(ns)/sum(s), pct=sum(mut)/sum(length), H3K4me1=sum(H3K4me1),nonH3K4me1=sum(H3K4me1==0),H3K4me1_pct=sum(H3K4me1)/sum(length), N=.N,mut_nonCT=sum(mut_nonCT), pct_nonCT=sum(mut_nonCT)/sum(length),mut_syn=sum(s), pct_mut_syn=sum(s)/sum(CDS), mut_noncoding=sum(noncoding), pct_mut_noncoding=sum(noncoding)/sum(length-CDS), CDS=sum(CDS), noncoding=sum(length-CDS)), by=.(grp=Hmisc::cut2(pnps, g=2))]
#all
chisq.test(gene_annotations_all_sums[,2:3])
(gene_annotations_all_sums$pct[2]-gene_annotations_all_sums$pct[1])/gene_annotations_all_sums$pct[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct", xlab="PnPs",ylab="Mutations/bp",ggtitle = "p < 2e-16")
#nonCT
chisq.test(gene_annotations_all_sums[,c("mut_nonCT","length"), with=F])
(gene_annotations_all_sums$pct_nonCT[2]-gene_annotations_all_sums$pct_nonCT[1])/gene_annotations_all_sums$pct_nonCT[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_nonCT", xlab="PnPs",ylab="Mutations/bp",ggtitle = "p = 2e-8")
#synonly
chisq.test(gene_annotations_all_sums[,c("mut_syn","CDS"), with=F])
(gene_annotations_all_sums$pct_mut_syn[2]-gene_annotations_all_sums$pct_mut_syn[1])/gene_annotations_all_sums$pct_mut_syn[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_mut_syn", xlab="PnPs",ylab="Synonymous Mutations/CDS bp",ggtitle = "p = 1e-12")
#non-coding
chisq.test(gene_annotations_all_sums[,c("mut_noncoding","noncoding"), with=F])
(gene_annotations_all_sums$pct_mut_noncoding[2]-gene_annotations_all_sums$pct_mut_noncoding[1])/gene_annotations_all_sums$pct_mut_noncoding[2]
plot_bars_rice(gene_annotations_all_sums, yvar="pct_mut_noncoding", xlab="PnPs",ylab="Mutations/bp",ggtitle = "p = 1e-12")
#N/S
chisq.test(gene_annotations_all_sums[,4:5])
plot_bars_rice(gene_annotations_all_sums, yvar="ns_s", xlab="PnPs",ylab="de novo MA N/S",ggtitle = "p = 0.64")
#H3K4me1
chisq.test(gene_annotations_all_sums[,c(8,3)])
plot_bars_rice(gene_annotations_all_sums, yvar="H3K4me1_pct", xlab="PnPs",ylab="H3K4me1 peaks (% of genes)",ggtitle = "p < 2e-16")

dev.off()


# checking non-genic H3K4me1 peaks
H3K4me1_non_genes<-foverlaps(H3K4me1, gene_annotations_basic,type="any")
H3K4me1_non_genes<-H3K4me1_non_genes[,.(marked=ifelse(sum(!is.na(locus))==0, "unmarked","marked")), by=.(chr, start=i.start, stop=i.stop)]
H3K4me1_non_genes$ID<-1:nrow(H3K4me1_non_genes)
H3K4me1_non_genes$length<-H3K4me1_non_genes$stop-H3K4me1_non_genes$start
table(H3K4me1_non_genes$marked)

out<-plot_peaks(H3K4me1_non_genes[marked=="unmarked"], "H3K4me1, non genic, SBS","Peak region", 10, ns_s)
out_means<-out[[2]][,.(mut=sum(mut), length=sum(length), pct=sum(mut)/sum(length)), by=region=="gene body"]
(out_means$pct[1]-out_means$pct[2])/out_means$pct[1]
chisq.test(out_means[,2:3])

out<-plot_peaks(H3K4me1, "H3K4me1, non genic, SBS","Peak region", 10, ns_s)
out_means<-out[[2]][,.(mut=sum(mut), length=sum(length), pct=sum(mut)/sum(length)), by=region=="gene body"]
(out_means$pct[1]-out_means$pct[2])/out_means$pct[1]
chisq.test(out_means[,2:3])
