# load  --------------------------------------------------------------

source("code/libraries_functions.R")
library(polymorphology2)
source("code/load_rice_data.R")

N<-nrow(ns_s_filt[WithinCDS=="Yes"])
reps=nrow(ns_s_filt[WithinCDS=="Yes"])

cds<-read.fasta("data/all.fa.cds")
genome<-read.fasta("data/Osativa_204_softmasked.fa.gz")
gff <- read.GFF("data/all.clean.gff")
CDS_regions<-gff[TYPE=="CDS"]
c<-names(genome)[1]
non_CDS_genome<-lapply(names(genome), function(c){
  chrom<-genome[[c]]
  CDS_regions_chrom<-CDS_regions[CHROM==c]
  CDS_regions_chrom_pos<-unlist(apply(CDS_regions_chrom, 1, function(r){
    as.numeric(r["START"]):as.numeric(r["STOP"])
  }))
  chrom_nonCDS<-chrom[-CDS_regions_chrom_pos]
})

tests_nonCDS<-nonsyn_syn_exp3(muts = ns_s_filt[WithinCDS=="No"], genome, CDS = cds, reps = reps)

tests_all<-nonsyn_syn_exp3(muts = ns_s_filt[], non_CDS_genome, CDS = cds, reps = reps)

tests_CDS<-nonsyn_syn_exp3(muts = ns_s_filt[WithinCDS=="Yes"], cds, CDS = cds, reps = reps)

bootN=1000
ratio_nonCDS<-ratio_calc(tests_nonCDS, bootN, N, 0)
ratio_all<-ratio_calc(tests_nonCDS, bootN, N, 0)
ratio_CDS<-ratio_calc(tests_nonCDS, bootN, N, 0)

ratio_nonCDS_s<-ratio_calc(tests_nonCDS, bootN, N, .2)
ratio_all_s<-ratio_calc(tests_nonCDS, bootN, N, .2)
ratio_CDS_s<-ratio_calc(tests_nonCDS, bootN, N, .2)

ratio_nonCDS_s2<-ratio_calc(tests_nonCDS, bootN, N, .2/.4)
ratio_all_s2<-ratio_calc(tests_nonCDS, bootN, N, .2/.4)
ratio_CDS_s2<-ratio_calc(tests_nonCDS, bootN, N, .2/.4)


ratios<-rbindlist(list(
  data.table(ratio=ratio_nonCDS, mut="non-CDS",sel="(Neutral)"),
  data.table(ratio=ratio_all, mut="All",sel="(Neutral)"),
  data.table(ratio=ratio_CDS, mut="CDS",sel="(Neutral)"),
  data.table(ratio=ratio_nonCDS_s, mut="non-CDS",sel="(-20% Selection CDS)"),
  data.table(ratio=ratio_all_s, mut="All",sel="(-20% Selection CDS)"),
  data.table(ratio=ratio_CDS_s, mut="CDS",sel="(-20% Selection CDS)"),
  data.table(ratio=ratio_nonCDS_s2, mut="non-CDS",sel="(-70% Selection CDS)"),
  data.table(ratio=ratio_all_s2, mut="All",sel="(-70% Selection CDS)"),
  data.table(ratio=ratio_CDS_s2,mut="CDS",sel="(-70% Selection CDS)")
  
))


pdf("figures/selection_simulation_NsS2.pdf", width=2.25, height=1.6)
sataka<-(c(
  .717/(1-.717),
  .712/(1-.712),
  .716/(1-.716),
  .723/(1-.723)))
PnPs<-(sum(genes_filt$Pn)+sum(genes_filt$stop_gained))/sum(genes_filt$Ps)

ggplot(ratios, aes(x=(ratio), fill=mut))+
  geom_histogram(col=NA, alpha=0.5, bins=100)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Non-Synonymous/Synonymous")+
  scale_fill_manual(values=c("gray20","gray40","gray60","gray80"), name="Sim. Mutation\nSpectra")+
  annotate(x=(2.56), y=500, geom="point", shape=25, fill="red")+
  annotate(x=PnPs, y=550, geom="point", shape=25, fill="purple")+
  annotate(x=sataka, y=0, geom="point", shape=24, fill="green4")+
  scale_y_continuous(name="N", limits=c(-50, 1000))+
  theme(legend.key.size = unit(0.25, "cm"), 
        legend.position = c(0.30, 0.30),
        legend.background = element_blank(),
        panel.grid = element_blank(),panel.background = element_blank(), plot.background = element_blank())

dev.off()


percentiles <- ratios[, .(mean=mean(ratio),
                          lower_95CI = quantile(ratio, probs = 0.025),
                          upper_95CI = quantile(ratio, probs = 0.975)), 
                      by = .(Selection=sel, Spectrum=mut)]
fwrite(percentiles, "tables/ratios_ns_s.txt")
library(Hmisc)
# Define a custom function to calculate the mean and confidence interval


ns<-nrow(ns_s_filt[MutationType=="non-synonymous"])
s<-nrow(ns_s_filt[MutationType=="synonymous"])

percentiles$chip<-sapply(percentiles$mean, function(m){
  test<-chisq.test(c(ns,s), p=prop.table(c(m, 1)))
  test$p.value
})

test<-chisq.test(c(ns,s), p=prop.table(c(PnPs, 1)))


L<-length(genome)
L=2
bp_freq<-rbindlist(lapply(1:L, function(i){
  bp_freq<-data.table(table(REF=genome[[i]]))[REF!="n"]
  bp_freq$REF<-toupper(bp_freq$REF)
  return(bp_freq)
}))
bp_freq<-bp_freq[,.(N=sum(N)), by="REF"]


