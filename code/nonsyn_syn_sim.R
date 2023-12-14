
library(polymorphology2)

cds<-read.fasta("data/all.fa.cds")
mutations<-fread("data/Mutations_Final_StopCodonAnnotated_IDAdded.csv")
mutations$REF<-mutations$Old
mutations$ALT<-mutations$New

ns_s_filt$REF<-ns_s_filt$Old
ns_s_filt$ALT<-ns_s_filt$New
N<-nrow(ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes"])

tests<-nonsyn_syn_exp(muts = ns_s_filt[MAPPABILITY==1], CDS = cds, reps = 100, n = N)
tests_new<-nonsyn_syn_exp_new(muts = ns_s_filt[MAPPABILITY==1], genome, CDS = cds, reps = 100, n = N)
tests_new2<-nonsyn_syn_exp_new(muts = ns_s_filt[MAPPABILITY==1], genome, CDS = cds, reps = 10, n = N)

t.test(tests$ratio, tests_new2$ratio)

tests_CDS<-nonsyn_syn_exp(muts = ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes"], CDS = cds, reps = 100, n = N)

tests_s<-nonsyn_syn_exp(muts = ns_s_filt[MAPPABILITY==1 & MutationType=="synonymous"], CDS = cds, reps = 100, n = N)

tests_nonCDS<-nonsyn_syn_exp(muts = ns_s_filt[MAPPABILITY==1 & WithinCDS=="No"], CDS = cds, reps = 100, n = N)

mean(tests_nonCDS$ratio)
mean(tests$ratio)
hist(tests_new$ratio, breaks=20)

mean(tests_s$ratio)
mean(tests_CDS$ratio)


ratio_s_selection<<-sapply(tests_nonCDS$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])-N*.20
  ratio<-(nonsyn)/syn
})

mean(ratio_s_selection)

table(ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes"]$REF, ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes"]$ALT)

table(REF=ns_s_filt[MAPPABILITY==1 & WithinCDS=="No"]$REF, ns_s_filt[MAPPABILITY==1 & WithinCDS=="No"]$ALT)

ns_s_filt$context<-tricontexts(ns_s_filt, genome)

plot_tricontexts(ns_s_filt[MAPPABILITY==1 & WithinCDS=="No" & genic==T & Genotype=="Homozygous"]$context, full = F)
plot_tricontexts(ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes" & genic==T& Genotype=="Homozygous"]$context, full = F)

plot_tricontexts(ns_s_filt[MutationType=="non-synonymous"]$context)
plot_tricontexts(ns_s_filt[MAPPABILITY==1 & WithinCDS=="Yes"& genic==T]$context)


gff_nmer<-Nmerfrequency()


ns_s_filt$ID<-1:nrow(ns_s_filt)
ns_s_filt_genic<-features_in_sites(gene_annotations_all, ns_s_filt)
ns_s_filt$genic<-ns_s_filt_genic$overlaps

table(ns_s_filt$genic, ns_s_filt$MutationType)

ratio_s<-sapply(tests_s$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])
  ratio<-(nonsyn)/syn
})

ratio_s_selection<<-sapply(tests_s$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])-nrow(mutations[WithinCDS=="Yes"])*.29
  ratio<-(nonsyn)/syn
})

ratio_s_selection2<<-sapply(tests_s$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])-nrow(mutations[WithinCDS=="Yes"])*.29/.41
  ratio<-(nonsyn)/syn
})

ratio_CDS<-sapply(tests_CDS$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])
  ratio<-(nonsyn)/syn
})

ratio_CDS_selection<-sapply(tests_CDS$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])-nrow(mutations[WithinCDS=="Yes"])*.29
  ratio<-(nonsyn)/syn
})

ratio_CDS_selection2<-sapply(tests_CDS$simmuts, function(i){
  stop<-nrow(i[mut_AA=="*"])
  syn<-nrow(i[non_syn==F])
  nonsyn<-nrow(i[non_syn==T])-nrow(mutations[WithinCDS=="Yes"])*.29/.41
  ratio<-(nonsyn)/syn
})

ratios<-rbindlist(list(
  data.table(ratio=ratio_s, mut="Syn",sel="(Neutral)"),
  data.table(ratio=ratio_s_selection, mut="Syn",sel="(-29% Selection)"),
  data.table(ratio=ratio_s_selection2, mut="Syn",sel="(-29% Selection all)"),
  data.table(ratio=ratio_CDS, mut="CDS",sel="(Neutral)"),
  data.table(ratio=ratio_CDS_selection, mut="CDS",sel="(-29% Selection)"),
  data.table(ratio=ratio_CDS_selection2, mut="CDS",sel="(-29% Selection all)")
  
))

ratios[,.(mean(ratio)), by=pasted]

fwrite(ratios, "~/Dropbox/Research/rice mutation paper/tables/Ns_s_sim.csv")
ratios<-fread("tables/Ns_s_sim.csv")
ratios$pasted<-factor(paste(ratios$sel, ratios$mut, sep="\n"))

pdf("figures/selection_simulation_NsS.pdf", width=3.5, height=1.55)
ggplot(ratios, aes(x=ratio, fill=pasted))+
  geom_histogram(col="black", alpha=0.5, bins=100)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Non-Synonymous (incl. stop) / Synonymous")+
  scale_fill_manual(values=rev(c("lightblue","blue","gold","goldenrod")), name="Simulation")+
  annotate(x=2.49, y=40, geom="point", shape=25, fill="red")+
  scale_y_continuous(name="")+
  theme(legend.key.size = unit(0.25, "cm"), 
        legend.position = "right",
        legend.background = element_blank(),
        panel.grid = element_blank(),panel.background = element_blank(), plot.background = element_blank())
dev.off()

table(ns_s_filt$MutationType)
table(ns_s_filt$`StopCodon?`)
(4852-290)/1862

PnPs<-(sum(genes_filt$Pn)+sum(genes_filt$stop_gained))/sum(genes_filt$Ps)

pdf("figures/selection_simulation_NsS2.pdf", width=2.25, height=1.5)
sataka<-c(
  .717/(1-.717),
  .712/(1-.712),
  .716/(1-.716),
  .723/(1-.723))

ggplot(ratios, aes(x=ratio, fill=mut))+
  geom_histogram(col=NA, alpha=0.5, bins=100)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Non-Synonymous (incl. stop) / Synonymous")+
  scale_fill_manual(values=c(Syn="lightblue",CDS="blue"), name="Spectra")+
  annotate(x=2.61, y=225, geom="point", shape=25, fill="red")+
  annotate(x=PnPs, y=225, geom="point", shape=25, fill="yellow")+
  annotate(x=sataka, y=0, geom="point", shape=24, fill="green4")+
  scale_y_continuous(name="N", limits=c(-50, 650))+
  theme(legend.key.size = unit(0.25, "cm"), 
        legend.position = c(0.25, 0.25),
        legend.background = element_blank(),
        panel.grid = element_blank(),panel.background = element_blank(), plot.background = element_blank())

dev.off()




nrow(mutations[WithinCDS=="Yes"])/0.71/.41

syn_muts<-nrow(mutations[MutationType=="synonymous"])

ns_missing<-nrow(mutations[WithinCDS=="Yes"])/0.71/.41-nrow(mutations[WithinCDS=="Yes"])

Nsyn_muts<-nrow(mutations[MutationType=="non-synonymous"])
(Nsyn_muts+ns_missing)/syn_muts

Nsyn_muts/nrow(mutations[WithinCDS=="Yes"])




