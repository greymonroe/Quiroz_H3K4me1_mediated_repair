source("code/libraries_functions.R")

mutations<-fread("data/Strelka2_mutations.csv")

somatic<-mutations[type=="somatic"]
msh6<-somatic[trt=="MSH6"]
wt<-somatic[trt=="WT"]

pdf("figures/SBS_somatic.pdf", width=6, height=1.5)
contexts_plot(msh6$context_SBS)[[2]]
contexts_plot(wt$context_SBS)[[2]]
dev.off()

pdf("figures/SBS_mutations.pdf", width=6, height=1.5)
msh6<-mutations[trt=="MSH6"]
wt<-mutations[trt=="WT"]
contexts_plot(msh6$context_SBS)[[2]]
contexts_plot(wt$context_SBS)[[2]]
dev.off()

pdf("figures/SBS_mutations_simple.pdf", width=1, height=2)
msh6<-mutations[trt=="MSH6"]
wt<-mutations[trt=="WT"]
contexts_plot(msh6$context_SBS, full = F)[[2]]+geom_bar(stat="identity", fill=NA, col="black", width=0.5)
contexts_plot(wt$context_SBS, full = F)[[2]]+geom_bar(stat="identity", fill=NA, col="black", width=0.5)
dev.off()




pdf("figures/SBS_somatic_simple.pdf", width=1, height=1.5)
contexts_plot(msh6$context_SBS, full = F)[[2]]
contexts_plot(wt$context_SBS, full = F)[[2]]
dev.off()