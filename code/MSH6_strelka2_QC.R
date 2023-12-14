
bp<- c("A","C","G","T")

strelka2_results<-fread("data/strelka2_results.csv")
strelka2_merged_results<-fread("data/strelka2_merged_results.csv")

strelka2_merged_results$normal_ref_depth<-as.numeric(apply(strelka2_merged_results, 1, function(x) unlist(strsplit(unlist(strsplit(x["NORMAL"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
strelka2_merged_results$normal_alt_depth<-as.numeric(apply(strelka2_merged_results, 1, function(x) unlist(strsplit(unlist(strsplit(x["NORMAL"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
strelka2_merged_results$normal_depth<-strelka2_merged_results$normal_ref_depth+strelka2_merged_results$normal_alt_depth
strelka2_merged_results$normal_depth_pct<-strelka2_merged_results$normal_alt_depth/strelka2_merged_results$normal_depth

strelka2_merged_results$tumor_ref_depth<-as.numeric(apply(strelka2_merged_results, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
strelka2_merged_results$tumor_alt_depth<-as.numeric(apply(strelka2_merged_results, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
strelka2_merged_results$tumor_depth<-strelka2_merged_results$tumor_ref_depth+strelka2_merged_results$tumor_alt_depth
strelka2_merged_results$depth_pct<-strelka2_merged_results$tumor_alt_depth/strelka2_merged_results$tumor_depth



strelka2_results$unique_file<-paste(strelka2_results$unique, strelka2_results$file)
strelka2_merged_results$unique_file<-paste(strelka2_merged_results$unique, strelka2_merged_results$sample)
mutations$unique_file<-paste(mutations$unique, mutations$file)


merged_normals<-strelka2_merged_results[,.(normal_depth_pct=mean(normal_depth_pct)), by=unique_file]

EVS<-strelka2_results[,.(EVS=mean(EVS)), by=unique_file]

mutations<-merge(mutations, EVS, by=c("unique_file"))
mutations<-merge(mutations, merged_normals, by=c("unique_file"))



MSH6_genotypes<-fread("data/MSH6_genotypes.csv")

TableS2<-mutations[,.(CHROM, POS, REF, ALT, EVS, tumor_alt_depth, tumor_depth, normal_depth_pct, loc, type, file, MSH6=trt)]
TableS2$genotype=MSH6_genotypes$genotype[match(TableS2$file, MSH6_genotypes$file)]
TableS2$line=MSH6_genotypes$line[match(TableS2$file, MSH6_genotypes$file)]

fwrite(TableS2,"tables/S1_Mutations_WT_MSH6.csv")


strelka2_results$PASS<-strelka2_results$unique_file %in% mutations$unique_file
EVS<-strelka2_results[,.(EVS=mean(EVS)), by=.(unique, PASS)]
EVS$trt<-mutations$trt[match(EVS$unique, mutations$unique)]
EVS$type<-mutations$type[match(EVS$unique, mutations$unique)]
EVS$mut<-mutations$mut[match(EVS$unique, mutations$unique)]
EVS$depth_pct<-mutations$depth_pct[match(EVS$unique, mutations$unique)]

# percent variants passing all filters: 0.4%
table(EVS$PASS)
721/(181736+721)*100

table(EVS$PASS)
EVS$somatic<-EVS$unique %in% mutations$unique
EVS$MSH6<-EVS$unique %in% mutations[trt=="MSH6"]$unique
EVS$group<-ifelse(EVS$PASS==F, "Filtered", ifelse(EVS$MSH6, "msh6", "WT"))
pdf("figures/strelka2_QC_EVS.pdf", width=2.5, height=2)
ggplot(EVS[PASS==F | somatic==T], aes(x=EVS))+
  geom_density(data=EVS[PASS==F ], alpha=0.5, fill="gray")+
  geom_density(data=EVS[somatic==T & MSH6==T ], alpha=0.5, fill="orange")+
  geom_density(data=EVS[somatic==T & MSH6==F ], alpha=0.5, fill="dodgerblue")+
  theme_classic(base_size = 6)+
  scale_x_log10()+
  theme(legend.position = c(0.5,0.5))+
  facet_grid()

EVS$group<-factor(EVS$group, levels=c("Filtered","WT","msh6"))
ggplot(EVS, aes(x=EVS, fill=group))+
  geom_histogram(alpha=0.5, bins=100)+
  scale_fill_manual(values=c("gray","dodgerblue","orange"))+
  theme_classic(base_size = 6)+
  scale_x_log10()+
  facet_grid(group~., scales="free")+
  theme(legend.position = "none")

dev.off()

pdf("figures/strelka2_QC_PASS.pdf", width=1, height=2)
  ggplot(EVS, aes(x=PASS, fill=PASS))+
  geom_bar(stat = "count", col="black")+
  theme_classic(base_size = 6)+
  theme(legend.position = "none")
dev.off()

ggplot(EVS, aes(x=PASS, fill=PASS))+
  geom_bar(stat = "count", col="black")+
  theme_classic(base_size = 6)


ggplot(EVS[PASS==T], aes(x=log(EVS+0.001), fill=trt))+
  geom_density(alpha=0.5)+
  facet_grid(~type)



# no significant different in EVS for WT and MSH6 lines
summary(lm(EVS~trt*type, EVS[PASS==T]))

summary(lm(EVS~trt, EVS[type=="somatic" & PASS==T]))

summary(lm(EVS~mut*trt, EVS[type=="somatic" & PASS==T]))

summary(lm(depth_pct~trt*mut, EVS[type=="somatic" & PASS==T]))

ggplot(EVS[PASS==T & type=="somatic"], aes(x=log(depth_pct), fill=trt))+
  geom_density(alpha=0.5)+
  facet_grid(~mut)

ggplot(EVS[PASS==T & type=="somatic"], aes(x=(EVS), fill=trt))+
  geom_density(alpha=0.5)+
  facet_grid(~mut)


