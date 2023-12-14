library(polymorphology)
library(vcfR)
library(seqinr)

genome<-read.fasta("data/TAIR10_chr_all.fas.gz")
names(genome)<-gsub("Chr","",names(genome))

#gff<-fread("data/A_thaliana_gff.txt")
genes<-fread("data/A_thal_genes_PDS5_enrich.csv")
#genes<-gff[type=="gene"]
genes$CHROM<-as.character(genes$chr)
genes$ID<-1:nrow(genes)
genes$START<-genes$start
genes$STOP<-genes$stop
setkey(genes, CHROM, start, stop)

TE<-fread("data/TAIR10_Transposable_Elements.txt")
TE$CHROM<-substr(TE$Transposon_Name, 3,3)
TE$ID<-1:nrow(TE)
TE$START<-TE$Transposon_min_Start
TE$STOP<-TE$Transposon_max_End
TE$start<-TE$Transposon_min_Start
TE$stop<-TE$Transposon_max_End
setkey(TE, CHROM, start, stop)

bp<- c("A","C","G","T")

files<-list.files("data/strelka_merged/strelka_out/", full.names = F)

norm<-files[1]
results<-rbindlist(lapply(files, function(norm){
  cat(norm)
if(file.exists(paste0("data/strelka_merged/strelka_out/",norm,"/results/variants/somatic.snvs.vcf.gz"))){
      dat<-read.vcfR(paste0("data/strelka_merged/strelka_out/",norm,"/results/variants/somatic.snvs.vcf.gz"), verbose = F)
      calls<-cbind(data.table(dat@fix),data.table(dat@gt))
      calls$sample<-gsub("_merged","", norm)
      calls$EVS<-as.numeric(gsub(".+SomaticEVS=","",calls$INFO))
      calls$POS<-as.numeric(calls$POS)
      calls$start<-calls$POS
      calls$stop<-calls$POS
      calls$CHROM<-as.character(calls$CHROM)
      calls$POSITION<-calls$POS
      mutations<-calls

      mutations$unique<-paste(mutations$CHROM, mutations$POS, sep = "_")
      
  setkey(mutations, CHROM, start, stop)
  mutations$genic<-features_overlap_mutation(features = genes, mutations = mutations)
  return(mutations)
}
}))

fwrite(results,"data/strel2_merged_results.csv")
PASS<-results[FILTER=="PASS"]
PASS_counts<-data.table(table(PASS$unique))
PASS$count<-PASS_counts$N[match(PASS$unique,PASS_counts$V1 )]
PASS<-PASS[count==1]

PASS$tumor_ref_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
PASS$tumor_alt_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
PASS$tumor_depth<-PASS$tumor_ref_depth+PASS$tumor_alt_depth
PASS$depth_pct<-PASS$tumor_alt_depth/PASS$tumor_depth

PASS$normal_ref_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["NORMAL"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
PASS$normal_alt_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["NORMAL"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
PASS$normal_depth<-PASS$normal_ref_depth+PASS$normal_alt_depth
PASS$normal_depth_pct<-PASS$normal_alt_depth/PASS$normal_depth

PASS$geno<-substr(PASS$sample,1,1)
PASS$trt<-ifelse(PASS$geno=="2", "WT","MSH6") # WT is wt/wt + wt/msh6 based on analyses of heterozygous SALK_037557
PASS$trt[PASS$sample %in% c("1G_S6","5J_S17","5B_S15")]<-"WT"

PASS$context_SBS<-contexts(vars = PASS ,fasta = genome)
PASS$mut <- paste(substr(PASS$context, 
                          2, 2), substr(PASS$context, 6, 6), sep = ">")
PASS$transversion<-!PASS$mut %in% c("C>T","T>C")
PASS$CtoA<-PASS$mut %in% c("C>A")
PASS$chi<-apply(PASS,1,  function(x){
  chisq.test(c(as.numeric(x["tumor_alt_depth"]),as.numeric(x["tumor_depth"])-as.numeric(x["tumor_alt_depth"])))$p.value
})

PASS$ID<-1:nrow(PASS)
PASS$start<-PASS$POS
PASS$stop<-PASS$POS
setkey(PASS, CHROM, start, stop)
PASS$trt<-factor(PASS$trt, levels=c("WT","MSH6"))
PASS$TE<-features_overlap_mutation(TE, PASS)
PASS$IG<-!PASS$genic & !PASS$TE
PASS$loc<-ifelse(PASS$genic, "genic","intergenic")
PASS$loc[PASS$TE]<-"TE"
View(PASS[confirmed==F])

ggplot(PASS ,aes(x=depth_pct, fill=chi<0.05))+
  geom_histogram()+
  facet_wrap(~trt)+
  theme_classic(base_size = 6)+
  scale_x_continuous(name="Alt frequency")+
  scale_y_continuous(name="Mutations")

PASS$type<-ifelse(PASS$depth_pct>0.5 & PASS$chi<0.05,"homozygous","somatic")
PASS$type[PASS$chi>(0.05)]<-"heterozygous"

PASS3<-fread("data/MSH6_KO_mutations.csv");PASS3$CHROM<-as.character(PASS3$CHROM)

PASS$confirmed<-PASS$unique %in% PASS3$unique
confirmed<-PASS[confirmed==T]
confirmed<-merge(PASS3, confirmed, by="unique")
confirmed<-confirmed[type.y=="somatic"]
depths<-data.table(sample=rep(c("Tumor","Normal"), each=nrow(confirmed)), depth=c(confirmed$depth_pct.y, confirmed$normal_depth_pct))
ggplot(depths, aes(x=log10(depth+0.0001), fill=sample))+
  geom_density(alpha=0.5)
ggplot(depths, aes(x=depth+0.0001, fill=sample))+
  geom_density(alpha=0.5)


