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

files<-list.files("data/strelka_out/", full.names = F)
normals<-unique(sapply(files, function(x) paste(unlist(strsplit(x,split = "_"))[1:2], collapse="_")))
tumor<-unique(sapply(files, function(x) paste(unlist(strsplit(x,split = "_"))[3:4], collapse="_")))

results<-rbindlist(lapply(normals, function(norm){
  cat(norm)
  mutations<-rbindlist(lapply(tumor, function(f){
    if(norm!=f){
      dat<-read.vcfR(paste0("data/strelka_out/",norm,"_",f,"/results/variants/somatic.snvs.vcf.gz"), verbose = F)
      calls<-cbind(data.table(dat@fix),data.table(dat@gt))
      calls$file<-f
      calls$norm<-norm
      calls$EVS<-as.numeric(gsub(".+SomaticEVS=","",calls$INFO))
      calls$POS<-as.numeric(calls$POS)
      calls$start<-calls$POS
      calls$stop<-calls$POS
      calls$CHROM<-as.character(calls$CHROM)
      calls$POSITION<-calls$POS
      return(calls)
    }
  }))
  
  mutations$unique<-paste(mutations$CHROM, mutations$POS, sep = "_")
  mutations$ID<-1:nrow(mutations)
  tab<-data.table(table(mutations$unique))
  mutations$Mut_N<-tab$N[match(mutations$unique, tab$V1)]
  setkey(mutations, CHROM, start, stop)
  mutations$genic<-features_overlap_mutation(features = genes, mutations = mutations)
  return(mutations)
}))

results$norm_geno<-substr(results$norm,1,1)
fwrite(results, "data/strelka2_results.csv")

setkey(results, CHROM, start, stop)
results$unique2<-paste(results$unique, results$file, sep = "_")
tab<-data.table(table(results$unique2))
results$results_N<-tab$N[match(results$unique2, tab$V1)]
PASS<-results[Mut_N==1 & CHROM %in% 1:5 & FILTER=="PASS"]

PASS$unique3<-paste(PASS$unique, PASS$file, sep = "_")
tab<-data.table(table(PASS$unique3))
PASS$N<-tab$N[match(PASS$unique3, tab$V1)]
PASS$geno<-substr(PASS$file,1,1)
PASS$trt<-ifelse(PASS$geno=="2", "WT","MSH6") # WT is wt/wt + wt/msh6 based on analyses of heterozygous SALK_037557
PASS$trt[PASS$file %in% c("1G_S6","5J_S17","5B_S15")]<-"WT"

PASS$tumor_ref_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["REF"])],split=","))[1]))
PASS$tumor_alt_depth<-as.numeric(apply(PASS, 1, function(x) unlist(strsplit(unlist(strsplit(x["TUMOR"], split = ":"))[4+which(bp==x["ALT"])],split=","))[1]))
PASS$tumor_depth<-PASS$tumor_ref_depth+PASS$tumor_alt_depth
PASS$depth_pct<-PASS$tumor_alt_depth/PASS$tumor_depth

PASS2<-PASS[N==16]

unique4<-unique(PASS2[,c("CHROM","POS","REF","ALT","unique","unique3", "tumor_alt_depth","tumor_depth","trt","genic","file","Mut_N","geno","POSITION","depth_pct","results_N"), with=F])
unique4$dup<-duplicated(unique4$unique3)

unique4$chi<-apply(unique4,1,  function(x){
  chisq.test(c(as.numeric(x["tumor_alt_depth"]),as.numeric(x["tumor_depth"])-as.numeric(x["tumor_alt_depth"])))$p.value
})

unique4$dist<-sapply(1:nrow(unique4), function(i){
  x<-unlist(unique4[i])
  chr<-x["CHROM"]
  POS<-as.numeric(x["POS"])
  f<-x["file"]
  dist<-PASS[CHROM==chr & file==f]$POS-POS
  dist<-min(abs(dist[dist!=0]))
  return(dist)
})

unique4$type<-ifelse(unique4$depth_pct>0.5 & unique4$chi<0.05,"homozygous","somatic")
unique4$type[unique4$chi>(0.001)]<-"heterozygous"
PASS3<-unique4[dup==F & dist>1000 & tumor_depth>150  & tumor_alt_depth>5]

PASS3$context_SBS<-contexts(vars = PASS3 ,fasta = genome)
PASS3$mut <- paste(substr(PASS3$context, 
                          2, 2), substr(PASS3$context, 6, 6), sep = ">")
PASS3$transversion<-!PASS3$mut %in% c("C>T","T>C")
PASS3$CtoA<-PASS3$mut %in% c("C>A")

PASS3$ID<-1:nrow(PASS3)
PASS3$start<-PASS3$POS
PASS3$stop<-PASS3$POS
setkey(PASS3, CHROM, start, stop)
PASS3$trt<-factor(PASS3$trt, levels=c("WT","MSH6"))
PASS3$TE<-features_overlap_mutation(TE, PASS3)
PASS3$IG<-!PASS3$genic & !PASS3$TE
PASS3$loc<-ifelse(PASS3$genic, "genic","intergenic")
PASS3$loc[PASS3$TE]<-"TE"
fwrite(PASS3, "data/MSH6_KO_mutations.csv")