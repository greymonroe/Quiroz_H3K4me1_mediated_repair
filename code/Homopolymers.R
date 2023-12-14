
# make homopolymers
homopolymer_rice<-make_homopolymer(genome, size=3)
save(homopolymer_rice, file = "data/homopolymer_rice.Rda")
load("data/homopolymer_rice.Rda")

x<-names(homopolymer_rice)[1]
homopolymers<-rbindlist(lapply(names(homopolymer_rice), function(x){
  tmp<-homopolymer_rice[[x]]
  data.table(tmp[,.(CHROM=x, START=start, STOP=end, LENGTH=length, BP=var)])
}))

homopolymer_context<-homopolymer_context(snps[1:10],genome, homopolymer_rice, 3)
vars<-snps
dist=1

snps$context10<-long_context(snps, 10)
snps$ID<-1:nrow(snps)
snps_homopolymer_neighbor<-homopolymer_var_annotate(snps, homopolymers, 3, 1)
snps$HP<-snps_homopolymer_neighbor$HP[match(snps$ID, snps_homopolymer_neighbor$ID)]
snps$Genotype[snps$Genotype=="Homozygouszygous"]<-"Homozygous"

View(snps[HP==T])
