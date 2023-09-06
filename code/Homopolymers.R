
# make homopolymers
homopolymer_rice<-make_homopolymer(genome, size=3)
save(homopolymer_rice, file = "data/homopolymer_rice.Rda")
homopolymers=homopolymer_rice

homopolymer_context<-homopolymer_context(snps[1:10],genome, homopolymer_rice, 3)
vars<-snps
dist=1

snps$context10<-long_context(snps, 10)
snps$homopolymer_neighbor<-homopolymer_var_annotate(snps, homopolymer_rice, 3, 1)

prop.table(table(snps$homopolymer_neighbor))

table(snps[homopolymer_neighbor==T]$ALT)
