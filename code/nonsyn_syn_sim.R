
library(polymorphology2)

muts <- fread("data/Mutations_Final_StopCodonAnnotated_IDAdded.csv")
table(muts$MutationType)
muts$REF<-muts$Old
muts$ALT<-muts$New
table(muts$ALT, muts$REF)
CDS <- read.fasta("data/all.cds")
result <- nonsyn_syn_exp(muts[MutationType %in% c("non-synonymous", "synonymous")], CDS, 10, 5370+2155)

mean(result$ratio)
hist((result$ratio))

sim<-rbindlist(result$simmuts)

chisq.test(table(sim$mut), p=prop.table(table(muts[MutationType %in% c("non-synonymous", "synonymous")]$`Single base substitution`)))

prop.table(table(sim$mut))-prop.table(table(muts[MutationType %in% c("non-synonymous", "synonymous")]$`Single base substitution`))

(5370+x)/2155=2.692239

x=(2.566965*2155)-5370
x

chisq.test(c(5370,2155), p=prop.table(c(2.566965, 1)))
5370/2155

sum(result$ratio>2.491879)

ggplot()+geom_histogram(aes(x=result$ratio))+
  geom_vline(xintercept = 5370/2155, col="red")
