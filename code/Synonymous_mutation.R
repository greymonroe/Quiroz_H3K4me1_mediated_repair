# Read mutations
library("readxl")
Mutations <- read_excel("data/Mutations-Identified-1504lines.xlsx")
Mutations <- as.data.frame(Mutations)
colnames(Mutations)<-Mutations[3,]
Mutations<-Mutations[-c(1:3),]
Mutations$Old<-gsub(">.*", "", Mutations$`Single base substitution`)
Mutations$New<-gsub(".*>", "", Mutations$`Single base substitution`)
Mutations <- Mutations[complete.cases(Mutations), ] 

# Read gff file
library(ape)
gff <- read.gff("data/all.gff3", GFF3 = TRUE)
CDS<- gff[which(gff$type=="CDS"),]
CDS$ID<-gsub(".*=","",CDS$attributes) #.*= means everything before and including =
CDS$CDSnum<-gsub(";.*", "", CDS$attributes)
CDS$CDSnum<-gsub(".*_", "", CDS$CDSnum)

# CDS Sequence
library(seqinr)
cds<-read.fasta("data/Rice_Mutation_Annotation/all.cds")
ID_list<-names(cds)

# See if the mutation is inside CDS
for (i in 1:nrow(Mutations)){
  chr<-Mutations$Chromosome[i]
  pos<-as.numeric(Mutations$Position[i])
  count<-0
  CDS_chr<-CDS[which(CDS$seqid==chr),]
  for (o in 1:nrow(CDS_chr)){
    if (CDS_chr$start[o]<=pos & CDS_chr$end[o]>=pos){
      count<-count+1
      rownum<-o
    }
  }
  if (count>0){
    Mutations[i,9]<-"Yes"
    ID<-CDS_chr$ID[rownum]
    CDSnum<-as.numeric(CDS_chr$CDSnum[rownum])
    if (CDS_chr$strand[rownum]=="+"){
      Length<-pos-CDS_chr$start[rownum]+1
    }
    if (CDS_chr$strand[rownum]=="-"){
      Length<-CDS_chr$end[rownum]-pos+1
    }
    if (CDSnum>1){
      for (e in 1:(CDSnum-1)){
        Length<-Length+CDS_chr$end[rownum-e]-CDS_chr$start[rownum-e]+1
      }
    }
    series<-which(ID_list==ID)
    WT<-unlist(cds[series])
    if (toupper(WT[Length])==Mutations$Old[i] & CDS_chr$strand[rownum]=="+"){
      Mutant<-WT
      Mutant[Length]<-tolower(Mutations$New[i])
      if (identical(translate(WT),translate(Mutant))){
        Mutations[i,10]<-"synonymous"
      }else{
        Mutations[i,10]<-"non-synonymous"
      }
    }
    if (toupper(WT[Length])==toupper(comp(Mutations$Old[i])) & CDS_chr$strand[rownum]=="-"){
      Mutant<-WT
      Mutant[Length]<-tolower(comp(Mutations$New[i]))
      if (identical(translate(WT),translate(Mutant))){
        Mutations[i,10]<-"synonymous"
      }else{
        Mutations[i,10]<-"non-synonymous"
      }
    }
    if ((toupper(WT[Length])!=Mutations$Old[i] & CDS_chr$strand[rownum]=="+") | (toupper(WT[Length])!=toupper(comp(Mutations$Old[i])) & CDS_chr$strand[rownum]=="-")){
      Mutations[i,10]<-"Mismatch"
    }
  }else{
    Mutations[i,9]<-"No"
  }
  cat(i)
}
colnames(Mutations)[9]<-"WithinCDS"
colnames(Mutations)[10]<-"MutationType"
save(Mutations, file = "Mutations_Final")
write.csv(Mutations,"Mutations_Final.csv", row.names = TRUE)

# Check if there's any mismatch
which(Mutations$V10=="Mismatch")

