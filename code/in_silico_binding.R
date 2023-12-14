
files <- c("arat_msh6", "arat_pds5c", "rice_msh6", "rice_pds5c")

load_scorefile <- function(file) {
  base_path <- paste0("data/Rosetta data/", file)
  
  scorefile.none <- fread(paste0(base_path, "/pdbs-none.list_distance_average.txt"), skip=0, header=F)
  scorefile.mono <- fread(paste0(base_path, "/pdbs-mono.list_distance_average.txt"), skip=0, header=F)
  scorefile.di   <- fread(paste0(base_path, "/pdbs-di.list_distance_average.txt"), skip=0, header=F)
  scorefile.tri  <- fread(paste0(base_path, "/pdbs-tri.list_distance_average.txt"), skip=0, header=F)
  
  scorefile.none$Lys <- "H3K4me0"
  scorefile.mono$Lys <- "H3K4me1"
  scorefile.di$Lys   <- "H3K4me2"
  scorefile.tri$Lys  <- "H3K4me3"
  
  scorefile <- rbind(scorefile.none, scorefile.mono, scorefile.di, scorefile.tri)
  scorefile$file <- file
  
  return(scorefile)
}

# Applying the function to each file
results_list <- lapply(files, load_scorefile)

# Combining results into one data table
combined_results <- rbindlist(results_list)

combined_results

means_in_silico<-combined_results[,.(meandist=mean(V2), dist5.5=sum(V2<5.5)/.N, dist6=sum(V2<6)/.N), by=.(Lys, file)]

means_in_silico<-combined_results[V2<5.5,.(H3K4me1=sum(Lys=="H3K4me1")/.N, 
                                     H3K4me2=sum(Lys=="H3K4me2")/.N, 
                                     H3K4me0=sum(Lys=="H3K4me0")/.N, 
                                     H3K4me3=sum(Lys=="H3K4me3")/.N
                                     ), by=.(file)]
means_in_silico_melt<-melt(means_in_silico, id.vars = c("file"))
pdf("figures/in_silico_tudors.pdf", width=3, height=1.75)
means_in_silico_melt$variable<-factor(means_in_silico_melt$variable, levels=c("H3K4me0", "H3K4me1", "H3K4me2", "H3K4me3"))
means_in_silico_melt$file<-factor(means_in_silico_melt$file, levels=c("arat_pds5c", "arat_msh6", "rice_pds5c", "rice_msh6"))

ggplot(means_in_silico_melt, aes(x=variable, y=value*100, fill=variable))+
  geom_bar(stat="identity", col="black", width=0.5)+
  facet_grid(~file)+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c(H3K4me0="gray",H3K4me1="red", H3K4me2="green4",H3K4me3="blue"), guide="none")+
  scale_y_continuous(name="% of models within the aromatic cage")+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()  


pdf("figures/in_silico_tudors_pds_comp.pdf", width=1.5, height=1.5)
tmp<-means_in_silico_melt[file=="arat_pds5c" & variable!="H3K4me0"]
tmp$y<-c(1.47, 30.9, 66.2)
ggplot(tmp, aes( x=y, y=value*100, fill=variable))+
  geom_point(shape=21)+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c(H3K4me1="red", H3K4me2="green4",H3K4me3="blue"), guide="none")+
  scale_y_continuous(name="% of models within\n the aromatic cage")+
  scale_x_continuous("ITC Kd (uM)")
dev.off()  




