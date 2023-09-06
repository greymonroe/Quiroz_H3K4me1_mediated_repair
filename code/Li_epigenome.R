
file="data/Li2019/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bedGraph"
read_Li<-function(file){
  Li_data<-fread(file)
  Li_data$V1<-gsub("chr","Chr", Li_data$V1)
  Li_data$chr<-gsub("chr","Chr", Li_data$V1)
  Li_data$ID<-1:nrow(Li_data)
  Li_data$start<-Li_data$V2
  Li_data$stop<-Li_data$V3
  Li_data$length<-Li_data$stop-Li_data$start+1
  
  setkey(Li_data, chr, start, stop)
  
  return(Li_data)
}

setkey(gene_windows, chr, start, stop)

Li_chr<-Li_data[chr=="Chr1"]
x<-unlist(Li_chr[1])
Li_chr_long<-rbindlist(apply(Li_chr[1:100000], 1, function(x){
  chr=x["chr"]
  pos=as.numeric(x["V2"]):as.numeric(x["V3"])
  depth=as.numeric(x["V4"])
  out<-data.table(chr, pos, end=pos, depth)
 
  return(out)
}))
setkey(Li_chr_long, chr, pos, end)
#Li_chr_long<-rbindlist(Li_chr_long)
Li_overlap<-foverlaps(Li_chr_long, gene_windows[chr=="Chr1"])
Li_mean<-Li_overlap[,.(depth=sum(depth, na.rm=T)/mean(length)), by=.(window_ID, region)]

Li_overlap2<-foverlaps(Li_chr, gene_windows[chr=="Chr1"])
Li_mean2<-Li_overlap2[,.(depth=sum(V4*i.length, na.rm=T)/mean(length)), by=.(window_ID, region)]

Li_overlap2[!is.na(i.length)]
Li_mean$depth2<-Li_mean2$depth[match(Li_mean$window_ID, Li_mean2$window_ID)]
Li_mean$diff<-Li_mean$depth-Li_mean$depth2

summary(lm(diff~region, Li_mean))

ggplot(Li_mean, aes(x=depth, y=depth2, col=region))+
  geom_point()

sum(Li_overlap2[window_ID==12]$V4*Li_overlap2[window_ID==12]$i.length)
Li_overlap2[window_ID==12]

file="data/Lu2019/GSE128434_RAW/GSM3674685_ChIP_Rice_7days_leaf_H3K4me1.bw"

bw_dt<-data.table(data.frame(import(con=file)))
bw_dt$score==Li_data$V4

gene_windows$start<-as.integer(gene_windows$start)
gene_windows$stop<-as.integer(gene_windows$stop)

gene_windows$start<-round(sprintf("%.2f", gene_windows$start))
gene_windows$stop<-round(sprintf("%.2f", gene_windows$stop))

write.table(gene_windows[,1:3], "~/Desktop/gene_windows.bed", col.names = F, sep="\t", quote=F, row.names = F)
bed_data <- import.bed("~/Desktop/gene_windows.bed")
asBED(gene_windows)

bedgraph_data <- import.bedGraph(file)

overlapping_counts <- findOverlaps(bed_data, bedgraph_data, type = "any")
window_depth <- tapply(bedgraph_data$score[subjectHits(overlapping_counts)], queryHits(overlapping_counts), sum)

length(queryHits(overlapping_counts))
length(window_depth)
nrow(gene_windows)
gene_windows$depth<-0
gene_windows$depth[queryHits(overlapping_counts)]<-window_depth
