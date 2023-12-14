library(readxl)
protein_conc<-rev(c(0,0.033, 0.066, 0.132, 0.264, 0.527, 1.055, 2.109, 4.219, 8.438, 16.875, 33.75, 67.5, 135, 270, 540))/1000000 #µM ->M
WT_data <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 WT#1.xlsx")
WT_data2 <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 WT#2.xlsx")
WT_data3 <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 WT#3.xlsx")

mut_data1 <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 3Amut#1.xlsx")
mut_data2 <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 3Amut#2.xlsx")
mut_data3 <- read_excel("data/FP/MSH6_FP_results_2023Aug_from_Aki/FP_His-MBP-GS-Pre-MSH6 Tudor_H3K4me WTvs3Amut 230807 3Amut#3.xlsx")

raw_of_anisotropy<-seq(from=87,to=267,by=12)-2
WT_an1<-as.data.frame(WT_data[raw_of_anisotropy,3:6]);
WT_an2<-as.data.frame(WT_data2[raw_of_anisotropy,3:6]);
WT_an3<-as.data.frame(WT_data3[raw_of_anisotropy,13:16]);

mut_an1<-as.data.frame(mut_data1[raw_of_anisotropy,8:11]);mut_an1
mut_an2<-as.data.frame(mut_data2[raw_of_anisotropy,8:11]);mut_an2
mut_an3<-as.data.frame(mut_data3[raw_of_anisotropy,18:21]);mut_an3

meanDF<-data.frame(protein_conc)
i<-1;meanDF$WT_m0<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,mean)
i<-2;meanDF$WT_m1<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,mean)
i<-3;meanDF$WT_m2<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,mean)
i<-4;meanDF$WT_m3<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,mean)
i<-1;meanDF$mut_m0<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,mean)
i<-2;meanDF$mut_m1<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,mean)
i<-3;meanDF$mut_m2<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,mean)
i<-4;meanDF$mut_m3<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,mean)

sdDF<-data.frame(protein_conc)
i<-1;sdDF$WT_m0<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,sd)
i<-2;sdDF$WT_m1<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,sd)
i<-3;sdDF$WT_m2<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,sd)
i<-4;sdDF$WT_m3<-apply(data.frame(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),1,sd)
i<-1;sdDF$mut_m0<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,sd)
i<-2;sdDF$mut_m1<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,sd)
i<-3;sdDF$mut_m2<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,sd)
i<-4;sdDF$mut_m3<-apply(data.frame(as.numeric(mut_an1[,i]),as.numeric(mut_an2[,i]),as.numeric(mut_an3[,i])),1,sd)

library(RColorBrewer)
spc<-brewer.pal(11, "Spectral");
cpl<-spc[c(1,3,9,11)]
#mut
plot(rep(meanDF$protein_conc,4),c(meanDF$mut_m0,meanDF$mut_m1,meanDF$mut_m2,meanDF$mut_m3),log='x',
     pch=rep(c(15,16,17,18),each=16),col=rep(cpl,each=16),ylim=c(40,95),
     xlab='receptor concentration (M)',ylab='Anisotropy')
segments(rep(meanDF$protein_conc,4),c(meanDF$mut_m0,meanDF$mut_m1,meanDF$mut_m2,meanDF$mut_m3)+c(sdDF$mut_m0,sdDF$mut_m1,sdDF$mut_m2,sdDF$mut_m3),
         rep(meanDF$protein_conc,4),c(meanDF$mut_m0,meanDF$mut_m1,meanDF$mut_m2,meanDF$mut_m3)-c(sdDF$mut_m0,sdDF$mut_m1,sdDF$mut_m2,sdDF$mut_m3),
         col=rep(cpl,each=16))

# WT
#FP sigmoig model?
# fit <- nls(y ~ a * x^2 + b * x + c, data = data, start = list(a = 1, b = 1, c = 1))
# LT = the total added concentration of ligand
# A = the experimental anisotropy
# Af = the anisotropy for the free ligand
# Ab = the anisotropy for the fully bound ligand
# RT = the total receptor concentration
# ligand = H3 peptide, probably.

# LT is known and A is measured for each RT. The equation can be solved for Kd, Ab, and Af
i<-2
data <- data.frame(
  A<-c(as.numeric(WT_an1[,i]),as.numeric(WT_an2[,i]),as.numeric(WT_an3[,i])),  # Measured A
  Lt=20/10^9, #nM->M 
  RT<-rep(protein_conc,3)  # Varying RT
)

fit <- nls(A ~ Af + (Ab - Af) * ((Lt + Kd + RT) - sqrt((Lt - Kd - RT)^2 - 4 * Lt * RT)) / (2 * Lt),
           data = data,algorithm = "port",
           start = list(Af = 50, Ab = 100, Kd =8.4380e-06))  # Initial guesses for Af, Ab, and Kd
summary(fit)


# MSH6 K4me1 model seem okay

plot(rep(meanDF$protein_conc,4),c(meanDF$WT_m0,meanDF$WT_m1,meanDF$WT_m2,meanDF$WT_m3),log='x',
     pch=rep(c(15,16,17,18),each=16),col=rep(cpl,each=16),ylim=c(40,95),
     xlab='receptor concentration (M)',ylab='Anisotropy')
segments(rep(meanDF$protein_conc,4),c(meanDF$WT_m0,meanDF$WT_m1,meanDF$WT_m2,meanDF$WT_m3)+c(sdDF$WT_m0,sdDF$WT_m1,sdDF$WT_m2,sdDF$WT_m3),
         rep(meanDF$protein_conc,4),c(meanDF$WT_m0,meanDF$WT_m1,meanDF$WT_m2,meanDF$WT_m3)-c(sdDF$WT_m0,sdDF$WT_m1,sdDF$WT_m2,sdDF$WT_m3),
         col=rep(cpl,each=16))

meanDF<-data.table(meanDF)
str(meanDF)

lines(RT, predict(fit, data = RT), col = cpl[2],xlab='',ylab='')
text(3.74e-05,85,"Kd=3.411e-05",col=cpl[2])
abline(v=3.74e-05,lty='dashed',col=cpl[2])




longDF <- melt(meanDF, id.vars = "protein_conc", variable.name = "trt_level", value.name = "value")

# Splitting the 'trt_level' column into 'trt' and 'level'
longDF[, c("trt", "level") := tstrsplit(trt_level, "_", type.convert = TRUE)]

# Dropping the original 'trt_level' column
longDF[, trt_level := NULL]

# Reordering columns if needed
setcolorder(longDF, c("protein_conc", "trt", "level", "value"))

# Melting the data.table
slongDF <- melt(data.table(sdDF), id.vars = "protein_conc", variable.name = "trt_level", value.name = "value")

# Splitting the 'trt_level' column into 'trt' and 'level'
slongDF[, c("trt", "level") := tstrsplit(trt_level, "_", type.convert = TRUE)]

# Dropping the original 'trt_level' column
slongDF[, trt_level := NULL]

# Reordering columns if needed
setcolorder(slongDF, c("protein_conc", "trt", "level", "value"))

longDF$sd<-slongDF$value
str(longDF)


pdf("figures/FP2.pdf", width=3.5, height=2)
longDF$trtf<-factor(longDF$trtf, levels=c("WT","mut"))
p <- ggplot(longDF, aes(x = protein_conc, y = value,  color = level, shape=level)) +
  geom_point() +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0) +
  labs(x = "Concentration (M)", y = "Ansiotropy") +
  scale_x_log10(breaks=c(1e-7, 1e-6, 1e-5, 1e-4), labels = polymorphology2::real_sci_format(c(1e-7, 1e-6, 1e-5, 1e-4))) +  # Use the custom function here
  facet_grid(~trt) +
  theme_bw(base_size = 8) +
  scale_color_manual(values=c(m0="gray",m1="red", m2="green4",m3="blue"))+
  scale_shape_manual(values=c(m0=15,m1=16, m2=17,m3=18))+
  geom_line(data=data.table(value=predict(fit, data = RT), protein_conc=RT, trt=factor("WT", levels=c("mut", "WT")), level="H3K4me1"), col="red")+
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )#+
  # geom_vline(xintercept = 3.74e-05, linetype="dashed", col="red", alpha=0.5)+
  # annotate(geom="text", x=3.74e-05,y=85,label="Kd=3.411e-05",col="red")

p
dev.off()
