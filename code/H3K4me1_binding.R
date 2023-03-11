library(data.table)
dat<-fread("data/Rosetta data/Book1.csv")[dist==6][,1:6]
dat<-melt(dat)
dat$Protein<-factor(dat$Protein, levels=c("PDS5C","MSH6"))

pdf("figures/Model_Tudor_H3K4.pdf", width=2.2, height=1.2)
ggplot(dat, aes(x=substr(variable, 5,7), y=value, fill=variable))+
  geom_bar(stat="identity", col="black")+
  facet_wrap(~Species+Protein, nrow=1)+
  theme_classic(base_size = 6)+
  theme(axis.text = element_text(angle=45, hjust=1), 
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  scale_x_discrete(name="H3K4")+
  scale_y_continuous(name="% of models\ndistance < 6 Å")+
  scale_fill_discrete(guide="none")
dev.off()

arat_pds5c<-dat[Species=="ARAT" & Protein=="PDS5C"]
arat_pds5c$Kd<-c(NA, 1.47, 30.9, 66.2)

pdf("figures/Kd_pDS5C_Arat.pdf", width=.8, height=1.2)
ggplot(arat_pds5c, aes(x=Kd, y=value, fill=variable))+
  geom_point(pch=21, col="black")+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="% of models\ndistance < 6 Å")+
  scale_x_continuous(name="ITC Kd (u M)", expand = c(.1,.1))+
  scale_fill_discrete(guide="none")
dev.off()


