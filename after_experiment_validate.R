library(gdata)
setwd("/mnt/data5/BGI/UCB/tangchao/for_vilidate/validated")

col.a <- "#800000" # col.annotated
col.u <- "#008080" # col.unannotated
col.h <- "#CD853F" # col.HFnovel
col.l <- "#EEE8AA" # col.LFnovel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron

validated_sj <- read.xls("./UCB experiment validated.xlsx", sheet=1, row.names=1, na.strings="NA")
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
sj <- sj[sj$exon_Shannon_entropy >= 1.8 & 
           !is.na(sj$exon_Relative_entropy) &
           sj$ReadSum >= 10 & 
           sj$overhang >= 10, ]
dim(sj)
# [1] 457282     37

sum(sj$novel=="Y" & sj$CellSum>1)
# [1] 54183

sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj[sj$novel=="Y" & sj$CellSum>1, ]$classification <- "HFNovel"
sj[sj$novel=="Y" & sj$CellSum==1, ]$classification <- "LFNovel"

sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","HFNovel","LFNovel"))
table(sj$classification)
#  Annotated Unannotated     HFNovel     LFNovel
#     155659       72684       54183      174756

which(!validated_sj$Splice.junction %in% sj$sj)
[1] 13 57
validated_sj[which(!validated_sj$Splice.junction %in% sj$sj),]

validated_sj <- validated_sj[validated_sj$Splice.junction %in% sj$sj, ]
validated_sj[,8] <- as.character(validated_sj[,8])
validated_sj[,9] <- as.character(validated_sj[,9])
validated_sj <- validated_sj[nchar(validated_sj[,8])!=0, ]

table(sj[sj$sj %in% validated_sj$Splice.junction, "classification"])
#  Annotated Unannotated     HFNovel     LFNovel
#          6          22          13          19
table(validated_sj$Category)




#### 1. MaxEntScan for Maximum Entropy of splicing sites ========================================================================================


library(MASS)
library(ggplot2)
library(viridis)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

annotated <- sj[sj$classification == "Annotated", c("donor_MaxEnt","acceptor_MaxEnt")]
row.names(annotated) <- sj[sj$classification == "Annotated", "sj"]
colnames(annotated) <- c("Donor","Acceptor")
annotated$density <- get_density(annotated$Donor, annotated$Acceptor)

unannotated <- sj[sj$classification == "Unannotated", c("donor_MaxEnt","acceptor_MaxEnt")]
row.names(unannotated) <- sj[sj$classification == "Unannotated", "sj"]
colnames(unannotated) <- c("Donor","Acceptor")
unannotated$density <- get_density(unannotated$Donor, unannotated$Acceptor)

HFnovel <- sj[sj$classification == "HFNovel", c("donor_MaxEnt","acceptor_MaxEnt")]
row.names(HFnovel) <- sj[sj$classification == "HFNovel", "sj"]
colnames(HFnovel) <- c("Donor","Acceptor")
HFnovel$density <- get_density(HFnovel$Donor, HFnovel$Acceptor)

LFnovel <- sj[sj$classification == "LFNovel", c("donor_MaxEnt","acceptor_MaxEnt")]
row.names(LFnovel) <- sj[sj$classification == "LFNovel", "sj"]
colnames(LFnovel) <- c("Donor","Acceptor")
LFnovel$density <- get_density(LFnovel$Donor, LFnovel$Acceptor)


p1 <- ggplot(annotated, aes(x = Donor, y = Acceptor, color = density))+
  geom_point(size = .6)+
  scale_color_viridis()+
  theme(legend.position = "none")+
  ggtitle("Annotated")+
  geom_point(data = annotated[row.names(annotated) %in% validated_sj$Splice.junction,], aes(x = Donor, y = Acceptor), color = "red", size = 1.2)

p2 <- ggplot(unannotated, aes(x = Donor, y = Acceptor,color = density))+
  geom_point(size = .6)+
  scale_color_viridis()+
  theme(legend.position = "none")+
  ggtitle("Unannotated")+
  geom_point(data = unannotated[row.names(unannotated) %in% validated_sj$Splice.junction,], aes(x = Donor, y = Acceptor), color = "red", size = 1.2)

p3 <- ggplot(HFnovel, aes(x = Donor, y = Acceptor,color = density))+
  geom_point(size = .6)+
  scale_color_viridis()+
  theme(legend.position = "none")+
  ggtitle("HF novel")+
  geom_point(data = HFnovel[row.names(HFnovel) %in% validated_sj$Splice.junction,], aes(x = Donor, y = Acceptor), color = "red", size = 1.2)

p4 <- ggplot(LFnovel, aes(x = Donor, y = Acceptor,color = density))+
  geom_point(size = .6)+
  scale_color_viridis()+
  theme(legend.position = "none")+
  ggtitle("LF novel")+
  geom_point(data = LFnovel[row.names(LFnovel) %in% validated_sj$Splice.junction,], aes(x = Donor, y = Acceptor), color = "red", size = 1.2)


library(cowplot)
pdf("/mnt/data5/BGI/UCB/tangchao/for_vilidate/validated/Maximum Entropy plot of validated sj.pdf")
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()



## Boxplot of novel junc distance

d1 <- log10((as.numeric(sj[sj$classification == "Annotated", ]$end)+1) - (as.numeric(sj[sj$classification == "Annotated", ]$start)-1))
d2 <- log10((as.numeric(sj[sj$classification == "Unannotated", ]$end)+1) - (as.numeric(sj[sj$classification == "Unannotated", ]$start)-1))
d3 <- log10((as.numeric(sj[sj$classification == "HFNovel", ]$end)+1) - (as.numeric(sj[sj$classification == "HFNovel", ]$start)-1))
d4 <- log10((as.numeric(sj[sj$classification == "LFNovel", ]$end)+1) - (as.numeric(sj[sj$classification == "LFNovel", ]$start)-1))

boxplot(d1,d2,d3,d4, names = F, xaxt = "n", ylab = "log10(bp)", col = c(col.a, col.u, col.h, col.l))
#dev.off()

names(d1) <- sj[sj$classification == "Annotated", "sj"]
names(d2) <- sj[sj$classification == "Unannotated", "sj"]
names(d3) <- sj[sj$classification == "HFNovel", "sj"]
names(d4) <- sj[sj$classification == "LFNovel", "sj"]

d1 <- data.frame(Distance = d1, Category = "Annotated", row.names = names(d1))
d2 <- data.frame(Distance = d2, Category = "Unannotated", row.names = names(d2))
d3 <- data.frame(Distance = d3, Category = "HFNovel", row.names = names(d3))
d4 <- data.frame(Distance = d4, Category = "LFNovel", row.names = names(d4))
dis_tab <- rbind(d1,d2,d3,d4)
dis_tab$Category <- as.factor(dis_tab$Category)

pdf("/mnt/data5/BGI/UCB/tangchao/for_vilidate/validated/Boxplot of validated junctions distance.pdf")

ggplot(data = dis_tab, aes(x = Category, y = Distance, fill = Category))+
	geom_boxplot()+
	scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  	theme(legend.position = "none")+
  	geom_point(data = dis_tab[row.names(dis_tab) %in% validated_sj$Splice.junction,], aes(x = Category, y = Distance), color = "red", size = 2)

dev.off()



#### GC Content


ggplot(sj, aes(x = classification, y = donor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  geom_point(data = sj[sj$sj %in% validated_sj$Splice.junction, ], aes(x = classification, y = donor_exon_GC), color = "red", size = 2)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")




p1 <- ggplot(sj, aes(x = classification, y = donor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  geom_point(data = sj[sj$sj %in% validated_sj$Splice.junction, ], aes(x = classification, y = donor_exon_GC), color = "red", size = 1)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")


p2 <- ggplot(sj, aes(x = classification, y = acceptor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  geom_point(data = sj[sj$sj %in% validated_sj$Splice.junction, ], aes(x = classification, y = acceptor_exon_GC), color = "red", size = 1)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Acceptor exon")+
  labs(y = "GC Content")



p3 <- ggplot(sj, aes(x = classification, y = donor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  geom_point(data = sj[sj$sj %in% validated_sj$Splice.junction, ], aes(x = classification, y = donor_intron_GC), color = "red", size = 1)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor intron")+
  labs(y = "GC Content")


p4 <- ggplot(sj, aes(x = classification, y = acceptor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
  geom_point(data = sj[sj$sj %in% validated_sj$Splice.junction, ], aes(x = classification, y = acceptor_intron_GC), color = "red", size = 1)+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Acceptor intron")+
  labs(y = "GC Content")

pdf("/mnt/data5/BGI/UCB/tangchao/for_vilidate/validated/Vioplot of validated junctions GC content.pdf")
library(cowplot)
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()





