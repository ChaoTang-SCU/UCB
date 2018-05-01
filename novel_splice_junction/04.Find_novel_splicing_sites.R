#### Find NSS
## Sat Mar 24 2018
## By Tang Chao


#### 0. Basic settings ======================================================================


depar <- par()
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS")

require(VennDiagram)

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/Find_NSS/RData/")


#### 1. Load data ===========================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_tu.RData")
mrna <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/mRNA_intron.sj", stringsAsFactors = F)
est <- read.table("/mnt/data1/projects/UCB/results_ucb/tangchao/GMAP/est_intron.sj", stringsAsFactors = F)
bulk <- read.table("/mnt/data1/projects/UCB/data_ucb/TXT/bulk_AS_events_dictionary.txt",header=T, stringsAsFactors = F)


#### 2. overview of four methods' SJ ========================================================


mRNA_SJ <- as.character(mrna$V1)
EST_SJ <- as.character(est$V1)
bulk_SJ <- as.character(bulk$event_id)
sc_SJ <- as.character(sj$sj)


pdf(paste(fo,"SJ_numble_of_different_splicing_methods.pdf",sep = ""), width = 10, height = 7)
x = barplot(c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ)), names.arg=c("SC","Bulk","EST","mRNA"), main="SJ numble of different splicing methods", ylim = c(0, 900000))
text(x, c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ))+15000, labels = c(length(sc_SJ),length(bulk_SJ),length(EST_SJ),length(mRNA_SJ)))
dev.off()


library(VennDiagram)
venn_plot <- venn.diagram(list(mRNA = mRNA_SJ, EST = EST_SJ, Bulk = bulk_SJ, SC = sc_SJ), fill = rainbow(4), filename = NULL, main.cex = 2, main = "Venn plot of splicing junction in different sequencing methods", height = 6, width = 10)

pdf(paste(fo,"venn_plot_of_sc_bulk_mrna_EST_SJ.pdf",sep = ""), width = 10, height = 7)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


#### 3. find sc specific SJ =================================================================


sum(!as.character(sj$sj) %in% unique(c(mRNA_SJ, EST_SJ, bulk_SJ)))
# [1] 620553
sum(as.character(sj$sj) %in% unique(c(mRNA_SJ, EST_SJ, bulk_SJ)))
# [1] 249498
sum(!as.character(sj$sj) %in% unique(c(mRNA_SJ, EST_SJ, bulk_SJ)))/length(as.character(sj$sj))
# [1] 0.7132375
## That is say 71% of SC sj havn't been validated in other database.

sc_sj <- sj[!sj$sj %in% unique(c(mRNA_SJ, EST_SJ, bulk_SJ)),]

table(sc_sj$annotation)


## annotation of sc sjific SJ
pdf(paste(fo,"annotated_percentage_of_sc_sjific_SJ.pdf",sep = ""), width = 10, height = 7)
x <- barplot(table(sc_sj$annotation)/sum(table(sc_sj$annotation))*100, #col = rainbow(2, alpha=.6, start = .6, end = .5), 
             main = paste("Annotated percentage of SC specific SJ\n(", sum(table(sc_sj$annotation)), ")"),
             names.arg = c("Unannotated", "Annotated"),
             ylim = c(0,118), ylab = "%")
p <- round(table(sc_sj$annotation)/sum(table(sc_sj$annotation))*100,2)
n <- table(sc_sj$annotation)
text(x, table(sc_sj$annotation)/sum(table(sc_sj$annotation))*100, 
     labels = paste(n,"\n",p,sep = ""), pos = 3, cex = .8)
dev.off()


## Barchart of single cell specific SS motif parsed strand
pdf(paste(fo,"Barchart of single cell specific SS motif parsed.pdf", sep = ""), width = 10, height = 7)
par(mfrow = c(2,1))
par(mar = c(1,4,1,1))
x = with(sc_sj, barplot(table(motif[annotation==1])/sum(table(motif[annotation==1]))*100, names.arg = F, ylim = c(0,112), main = "SC specific SS motif", ylab = "%"))
p <- round(with(sc_sj, table(motif[annotation==1]))/sum(with(sc_sj, table(motif[annotation==1])))*100,2)
n <- with(sc_sj, table(motif[annotation==1]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Annotated", bty = "n")

par(mar = c(3,4,1,1))
x = with(sc_sj, barplot(table(motif[annotation==0])/sum(table(motif[annotation==0]))*100, ylim = c(0,100), ylab = "%"))
p <- round(with(sc_sj, table(motif[annotation==0]))/sum(with(sc_sj, table(motif[annotation==0])))*100,2)
n <- with(sc_sj, table(motif[annotation==0]))
text(x = x, y = p, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)
legend("topright", legend = "Unannotated", bty = "n")
par(depar)
dev.off()


#### 4. find novel splicing junction ========================================================


sj$SCSpe <- "N"
sj$novel <- "N"
sj[sj$sj %in% sc_sj$sj, ]$SCSpe <- "Y"
sj[sj$sj %in% sc_sj[sc_sj$annotation==0,]$sj, ]$novel <- "Y"

sj$SCSpe <- as.factor(sj$SCSpe)
sj$novel <- as.factor(sj$novel)

head(sj)
table(sj$SCSpe)
table(sj$novel)
dim(sj)


## Boxplot of cells and reads
pdf(paste(fo,"Boxplot of cells and reads of novel splice junction.pdf", sep = ""))
par(mfrow = c(1,2))
with(sj,boxplot(log10(CellSum[novel=="Y"]),log10(CellSum[novel=="N"]),
                names = c("Novel", "Other"), main = "No. of cell", 
                ylab = "Log10(cells)"))

with(sj,boxplot(log10(ReadSum[novel=="Y"]),log10(ReadSum[novel=="N"]),
                names = c("Novel", "Other"), main = "No. of reads", 
                ylab = "Log10(reads)"))
par(depar)
dev.off()

pdf(paste(fo,"Vioplot of cells and reads of novel splice junction.pdf", sep = ""))
library(vioplot)
with(sj,vioplot(log10(CellSum[novel=="Y"]),log10(CellSum[novel=="N"]), names = c("Novel","Other")))
with(sj,vioplot(log10(ReadSum[novel=="Y"]),log10(ReadSum[novel=="N"]), names = c("Novel","Other")))
head(sj)
library(ggplot2)
ggplot(data = sj, aes(x = novel, y = log10(CellSum)))+
  geom_boxplot()

ggplot(data = sj, aes(x = novel, y = log10(CellSum)))+
  geom_violin()

ggplot(data = sj, aes(x = novel, y = log10(ReadSum)))+
  geom_boxplot()

ggplot(data = sj, aes(x = novel, y = log10(ReadSum)))+
  geom_violin()

dev.off()


## Boxplot of novel junc distance
pdf(paste(fo,"Boxplot of novel junc distance.pdf", sep = ""))
d1 <- log10((as.numeric(sj[sj$novel == "Y", ]$end)+1) - (as.numeric(sj[sj$novel == "Y", ]$start)-1))
d2 <- log10((as.numeric(sj[sj$novel == "N", ]$end)+1) - (as.numeric(sj[sj$novel == "N", ]$start)-1))
boxplot(d1,d2, names = c("Novel", "Other"), main = "Distance", ylab = "log10(bp)")
dev.off()


save(sj, file = paste(ro, "NSS.RData", sep = ""))
write.table(sj, file = paste(ro, "NSS.txt", sep = ""), sep = "\t", row.names = F, quote = F)
