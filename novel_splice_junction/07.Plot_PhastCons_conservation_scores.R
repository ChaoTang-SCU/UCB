#### Plot PhastCons conservation scores
## Sat Mar 25 2018
## By Tang Chao

#### 0. Basic settings ======================================================================


depar <- par()
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/")

require(ggplot2)
library(reshape2)

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/")


#### 1. Load data ===========================================================================


AEM <- read.table("Annotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
AEA <- read.table("Annotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
ASM <- read.table("Annotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
ASA <- read.table("Annotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)
NEM <- read.table("Novel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
NEA <- read.table("Novel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
NSM <- read.table("Novel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
NSA <- read.table("Novel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)


## PhastCons conservation scores agg:
pdf(paste(fo, "PhastCons conservation scores agg.pdf",sep = ""), width = 10, height = 7)
ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()

ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))

dev.off()


## PhastCons conservation scores matrix:
pdf(paste(fo, "PhastCons conservation scores matrix.pdf",sep = ""), width = 10, height = 7)
ggplot(melt(ASM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(AEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NSM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
dev.off()





#### cutoff: 10 reads

setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffReads10")


## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffReads10/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffReads10/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffReads10/")


#### 1. Load data ===========================================================================


AEM <- read.table("Annotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
AEA <- read.table("Annotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
ASM <- read.table("Annotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
ASA <- read.table("Annotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)
NEM <- read.table("Novel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
NEA <- read.table("Novel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
NSM <- read.table("Novel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
NSA <- read.table("Novel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)


## PhastCons conservation scores agg:
pdf(paste(fo, "PhastCons conservation scores agg reads10.pdf",sep = ""), width = 10, height = 7)
ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()

ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))

dev.off()


## PhastCons conservation scores matrix:
pdf(paste(fo, "PhastCons conservation scores matrix reads10.pdf",sep = ""), width = 10, height = 7)
ggplot(melt(ASM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(AEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NSM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
dev.off()





#### cutoff: 10 cells


setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffCells10")


## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffCells10/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffCells10/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/conservation/phastCons/cutoffCells10/")


#### 1. Load data ===========================================================================


AEM <- read.table("Annotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
AEA <- read.table("Annotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
ASM <- read.table("Annotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
ASA <- read.table("Annotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)
NEM <- read.table("Novel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
NEA <- read.table("Novel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
NSM <- read.table("Novel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
NSA <- read.table("Novel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)


## PhastCons conservation scores agg:
pdf(paste(fo, "PhastCons conservation scores agg cells10.pdf",sep = ""), width = 10, height = 7)
ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()

ggplot(ASA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(AEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NSA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))
ggplot(NEA, aes(x = Position, y = Mean_1))+
  geom_line()+
  coord_cartesian(ylim=c(0,1))

dev.off()


## PhastCons conservation scores matrix:
pdf(paste(fo, "PhastCons conservation scores matrix cells10.pdf",sep = ""), width = 10, height = 7)
ggplot(melt(ASM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(AEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NSM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
ggplot(melt(NEM), aes(x = variable, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  theme(axis.title = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
dev.off()



