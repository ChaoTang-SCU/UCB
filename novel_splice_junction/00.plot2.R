#===============================================================================================================================================#
#********************************************************************** Plot *******************************************************************#
#===============================================================================================================================================#


#### Novel splice junction result illustrate
## Wed Mar 28 2018
## By Tang Chao

#### Basic settings =============================================================================================================================
depar <- par()

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/")


col.a <- "#800000" # col.annotated
col.u <- "#008080" # col.unannotated
col.n <- "#CD853F" # col.novel
col.h <- "#EEE8AA" # col.HFnovel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron


require(VennDiagram)

#### Load data ==================================================================================================================================
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/EST_mRNA_bulk_junction.RData")

#### cutoff ----
sj_raw <- sj

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
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"

sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)
#   Annotated Unannotated       Novel
#      155659       72684      228939

sum(sj$classification=="Novel" & sj$CellSum>1)
#[1] 54183
sj_hf <- sj[sj$classification=="Novel" & sj$CellSum>1,]
sj_hf$classification <- "HFNovel"
sj_tu <- rbind(sj, sj_hf) 



#### 2. Basic statistic of novel splicing junction ==============================================================================================

## No. and percentage of SC SJ
pdf(paste(fo, "No. and percentage of SC SJ.pdf", sep = ""))
x = barplot(table(sj_tu$classification), col = c(col.a, col.u, col.n, col.h),
            ylim = c(0,max(table(sj_tu$classification))*1.2))
n = table(sj$classification)
p = round(n/sum(n)*100,2)
text(x = x[1:3], y = n, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)

n = table(sj_tu$classification)[4]
p = round(table(sj_tu$classification)[4]/table(sj_tu$classification)[3]*100, 2)
text(x = x[4], y = n, labels = paste(n, "\n(", p, ")", sep = ""), cex = .8, pos = 3, col = 2)
dev.off()


## Boxplot of cells and reads
pdf(paste(fo,"Boxplot of cells and reads of novel splice junction.pdf", sep = ""))
par(mfrow = c(1,2))
boxplot(log10(CellSum) ~ classification, data = sj_tu, names = F, xaxt="n",
        ylab = "Log10(cells)", col = c(col.a, col.u, col.n, col.h), pch = 20)
boxplot(log10(ReadSum) ~ classification, data = sj_tu, names = F, xaxt="n", 
        ylab = "Log10(Reads)", col = c(col.a, col.u, col.n, col.h), pch = 20)
par(depar)
dev.off()

library(ggplot2)

## Vioplot of cells and reads
f1 <- ggplot(data = sj_tu, aes(x = classification, y = log10(CellSum), fill = classification))+
  geom_violin(adjust = 2)+
  scale_fill_manual(values=c(col.a, col.u, col.n, col.h))+
  theme(legend.position="none",
  		axis.title.x = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())

f2 <- ggplot(data = sj_tu, aes(x = classification, y = log10(ReadSum), fill = classification))+
  geom_violin()+
  scale_fill_manual(values=c(col.a, col.u, col.n, col.h))+
  theme(legend.position="none",
  		axis.title.x = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())
       
pdf(paste(fo,"Vioplot of cells and reads of novel splice junction.pdf", sep = ""), width = 10, height = 5)
library(cowplot)
plot_grid(f1,f2, ncol = 2)
dev.off()


## Boxplot of novel junc distance

pdf(paste(fo,"Boxplot of novel junc distance.pdf", sep = ""))
d1 <- log10((as.numeric(sj_tu[sj_tu$classification == "Annotated", ]$end)+1) - (as.numeric(sj_tu[sj_tu$classification == "Annotated", ]$start)-1))
d2 <- log10((as.numeric(sj_tu[sj_tu$classification == "Unannotated", ]$end)+1) - (as.numeric(sj_tu[sj_tu$classification == "Unannotated", ]$start)-1))
d3 <- log10((as.numeric(sj_tu[sj_tu$classification == "Novel", ]$end)+1) - (as.numeric(sj_tu[sj_tu$classification == "Novel", ]$start)-1))
d4 <- log10((as.numeric(sj_tu[sj_tu$classification == "HFNovel", ]$end)+1) - (as.numeric(sj_tu[sj_tu$classification == "HFNovel", ]$start)-1))

boxplot(d1,d2,d3,d4, names = F, xaxt = "n", ylab = "log10(bp)", col = c(col.a, col.u, col.n, col.h))
dev.off()


## Boxplot of GC content

p1 <- ggplot(sj_tu, aes(x = classification, y = donor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")

p2 <- ggplot(sj_tu, aes(x = classification, y = acceptor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor exon")


p3 <- ggplot(sj_tu, aes(x = classification, y = donor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor intron")+
  labs(y = "GC Content")


p4 <- ggplot(sj_tu, aes(x = classification, y = acceptor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor intron")

pdf(paste(fo,"Vioplot of GC content.pdf", sep = ""))
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()



p1 <- ggplot(sj_tu, aes(x = donor_exon_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")

p2 <- ggplot(sj_tu, aes(x = acceptor_exon_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Acceptor exon")+
  labs(y = "GC Content")

p3 <- ggplot(sj_tu, aes(x = donor_intron_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor intron")+
  labs(y = "GC Content")

p4 <- ggplot(sj_tu, aes(x = acceptor_intron_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Acceptor intron")+
  labs(y = "GC Content")

pdf(paste(fo,"Density of GC content.pdf", sep = ""))
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()



#### 3. PhastCons conservation scores of annotated and novel SS =================================================================================


#### shell
cd /mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/

## Prepare bed file:
awk '$7==1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Annotated_sj_start_site.bed
awk '$7==1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Annotated_sj_end_site.bed

awk '$7==0&&$11=="N"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Uuannotated_sj_start_site.bed
awk '$7==0&&$11=="N"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Unannotated_sj_end_site.bed

awk '$7==0&&$11=="Y"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > Novel_sj_start_site.bed
awk '$7==0&&$11=="Y"' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > Novel_sj_end_site.bed

awk '$7==0&&$11=="Y"&&$9>1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$3-50"\t"$3+49}' > HFNovel_sj_start_site.bed
awk '$7==0&&$11=="Y"&&$9>1' sj_for_PhastCons.txt | awk '{print "chr"$2"\t"$4-49"\t"$4+50}' > HFNovel_sj_end_site.bed

## Extract PhastCons conservation scores from bigwig file:
# Annotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_start_site_matrix.txt

# Annotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Annotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Annotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Annotated_sj_end_site_matrix.txt


# Unannotated_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Uuannotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Unannotated_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Uuannotated_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Unannotated_sj_start_site_matrix.txt

# Unannotated_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Unannotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Unannotated_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Unannotated_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Unannotated_sj_end_site_matrix.txt


# Novel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_start_site_matrix.txt

# Novel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > Novel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 Novel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > Novel_sj_end_site_matrix.txt


# HFNovel_sj_start_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 HFNovel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > HFNovel_sj_start_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 HFNovel_sj_start_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > HFNovel_sj_start_site_matrix.txt

# HFNovel_sj_end_site
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool agg 50:50 HFNovel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw -header -expanded /dev/stdout > HFNovel_sj_end_site_agg.txt
/home/xiaqiuqi/bigwig_tools/bwtool/bwtool matrix 50:50 HFNovel_sj_end_site.bed /mnt/data8/xiaqiuqi/hg38_phyloP100way/hg38.100way.phyloP100way/hg38_100way.phastCons/phastCons100_23th_chr_fixed_step.bw /dev/stdout > HFNovel_sj_end_site_matrix.txt


library(reshape2)


AEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
AEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
ASM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
ASA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Annotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)

UEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
UEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
USM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
USA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Unannotated_sj_start_site_agg.txt", stringsAsFactors = F, header = T)

NEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
NEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
NSM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
NSA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Novel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)

HEM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_sj_end_site_matrix.txt", stringsAsFactors = F, header = F)
HEA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_sj_end_site_agg.txt", stringsAsFactors = F, header = T)
HSM <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_sj_start_site_matrix.txt", stringsAsFactors = F, header = F)
HSA <- read.table("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/HFNovel_sj_start_site_agg.txt", stringsAsFactors = F, header = T)



## PhastCons conservation scores agg:
#### Plot in one figure
rbind(ASA, USA, NSA, HSA) -> test
test$Type <- factor(rep(c("Annotated","Unannotated","Novel","HFNovel"), c(nrow(ASA),nrow(USA),nrow(NSA),nrow(HSA))),
                    levels = c("Annotated","Unannotated","Novel","HFNovel"))

rbind(AEA, UEA, NEA, HEA) -> test2
test2$Type <- factor(rep(c("Annotated","Unannotated","Novel","HFNovel"), c(nrow(AEA),nrow(UEA),nrow(NEA),nrow(HEA))),
                     levels = c("Annotated","Unannotated","Novel","HFNovel"))

rbind(test2,test) -> test3
test3$Role <- factor(rep(c("Acceptor","Donor"),c(nrow(test2),nrow(test))),
                     levels = c("Acceptor","Donor"))

pdf(paste(fo, "PhastCons conservation scores agg.pdf",sep = ""), width = 10, height = 7)
ggplot(test3, aes(x = Position, y = Mean_1, colour = Type)) + 
  geom_path(size = 1.5)+
  facet_grid(. ~ Role)+
  scale_color_manual(values = c(col.a, col.u, col.n, col.h))+
  theme(panel.grid = element_blank(),
  		legend.position = "none")+
  labs(y = "PhastCons Score")

dev.off()


## PhastCons conservation scores matrix:
library(reshape2)
p1 <- ggplot(melt(ASM), aes(x = variable, y = value, fill = col.a))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank())+
  scale_fill_manual(values = col.a)

p2 <- ggplot(melt(AEM), aes(x = variable, y = value, fill = col.a))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())+
  scale_fill_manual(values = col.a)

p3 <- ggplot(melt(USM), aes(x = variable, y = value, fill = col.u))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks = element_blank())+
  scale_fill_manual(values = col.u)

p4 <- ggplot(melt(UEM), aes(x = variable, y = value, fill = col.u))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())+
  scale_fill_manual(values = col.u)

p5 <- ggplot(melt(NSM), aes(x = variable, y = value, fill = col.n))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.y = element_blank())+
  scale_fill_manual(values = col.n)

p6 <- ggplot(melt(NEM), aes(x = variable, y = value, fill = col.n))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none")+
  scale_fill_manual(values = col.n)

p7 <- ggplot(melt(HSM), aes(x = variable, y = value, fill = col.h))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            axis.line = element_blank(),
            axis.ticks.y = element_blank())+
  scale_fill_manual(values = col.h)

p8 <- ggplot(melt(HEM), aes(x = variable, y = value, fill = col.h))+
      geom_boxplot(outlier.alpha = 0)+
      theme(axis.title = element_blank(), 
            panel.grid = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none")+
  scale_fill_manual(values = col.h)

library(cowplot)
pdf(paste(fo, "PhastCons conservation scores matrix.pdf",sep = ""), width = 10, height = 7)
plot_grid(p2, p1, p4, p3, p6, p5, p8, p7, ncol = 2)
dev.off()



#### 4. Seqlogo of splice sites =================================================================================================================

## plot seqLogo plot:
seqLogoFromMatrix <- function(tableResult, plot = TRUE, ic.scale=FALSE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15){
  ## 
  N = ncol(tableResult)
  
  if("seqLogo" %in% installed.packages()){
    library(seqLogo)
  }else{
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqLogo")
    library(seqLogo)
  }
  ## 
  if(is.matrix(tableResult)){
    data_tu <- matrix(rep(NA,4*N),nrow = 4)
    rownames(data_tu) <- c("A", "C", "G", "T")
    
    for(j in 1:N){
      for(i in c("A", "C", "G", "T")){
        data_tu[i,j] <- sum(tableResult[,j] == i)
      }
    }
    data_tu_p <- data_tu/colSums(data_tu)
    if(plot == FALSE){
      return(data_tu_p)
    }else{
      pwm <- makePWM(data_tu_p)
      
      seqLogo(pwm, ic.scale = ic.scale, xaxis = xaxis, yaxis = yaxis, xfontsize = xfontsize, yfontsize = yfontsize)
    }
  }else{
    stop("The input data of tableResult must be a matrix")
  }
}



seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 48, 56), split="")), xaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 48, 56), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 145, 153), split="")), xaxis = F, yaxis = F)
seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 145, 153), split="")), xaxis = F, yaxis = F)

dev.off()

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F)
grid.text(1:9, x = seq(0.25,0.87,0.0775),y = .08)

## function:
MyEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log2(x))
  Reduce("+", -tmp)
}

round(apply(annotated_donor_percentage, 2, MyEntropy),2)


MyRelativeEntropy <- function(input){
  tmp <- sapply(input,function(x) x*log(x/0.25))
  Reduce("+", tmp)
}


pdf(paste(fo,"Seqlogo_Shannon_entropy.pdf",sep=""), width = 7, height = 3)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

dev.off()


pdf(paste(fo,"Seqlogo_Relative_entropy.pdf",sep=""), width = 7, height = 3)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 48, 56), split="")), xaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 48, 56), split="")), xaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.25,0.87,0.0775),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Annotated"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Unannotated"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "Novel"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 145, 153), split="")), xaxis = F, yaxis = F)
x <- seqLogoFromMatrix(do.call(rbind,strsplit(substr(with(sj_tu, ParsedFasta[classification == "HFNovel"]), 145, 153), split="")), xaxis = F, yaxis = F, plot = F)
grid.text(round(apply(x, 2, MyRelativeEntropy),2), x = seq(0.14,0.86,0.09),y = .08)

dev.off()



#### 5. MaxEntScan for Maximum Entropy of splicing sites ========================================================================================


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

annotated <- sj_tu[sj_tu$classification == "Annotated", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(annotated) <- c("Donor","Acceptor")
annotated$density <- get_density(annotated$Donor, annotated$Acceptor)

unannotated <- sj_tu[sj_tu$classification == "Unannotated", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(unannotated) <- c("Donor","Acceptor")
unannotated$density <- get_density(unannotated$Donor, unannotated$Acceptor)

novel <- sj_tu[sj_tu$classification == "Novel", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(novel) <- c("Donor","Acceptor")
novel$density <- get_density(novel$Donor, novel$Acceptor)

HFnovel <- sj_tu[sj_tu$classification == "HFNovel", c("donor_MaxEnt","acceptor_MaxEnt")]
colnames(HFnovel) <- c("Donor","Acceptor")
HFnovel$density <- get_density(HFnovel$Donor, HFnovel$Acceptor)


p1 <- ggplot(annotated, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
p2 <- ggplot(unannotated, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
p3 <- ggplot(novel, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")
p4 <- ggplot(HFnovel, aes(x = Donor, y = Acceptor,color = density))+
  		geom_point(size = .6)+
  		scale_color_viridis()+
  		theme(legend.position = "none")

library(cowplot)
pdf(paste(fo, "Maximum Entropy plot of novel and other sj.pdf", sep = ""))
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()


#### 6. Stop Codon ==============================================================================================================================


pdf(paste(fo,"Boxplot of stop codons.pdf", sep = ""))
par(mfrow = c(1,2), las = 3)
### Donor
boxplot(donor_exon_sum_stop_codon ~ classification, data = sj_tu,
        boxwex = 0.25, at = 1:4 - 0.2, ylim = c(0,8.3),
        col = c(col.a, col.u, col.n, col.h), outline = F,
        main = "Donor", names = c("Exon", "Exon", "Exon", "Exon"),
        ylab = "Sum of Stop Codon")
boxplot(donor_intron_sum_stop_codon ~ classification, data = sj_tu, add = TRUE,
        boxwex = 0.25, at = 1:4 + 0.2, outline = F, names = c("Intron", "Intron", "Intron", "Intron"),
        col = c(col.a, col.u, col.n, col.h))
### Acceptor
boxplot(acceptor_exon_sum_stop_codon ~ classification, data = sj_tu,
        boxwex = 0.25, at = 1:4 - 0.2, ylim = c(0,8.3),
        col = c(col.a, col.u, col.n, col.h), outline = F,
        main = "Acceptor", names = c("Exon", "Exon", "Exon", "Exon"),
        ylab = "Sum of Stop Codon")
boxplot(acceptor_intron_sum_stop_codon ~ classification, data = sj_tu, add = TRUE,
        boxwex = 0.25, at = 1:4 + 0.2, outline = F, names = c("Intron", "Intron", "Intron", "Intron"),
        col = c(col.a, col.u, col.n, col.h))
par(depar)
dev.off()



####
#### 7. GC content trend ========================================================================================================================
####

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
dim(sj)
# [1] 870046     37
plot(density(sj$donor_exon_GC, adjust = 1.5), main = "Annotated junction", xlab = "GC Content")
lines(density(sj[sj$ReadSum>=10,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$ReadSum>=100,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$ReadSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$ReadSum>=10000,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$ReadSum>=100000,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$ReadSum>=1000000,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$ReadSum>=10000000,]$donor_exon_GC, adjust = 1.5), col = 8)

All_junction_donor_exon_GC_reads <- rbind(
data.frame(GC = sj$donor_exon_GC, ReadSum = 2),
data.frame(GC = sj[sj$ReadSum>=10,]$donor_exon_GC, ReadSum = 10),
data.frame(GC = sj[sj$ReadSum>=100,]$donor_exon_GC, ReadSum = 100),
data.frame(GC = sj[sj$ReadSum>=1000,]$donor_exon_GC, ReadSum = 1000),
data.frame(GC = sj[sj$ReadSum>=10000,]$donor_exon_GC, ReadSum = 10000),
data.frame(GC = sj[sj$ReadSum>=100000,]$donor_exon_GC, ReadSum = 100000),
data.frame(GC = sj[sj$ReadSum>=1000000,]$donor_exon_GC, ReadSum = 1000000),
data.frame(GC = sj[sj$ReadSum>=10000000,]$donor_exon_GC, ReadSum = 10000000)
)
All_junction_donor_exon_GC_reads$ReadSum <- as.factor(All_junction_donor_exon_GC_reads$ReadSum)

p1 <- ggplot(All_junction_donor_exon_GC_reads, aes(x = GC, fill = ReadSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ ReadSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Junctions")+
    labs(x = "GC Content", y = "Density")







plot(density(sj$donor_exon_GC, adjust = 1.5))
lines(density(sj[sj$CellSum>=2,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$CellSum>=10,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$CellSum>=50,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$CellSum>=100,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$CellSum>=500,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$CellSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$CellSum>=2000,]$donor_exon_GC, adjust = 1.5), col = 8)




All_junction_donor_exon_GC_cells <- rbind(
data.frame(GC = sj$donor_exon_GC, CellSum = 1),
data.frame(GC = sj[sj$CellSum>=2,]$donor_exon_GC, CellSum = 2),
data.frame(GC = sj[sj$CellSum>=10,]$donor_exon_GC, CellSum = 10),
data.frame(GC = sj[sj$CellSum>=50,]$donor_exon_GC, CellSum = 50),
data.frame(GC = sj[sj$CellSum>=100,]$donor_exon_GC, CellSum = 100),
data.frame(GC = sj[sj$CellSum>=500,]$donor_exon_GC, CellSum = 500),
data.frame(GC = sj[sj$CellSum>=1000,]$donor_exon_GC, CellSum = 1000),
data.frame(GC = sj[sj$CellSum>=2000,]$donor_exon_GC, CellSum = 2000)
)
All_junction_donor_exon_GC_cells$CellSum <- as.factor(All_junction_donor_exon_GC_cells$CellSum)

p2 <- ggplot(All_junction_donor_exon_GC_cells, aes(x = GC, fill = CellSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ CellSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Junctions")+
    labs(x = "GC Content", y = "Density")







plot(density(sj[sj$annotation==1, ]$donor_exon_GC, adjust = 1.5), ylim = c(0,4))
lines(density(sj[sj$annotation==1 & sj$ReadSum>=10,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=100,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=10000,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=100000,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=1000000,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$annotation==1 & sj$ReadSum>=10000000,]$donor_exon_GC, adjust = 1.5), col = 8)


Anotated_junction_donor_exon_GC_reads <- rbind(
data.frame(GC = sj[sj$annotation==1, ]$donor_exon_GC, ReadSum = 2),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=10,]$donor_exon_GC, ReadSum = 10),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=100,]$donor_exon_GC, ReadSum = 100),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=1000,]$donor_exon_GC, ReadSum = 1000),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=10000,]$donor_exon_GC, ReadSum = 10000),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=100000,]$donor_exon_GC, ReadSum = 100000),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=1000000,]$donor_exon_GC, ReadSum = 1000000),
data.frame(GC = sj[sj$annotation==1 & sj$ReadSum>=10000000,]$donor_exon_GC, ReadSum = 10000000)
)
Anotated_junction_donor_exon_GC_reads$ReadSum <- as.factor(Anotated_junction_donor_exon_GC_reads$ReadSum)

p3 <- ggplot(Anotated_junction_donor_exon_GC_reads, aes(x = GC, fill = ReadSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ ReadSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Annotated Junctions")+
    labs(x = "GC Content", y = "Density")




plot(density(sj[sj$annotation==1, ]$donor_exon_GC, adjust = 1.5), ylim = c(0,4))
lines(density(sj[sj$annotation==1 & sj$CellSum>=2,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==1 & sj$CellSum>=10,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==1 & sj$CellSum>=50,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==1 & sj$CellSum>=100,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==1 & sj$CellSum>=500,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==1 & sj$CellSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$annotation==1 & sj$CellSum>=2000,]$donor_exon_GC, adjust = 1.5), col = 8)

Anotated_junction_donor_exon_GC_cells <- rbind(
data.frame(GC = sj[sj$annotation==1, ]$donor_exon_GC, CellSum = 1),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=2,]$donor_exon_GC, CellSum = 2),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=10,]$donor_exon_GC, CellSum = 10),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=50,]$donor_exon_GC, CellSum = 50),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=100,]$donor_exon_GC, CellSum = 100),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=500,]$donor_exon_GC, CellSum = 500),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=1000,]$donor_exon_GC, CellSum = 1000),
data.frame(GC = sj[sj$annotation==1 & sj$CellSum>=2000,]$donor_exon_GC, CellSum = 2000)
)
Anotated_junction_donor_exon_GC_cells$CellSum <- as.factor(Anotated_junction_donor_exon_GC_cells$CellSum)
p4 <- ggplot(Anotated_junction_donor_exon_GC_cells, aes(x = GC, fill = CellSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ CellSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Annotated Junctions")+
    labs(x = "GC Content", y = "Density")







plot(density(sj[sj$annotation==0, ]$donor_exon_GC, adjust = 1.5), ylim = c(0,5))
lines(density(sj[sj$annotation==0 & sj$ReadSum>=10,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==0 & sj$ReadSum>=100,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==0 & sj$ReadSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==0 & sj$ReadSum>=10000,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==0 & sj$ReadSum>=100000,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==0 & sj$ReadSum>=500000,]$donor_exon_GC, adjust = 1.5), col = 7)

Unannotated_junction_donor_exon_GC_reads <- rbind(
data.frame(GC = sj[sj$annotation==0, ]$donor_exon_GC, ReadSum = 2),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=10,]$donor_exon_GC, ReadSum = 10),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=100,]$donor_exon_GC, ReadSum = 100),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=1000,]$donor_exon_GC, ReadSum = 1000),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=10000,]$donor_exon_GC, ReadSum = 10000),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=100000,]$donor_exon_GC, ReadSum = 100000),
data.frame(GC = sj[sj$annotation==0 & sj$ReadSum>=500000,]$donor_exon_GC, ReadSum = 500000)
)
Unannotated_junction_donor_exon_GC_reads$ReadSum <- as.factor(Unannotated_junction_donor_exon_GC_reads$ReadSum)
p5 <- ggplot(Unannotated_junction_donor_exon_GC_reads, aes(x = GC, fill = ReadSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ ReadSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Unannotated Junctions")+
    labs(x = "GC Content", y = "Density")





plot(density(sj[sj$annotation==0, ]$donor_exon_GC, adjust = 1.5), ylim = c(0,5))
lines(density(sj[sj$annotation==0 & sj$CellSum>=2,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==0 & sj$CellSum>=10,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==0 & sj$CellSum>=50,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==0 & sj$CellSum>=100,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==0 & sj$CellSum>=500,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==0 & sj$CellSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$annotation==0 & sj$CellSum>=2000,]$donor_exon_GC, adjust = 1.5), col = 8)

Unannotated_junction_donor_exon_GC_cells <- rbind(
data.frame(GC = sj[sj$annotation==0, ]$donor_exon_GC, CellSum = 1),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=2,]$donor_exon_GC, CellSum = 2),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=10,]$donor_exon_GC, CellSum = 10),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=50,]$donor_exon_GC, CellSum = 50),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=100,]$donor_exon_GC, CellSum = 100),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=500,]$donor_exon_GC, CellSum = 500),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=1000,]$donor_exon_GC, CellSum = 1000),
data.frame(GC = sj[sj$annotation==0 & sj$CellSum>=2000,]$donor_exon_GC, CellSum = 2000)
)
Unannotated_junction_donor_exon_GC_cells$CellSum <- as.factor(Unannotated_junction_donor_exon_GC_cells$CellSum)
p6 <- ggplot(Unannotated_junction_donor_exon_GC_cells, aes(x = GC, fill = CellSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ CellSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("All Unannotated Junctions")+
    labs(x = "GC Content", y = "Density")
	




plot(density(sj[sj$annotation==0 & sj$SCSpe == "Y", ]$donor_exon_GC, adjust = 1.5), ylim = c(0,5))
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=10,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=100,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=2000,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=5000,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=6000,]$donor_exon_GC, adjust = 1.5), col = 7)

Novel_junction_donor_exon_GC_reads <- rbind(
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y", ]$donor_exon_GC, ReadSum = 2),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=10,]$donor_exon_GC, ReadSum = 10),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=50,]$donor_exon_GC, ReadSum = 50),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=100,]$donor_exon_GC, ReadSum = 100),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=500,]$donor_exon_GC, ReadSum = 500),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=1000,]$donor_exon_GC, ReadSum = 1000),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=5000,]$donor_exon_GC, ReadSum = 5000),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$ReadSum>=10000,]$donor_exon_GC, ReadSum = 10000)
)
Novel_junction_donor_exon_GC_reads$ReadSum <- as.factor(Novel_junction_donor_exon_GC_reads$ReadSum)
p7 <- ggplot(Novel_junction_donor_exon_GC_reads, aes(x = GC, fill = ReadSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ ReadSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("Novel Junctions")+
    labs(x = "GC Content", y = "Density")






plot(density(sj[sj$annotation==0 & sj$SCSpe == "Y", ]$donor_exon_GC, adjust = 1.5), ylim = c(0,5))
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=2,]$donor_exon_GC, adjust = 1.5), col = 2)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=10,]$donor_exon_GC, adjust = 1.5), col = 3)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=50,]$donor_exon_GC, adjust = 1.5), col = 4)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=100,]$donor_exon_GC, adjust = 1.5), col = 5)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=500,]$donor_exon_GC, adjust = 1.5), col = 6)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=800,]$donor_exon_GC, adjust = 1.5), col = 7)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=1000,]$donor_exon_GC, adjust = 1.5), col = 8)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=1200,]$donor_exon_GC, adjust = 1.5), col = 9)
lines(density(sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=1500,]$donor_exon_GC, adjust = 1.5), col = 10)

Novel_junction_donor_exon_GC_cells <- rbind(
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y", ]$donor_exon_GC, CellSum = 1),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=2,]$donor_exon_GC, CellSum = 2),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=10,]$donor_exon_GC, CellSum = 10),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=50,]$donor_exon_GC, CellSum = 50),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=100,]$donor_exon_GC, CellSum = 100),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=500,]$donor_exon_GC, CellSum = 500),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=800,]$donor_exon_GC, CellSum = 800),
data.frame(GC = sj[sj$annotation==0 & sj$SCSpe == "Y" & sj$CellSum>=1000,]$donor_exon_GC, CellSum = 1000)
)
Novel_junction_donor_exon_GC_cells$CellSum <- as.factor(Novel_junction_donor_exon_GC_cells$CellSum)
p8 <- ggplot(Novel_junction_donor_exon_GC_cells, aes(x = GC, fill = CellSum))+
	geom_density(adjust = 1.5)+
	facet_wrap(~ CellSum, ncol = 1, strip.position = "right")+
	theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())+
	guides(fill=FALSE)+
	ggtitle("Novel Junctions")+
    labs(x = "GC Content", y = "Density")

pdf(paste(fo,"GC Content trend.pdf", sep = ""), width = 7, height = 28)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2)
dev.off()


