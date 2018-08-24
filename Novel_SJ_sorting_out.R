#### Novel splice junction new version
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table_v2/Figure3/")

col.a <- "#800000" # col.annotated
col.u <- "#008080" # col.unannotated
col.h <- "#CD853F" # col.HFnovel
col.l <- "#EEE8AA" # col.LFnovel

col.e <- "#00CCCC" # col.exon
col.i <- "#FF8C00" # col.intron

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
dim(sj)
#[1] 870046     45
table(sj$annotation)
#      0      1
# 693517 176529
sum(sj$annotation == 0 & is.na(sj[,42]) & is.na(sj[,43]))
# [1] 594054
sum(sj$annotation == 0 & is.na(sj[,42]) & !is.na(sj[,43]))
# [1] 38209
sum(sj$annotation == 0 & !is.na(sj[,42]) & is.na(sj[,43]))
# [1] 38287
sum(sj$annotation == 0 & !is.na(sj[,42]) & !is.na(sj[,43]))
# [1] 22967


sj <- sj[sj$ReadSum >= 10 & sj$overhang >= 10, ]
dim(sj)
# [1] 477114     45
table(sj$annotation)
#      0      1
# 319334 157780
sum(sj$annotation == 0 & is.na(sj[,42]) & is.na(sj[,43]))
# [1] 231476
sum(sj$annotation == 0 & is.na(sj[,42]) & !is.na(sj[,43]))
# [1] 33592
sum(sj$annotation == 0 & !is.na(sj[,42]) & is.na(sj[,43]))
# [1] 33628
sum(sj$annotation == 0 & !is.na(sj[,42]) & !is.na(sj[,43]))
# [1] 20638



sj$annotation_type <- NA
sj[sj$annotation == 1, "annotation_type"] <- "Known junction"
sj[sj$annotation == 0 & is.na(sj[,42]) & !is.na(sj[,43]), "annotation_type"] <- "Unannotated donor"
sj[sj$annotation == 0 & !is.na(sj[,42]) & is.na(sj[,43]), "annotation_type"] <- "Unannotated acceptor"
sj[sj$annotation == 0 & !is.na(sj[,42]) & !is.na(sj[,43]), "annotation_type"] <- "Unannotated pattern"
sj[sj$annotation == 0 & is.na(sj[,42]) & is.na(sj[,43]), "annotation_type"] <- "Unannotated junction"
sj$annotation_type <- factor(sj$annotation_type, levels = c("Known junction", "Unannotated acceptor", "Unannotated donor", "Unannotated pattern", "Unannotated junction"))
table(sj$annotation_type)
#      Known junction Unannotated acceptor    Unannotated donor
#              157780                33628                33592
# Unannotated pattern Unannotated junction
#               20638               231476

pie(c(231476,157780,33628,33592,20638), 
    labels = c("231,476","157,780","33,628","33,592","20,638"), 
    col = c("purple", "red", "green", "blue","orange"))


library(ggplot2)
pdf(paste(fo,"annotation type comparison.pdf", sep = ""))

ggplot(sj, aes(x = annotation_type, y = log10(ReadSum+1)))+
    geom_violin()

ggplot(sj, aes(x = annotation_type, y = log2(CellSum+1)))+
    geom_violin()

ggplot(sj, aes(x = annotation_type, y = log2(CellSum+1)))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log2(abs(end-start)+1)))+
    geom_boxplot()+
    ylab("log2(Distance)")

ggplot(sj, aes(x = annotation_type, y = donor_MaxEnt))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_MaxEnt))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = donor_exon_sum_stop_codon))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = donor_exon_min_stop_codon))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = exon_Shannon_entropy))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = intron_Relative_entropy))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_exon_GC))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_intron_GC))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log10(upstream_exon_length+1)))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log10(downstream_exon_length+1)))+
    geom_boxplot()

dev.off()


sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < 2)
[1] 168842
sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < 5)
[1] 211973
sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < 10)
[1] 223440

mapply(function(x) {sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < x)}, 2:10)

pdf(paste(fo,"cellsum cutoff of unannotated junction.pdf", sep = ""))
plot(x = 2:20, y = mapply(function(x) {sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < x)}, 2:20), 
     ylab = "No. Low Quality Splice Junction", xlab = "No. Min Cells")
lines(x = 2:20, y = mapply(function(x) {sum(sj$annotation_type == "Unannotated junction" & sj$CellSum < x)}, 2:20))
dev.off()

#sj <- sj[!(sj$annotation_type == "Unannotated junction" & sj$CellSum < 10), ]
#table(sj$annotation_type)
##      Known junction Unannotated acceptor    Unannotated donor
##              157780                33628                33592
## Unannotated pattern Unannotated junction
##               20638                 8036

pdf(paste(fo,"annotation type comparison after cellsum cutoff.pdf", sep = ""))

ggplot(sj, aes(x = annotation_type, y = log10(ReadSum+1)))+
    geom_violin()

ggplot(sj, aes(x = annotation_type, y = log2(CellSum+1)))+
    geom_violin()

ggplot(sj, aes(x = annotation_type, y = log2(CellSum+1)))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log2(abs(end-start)+1)))+
    geom_boxplot()+
    ylab("log2(Distance)")

ggplot(sj, aes(x = annotation_type, y = donor_MaxEnt))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_MaxEnt))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = donor_exon_sum_stop_codon))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = donor_exon_min_stop_codon))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = exon_Shannon_entropy))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = intron_Relative_entropy))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_exon_GC))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = acceptor_intron_GC))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log10(upstream_exon_length+1)))+
    geom_boxplot()

ggplot(sj, aes(x = annotation_type, y = log10(downstream_exon_length+1)))+
    geom_boxplot()

dev.off()





sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj[sj$novel=="Y" & sj$CellSum>1, ]$classification <- "HFNovel"
sj[sj$novel=="Y" & sj$CellSum==1, ]$classification <- "LFNovel"

sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","HFNovel","LFNovel"))
table(sj$classification)
#  Annotated Unannotated     HFNovel     LFNovel
#     157780       73771       60985      184578


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


sc_SJ <- as.character(sj[sj$annotation == 0, ]$sj)
venn_plot <- venn.diagram(list(mRNA = mRNA_SJ, EST = EST_SJ, Bulk = bulk_SJ, SC = sc_SJ), fill = rainbow(4), filename = NULL, main.cex = 2, main = "Venn plot of splicing junction in different sequencing methods", height = 6, width = 10)
pdf(paste(fo,"venn_plot_of_sc_annotated_validated_by_bulk_mrna_EST_SJ.pdf",sep = ""), width = 10, height = 7)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


sum(sj$annotation == 0)
#[1] 319334
table(sj[sj$annotation == 0, "classification"])
#  Annotated Unannotated     HFNovel     LFNovel
#          0       73771       60985      184578

dim(sj[(sj$annotation == 0 & !is.na(sj[25]) & is.na(sj[27])) | (sj$annotation == 0 & !is.na(sj[26]) & is.na(sj[27])), ])
#[1] 1376   47
dim(sj[(sj$annotation == 0 & is.na(sj[25]) & is.na(sj[26]) & !is.na(sj[27])), ])
#[1] 66099    47
dim(sj[sj$annotation == 0 & !is.na(sj[27]) & (!is.na(sj[25]) | !is.na(sj[26])), ])
#[1] 6296   47


pie(c(Bulk = 66099, "EST/mRNA" = 1376, "Bulk|EST/mRNA" = 6296, SC = 245563))




#### 2. Basic statistic of novel splicing junction ==============================================================================================



## No. and percentage of SC SJ
pdf(paste(fo, "No. and percentage of SC SJ.pdf", sep = ""))
x = barplot(table(sj$classification), col = c(col.a, col.u, col.h, col.l),
            ylim = c(0,max(table(sj$classification))*1.2))
n = table(sj$classification)
p = round(n/sum(n)*100,2)
text(x = x, y = n, labels = paste(n,"\n",p,sep = ""), cex = .8, pos = 3)

dev.off()


## Boxplot of cells and reads
pdf(paste(fo,"Boxplot of cells and reads of novel splice junction.pdf", sep = ""))
par(mfrow = c(1,2))
boxplot(log10(CellSum) ~ classification, data = sj, names = F, xaxt="n",
        ylab = "Log10(cells)", col = c(col.a, col.u, col.h, col.l), pch = 20)
boxplot(log10(ReadSum) ~ classification, data = sj, names = F, xaxt="n", 
        ylab = "Log10(Reads)", col = c(col.a, col.u, col.h, col.l), pch = 20)
par(depar)
dev.off()



f1 <- ggplot(data = sj, aes(x = classification, y = log10(CellSum), fill = classification))+
  geom_violin(adjust = 2)+
  scale_fill_manual(values=c(col.a, col.u, col.h, col.l))+
  theme(legend.position="none",
        axis.title.x = element_blank(), 
        panel.grid = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())

f2 <- ggplot(data = sj, aes(x = classification, y = log10(ReadSum), fill = classification))+
  geom_violin()+
  scale_fill_manual(values=c(col.a, col.u, col.h, col.l))+
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
ggplot(sj, aes(x = classification, y = log2(abs(end-start)+1), fill = classification))+
    geom_violin()+
    scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
    ylab("log2(Distance)")+
    theme(legend.position = "none")
dev.off()



## Boxplot of GC content



p1 <- ggplot(sj, aes(x = classification, y = donor_exon_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
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
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor exon")


p3 <- ggplot(sj, aes(x = classification, y = donor_intron_GC, fill = classification))+
  geom_violin(adjust = 1.5)+
  scale_fill_manual(values = c(col.a, col.u, col.h, col.l))+
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
  theme(legend.position = "none",
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  ggtitle("Acceptor intron")

pdf(paste(fo,"Vioplot of GC content.pdf", sep = ""))
plot_grid(p1,p2,p3,p4, ncol = 2)
dev.off()



p1 <- ggplot(sj, aes(x = donor_exon_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.h, col.l))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor exon")+
  labs(y = "GC Content")

p2 <- ggplot(sj, aes(x = acceptor_exon_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.h, col.l))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Acceptor exon")+
  labs(y = "GC Content")

p3 <- ggplot(sj, aes(x = donor_intron_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.h, col.l))+
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Donor intron")+
  labs(y = "GC Content")

p4 <- ggplot(sj, aes(x = acceptor_intron_GC, colour = classification))+
  geom_density(adjust = 1.5)+
  scale_colour_manual(values = c(col.a, col.u, col.h, col.l))+
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





















