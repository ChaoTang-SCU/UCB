#### splicing junction basic statistic
## Sat Mar 24 16:09:04 CST 2018
## By Tang Chao


#### 0. Basic settings ======================================================================


depar <- par()
setwd("/mnt/data5/BGI/UCB/tangchao/novel_SS")

## Figure output to:
fo <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/Figure/")

## RData output to:
ro <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/")

## Table output to:
to <- file.path("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/")


#### 1. Load data ===========================================================================


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/SJ/SJ_tu.RData")


#### 2. Basic statistic =====================================================================


with(sj,summary(ReadSum[annotation==0]))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#      2.0       3.0       7.0     286.2      74.0 2709681.0
with(sj,summary(ReadSum[annotation==1]))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#        2       328      2609     65350     15490 189450204
with(sj,summary(CellSum[annotation==1]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.00    3.00   14.00   96.17   64.00 3574.00
with(sj,summary(CellSum[annotation==0]))
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   1.000    1.000    1.000    2.131    1.000 2696.000 


## Bar chart of annotation
pdf(paste(fo,"Bar chart of annotation.pdf", sep = ""))
x = barplot(table(sj$annotation)/nrow(sj)*100, names.arg = c("Unannotated","Annotated"), ylim = c(0,100), main = "Annotation", ylab = "%")
p <- round(table(sj$annotation)/nrow(sj)*100,2)
n <- table(sj$annotation)
text(x = x, y = p+8, labels = paste(n,"\n",p,sep = ""), cex = .8)
dev.off()


## Boxplot of cells and reads
pdf(paste(fo,"Boxplot of cells and reads.pdf", sep = ""))
par(mfrow = c(1,2))
with(sj,boxplot(log10(CellSum[annotation==1]),log10(CellSum[annotation==0]),
                names = c("Annotated", "Unannotated"), main = "No. of cell", 
                ylab = "Log10(cells)"))

with(sj,boxplot(log10(ReadSum[annotation==1]),log10(ReadSum[annotation==0]),
                names = c("Annotated", "Unannotated"), main = "No. of reads", 
                ylab = "Log10(reads)"))
par(depar)
dev.off()


## Barchart of motif
pdf(paste(fo,"Barchart of motif.pdf", sep = ""), width = 10, height = 7)
par(mfrow = c(2,1))
par(mar = c(1,4,1,1))
x = with(sj, barplot(table(motif[annotation==1])/sum(table(motif[annotation==1]))*100, names.arg = F, ylim = c(0,60), main = "SS motif", ylab = "%"))
p <- round(with(sj, table(motif[annotation==1]))/sum(with(sj, table(motif[annotation==1])))*100,2)
n <- with(sj, table(motif[annotation==1]))
text(x = x, y = p+5, labels = paste(n,"\n",p,sep = ""), cex = .8)
legend("topright", legend = "Annotated", bty = "n")

par(mar = c(3,4,1,1))
x = with(sj, barplot(table(motif[annotation==0])/sum(table(motif[annotation==0]))*100, ylim = c(0,42), ylab = "%",
                     names.arg = c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT")))
p <- round(with(sj, table(motif[annotation==0]))/sum(with(sj, table(motif[annotation==0])))*100,2)
n <- with(sj, table(motif[annotation==0]))
text(x = x, y = p+4, labels = paste(n,"\n",p,sep = ""), cex = .8)
legend("topright", legend = "Unannotated", bty = "n")
par(depar)
dev.off()

## Boxplot of junc distance
pdf(paste(fo,"Boxplot of junc distance.pdf", sep = ""))
d1 <- log10((as.numeric(sj[sj$annotation == 1, ]$end)+1) - (as.numeric(sj[sj$annotation == 1, ]$start)-1))
d2 <- log10((as.numeric(sj[sj$annotation == 0, ]$end)+1) - (as.numeric(sj[sj$annotation == 0, ]$start)-1))
boxplot(d1,d2, names = c("Annotated", "Unannotated"), main = "Distance", ylab = "log10(bp)")
dev.off()


#### 3. Parse strand ========================================================================


head(sj)
table(sj$strand)
# strand (0: undefined, 1: +, 2: -)
# intron motif: 0: non-canonical; 
# 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
table(sj[sj$strand==1,]$motif)
#      1      3      5 
# 304093  47719   1325 
table(sj[sj$strand==2,]$motif)
#      2      4      6 
# 312939  51029   1166 
table(sj[sj$strand==0,]$motif)
#      0 
# 151780 
sj -> sj_sp
sj_sp[sj_sp$motif==1,]$motif <- "GT/AG"
sj_sp[sj_sp$motif==2,]$motif <- "GT/AG"
sj_sp[sj_sp$motif==3,]$motif <- "GC/AG"
sj_sp[sj_sp$motif==4,]$motif <- "GC/AG"
sj_sp[sj_sp$motif==5,]$motif <- "AT/AC"
sj_sp[sj_sp$motif==6,]$motif <- "AT/AC"
sj_sp[sj_sp$motif==0,]$motif <- "non-canonical"

sj_sp$motif <- as.factor(sj_sp$motif)


## Barchart of motif Parse strand
pdf(paste(fo,"Barchart of motif parsed.pdf", sep = ""), width = 10, height = 7)
par(mfrow = c(2,1))
par(mar = c(1,4,1,1))
x = with(sj_sp, barplot(table(motif[annotation==1])/sum(table(motif[annotation==1]))*100, names.arg = F, ylim = c(0,112), main = "SS motif", ylab = "%"))
p <- round(with(sj_sp, table(motif[annotation==1]))/sum(with(sj_sp, table(motif[annotation==1])))*100,2)
n <- with(sj_sp, table(motif[annotation==1]))
text(x = x, y = p+5, labels = paste(n,"\n",p,sep = ""), cex = .8)
legend("topright", legend = "Annotated", bty = "n")

par(mar = c(3,4,1,1))
x = with(sj_sp, barplot(table(motif[annotation==0])/sum(table(motif[annotation==0]))*100, ylim = c(0,100), ylab = "%"))
p <- round(with(sj_sp, table(motif[annotation==0]))/sum(with(sj_sp, table(motif[annotation==0])))*100,2)
n <- with(sj_sp, table(motif[annotation==0]))
text(x = x, y = p+5, labels = paste(n,"\n",p,sep = ""), cex = .8)
legend("topright", legend = "Unannotated", bty = "n")
par(depar)
dev.off()


