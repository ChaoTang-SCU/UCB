
## investigate the alternative of CD45
## CD45: 1:198638671-198757283_+    gene_id "ENSG00000081237"; gene_name "PTPRC"

#### IR ========================================================================================================================================================

load(file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_GTF_merge/IR_known_introns_match.RData")
load(file = "/mnt/data5/BGI/UCB/tangchao/IR/IR_merge2/IR_irl_tu.RData")

sum(IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283)
# [1] 50
table(IR_GTF[IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283, ]$gene_name)
# PTPRC
#    50

IR_GTF[IR_GTF$chr == 1 & IR_GTF$start >= 198638671 & IR_GTF$end <= 198757283, ] -> IR_CD45
length(unique(IR_CD45$event_id))
# [1] 13
length(unique(IR_CD45$transcript_id))
# [1] 13

IR_irl_CD45 <- IR_irl_tu[unique(IR_CD45$event_id),]
rowSums(!is.na(IR_irl_CD45))
#1:198639138-198639225 1:198696910-198699563 1:198699705-198702386
#                    2                    23                     7
#1:198699705-198703297 1:198702531-198703297 1:198706953-198708132
#                    6                    33                    33
#1:198709825-198712952 1:198731727-198732299 1:198734426-198735126
#                   17                     4                     3
#1:198742027-198742231 1:198752372-198752593 1:198752773-198754268
#                   11                     2                     4
#1:198754405-198755905
#                   15

library(reshape2)
summary(melt(IR_irl_CD45)$value)
summary(na.omit(melt(IR_irl_CD45)$value))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.004451 0.083175 0.201347 0.361122 0.600254 1.579303


apply(IR_irl_CD45,1,function(x) which(!is.na(x)))
apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))])
table(table(as.character(unlist(apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))])))))
#   1   2   3
# 133   9   3
## Only 3 cells have 3 IR events, and 133 cells only have one IR events.

apply(IR_irl_CD45,1,function(x) x[!is.na(x)])



CD45_IR <- data.frame(event_id = rep(names(rowSums(!is.na(IR_irl_CD45))), as.numeric(rowSums(!is.na(IR_irl_CD45)))),
	cell = as.character(unlist(apply(IR_irl_CD45,1,function(x) colnames(IR_irl_CD45)[which(!is.na(x))]))),
	psi = as.character(unlist(apply(IR_irl_CD45,1,function(x) x[!is.na(x)]))))


gtf <- read.table("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf",sep="\t",header=F, stringsAsFactors = F)
transcript <- subset(gtf,gtf$V3=="transcript")
test <- as.data.frame(do.call(rbind, strsplit(transcript$V9, split = "[; ]")))

valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

transcript_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

merge(x = IR_CD45, y = transcript_infor, by = "transcript_id", all)
table(IR_CD45$transcript_biotype)
#      non_stop_decay processed_transcript       protein_coding
#                   7                    3                   37
#     retained_intron
#                   3
table(IR_CD45[!duplicated(IR_CD45$transcript_id), c("transcript_id", "transcript_biotype")]$transcript_biotype)
#      non_stop_decay processed_transcript       protein_coding
#                   1                    2                    8
#     retained_intron
#                   2




#### SE ========================================================================================================================================================

load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_GTF.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_type_and_PSI.RData")

sum(do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,1] == 1 & 
	do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,2] >= 198638671 & 
	do.call(rbind,strsplit(all_cells_SE_type_and_PSI$ASrange, split="[:-]"))[,3] <= 198757283)
# [1] 0

sum(SE_GTF$V1 == 1 & SE_GTF$V4 >= 198638671 & SE_GTF$V5 <= 198757283, na.rm=T)
# [1] 72

SE_GTF_CD45 <- SE_GTF[which(SE_GTF$V1 == 1 & SE_GTF$V4 >= 198638671 & SE_GTF$V5 <= 198757283), ]
dim(SE_GTF_CD45)
# [1] 72 21
table(SE_GTF_CD45$AStype)
#      exon-skipping_exactly     exon-skipping_half-exon
#                         29                           2
#exon-skipping_multipul-exon
#                         41


length(unique(SE_GTF_CD45$loci))
# [1] 16


test <- as.data.frame(do.call(rbind, strsplit(SE_GTF_CD45$V9, split = "[; ]")))

valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
           transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
           gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
           transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
           transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

dim(gtf_infor)
# [1] 72  6

SE_GTF_CD45 <- cbind(SE_GTF_CD45,gtf_infor)
table(SE_GTF_CD45[!duplicated(SE_GTF_CD45$transcript_id),]$transcript_biotype)
# protein_coding retained_intron
#              5               1

SE_PSI_CD45 <- all_cells_SE_type_and_PSI[all_cells_SE_type_and_PSI$loci %in% SE_GTF_CD45$loci, ]
dim(SE_PSI_CD45)
# [1]   16 3585

row.names(SE_PSI_CD45) <- SE_PSI_CD45$loci
SE_PSI_CD45 <- SE_PSI_CD45[,12:ncol(SE_PSI_CD45)]
rowSums(!is.na(SE_PSI_CD45))
#1:198692374-198696711-198696910-198699563
#                                      600
#1:198699705-198702386-198702531-198703297
#                                      735
#1:198708262-198709686-198709825-198712952
#                                       28
#1:198692374-198699563-198699705-198703297
#                                       19
#1:198692374-198694017-198696910-198699563
#                                       58
#1:198692374-198699563-198702531-198703297
#                                       15
#1:198692374-198694017-198694141-198696711
#                                       14
#1:198692374-198696711-198699705-198703297
#                                        2
#1:198692374-198696711-198698932-198699563
#                                        2
#1:198692374-198696711-198702531-198703297
#                                       12
#1:198748200-198749415-198749550-198750491
#                                        2
#1:198752773-198754268-198754405-198755905
#                                        1
#1:198699705-198702349-198702531-198703297
#                                        4
#1:198706953-198708132-198709825-198712952
#                                        1
#1:198692374-198694017-198694141-198699563
#                                        1
#1:198699705-198701774-198702531-198703297
#                                        1

summary(na.omit(melt(SE_PSI_CD45)$value))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0006127 0.1656313 0.5384609 0.5303473 0.9459642 0.9995992

table(table(as.character(unlist(apply(SE_PSI_CD45,1,function(x) colnames(SE_PSI_CD45)[which(!is.na(x))])))))
#   1   2   3
# 732 365  11
## Only 11 cells have 3 IR events, and 732 cells only have one IR events.


(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",!is.na()])
intersect(colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",])],colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198699705-198702386-198702531-198703297",])])
length(intersect(colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198692374-198696711-198696910-198699563",])],colnames(SE_PSI_CD45)[!is.na(SE_PSI_CD45["1:198699705-198702386-198702531-198703297",])]))
# [1] 324
## 600 intersect 735 = 324

SE_GTF_CD45[SE_GTF_CD45$loci == "1:198692374-198696711-198696910-198699563", ]$AStype
#  "exon-skipping_exactly"
SE_GTF_CD45[SE_GTF_CD45$loci == "1:198699705-198702386-198702531-198703297", ]$AStype
#  "exon-skipping_exactly"






#### A3SS ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A3SS/RData/Parsed_A3SS_PSI_GTF.RData")

A3SS_GTF <- na.omit(A3SS_GTF)
sum(A3SS_GTF$chr == 1 & A3SS_GTF$start >= 198638671 & A3SS_GTF$end <= 198757283, na.rm=T)
# [1] 0





#### A5SS ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/A5SS/RData/Parsed_A5SS_PSI_GTF.RData")

A5SS_GTF <- na.omit(A5SS_GTF)
sum(A5SS_GTF$chr == 1 & A5SS_GTF$start >= 198638671 & A5SS_GTF$end <= 198757283, na.rm=T)
# [1] 0









#### MXE ========================================================================================================================================================




load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/MXE/RData/all_cells_MXE_type_and_existing.RData")

A5SS_GTF <- na.omit(A5SS_GTF)
sum(all_cells_MXE_type_and_existing$chr == 1 & all_cells_MXE_type_and_existing$start >= 198638671 & all_cells_MXE_type_and_existing$end <= 198757283, na.rm=T)
# [1] 2
## There are two MXE in CD45

MXE_CD45 <- all_cells_MXE_type_and_existing[all_cells_MXE_type_and_existing$chr == 1 & all_cells_MXE_type_and_existing$start >= 198638671 & all_cells_MXE_type_and_existing$end <= 198757283,]
dim(MXE_CD45)
# [1]    2 3590

rowSums(MXE_CD45[,17:ncol(MXE_CD45)] == 1)
# 1206 1701
#    2    1
## It makes no sense








#### plot in tsne ======================================================================================================================




load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

dim(IR_irl_CD45)
dim(SE_PSI_CD45)



identical(colnames(SE_PSI_CD45), colnames(IR_irl_CD45))
rbind(SE_PSI_CD45,IR_irl_CD45) -> PSI
PSI <- PSI[,row.names(tsne)]
dim(PSI)

data.frame(rowSums(!is.na(PSI)))
data.frame(rowSums(!is.na(PSI)), row.names = 1:nrow(PSI), name = row.names(PSI))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells.pdf")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(x = tsne$tSNE_1, y = tsne$tSNE_2, pch = 19, cex = .6, 
     col = rainbow(1, start =.1, alpha = .1), xlab = "tSNE_1", ylab = "tSNE_2")

points(x = tsne[colnames(PSI)[!is.na(PSI[1,])],]$tSNE_1, 
       y = tsne[colnames(PSI)[!is.na(PSI[1,])],]$tSNE_2, 
       col = rainbow(1), pch = 19, cex = .6)
points(x = tsne[colnames(PSI)[!is.na(PSI[2,])],]$tSNE_1, 
       y = tsne[colnames(PSI)[!is.na(PSI[2,])],]$tSNE_2, 
       col = rainbow(1,start = .2), pch = 19, cex = .6)
# Add legend to top right, outside plot region
legend("topright", inset=c(-0.2,.4), legend=c("SE1","SE2"), pch=19, cex = 1, bty = "n",
       col = c(rainbow(1),rainbow(1,start = .2)))
dev.off()



read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", sep = "\t", header = F, row.names = 1) -> celltype
colnames(celltype) <- "Cell type"
PSI[1:2,] -> psitu
psitu[is.na(psitu)] <- 0

heatmap(as.matrix(psitu), ColSideColors = as.character(celltype$V2), labRow = "")

identical(row.names(celltype), colnames(psitu))

library(pheatmap)

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells_heatmap.pdf")

pheatmap(mat = as.matrix(psitu), annotation_col = celltype, show_colnames = F, show_rownames = F)

test <- psitu[,colSums(psitu) > 0]

pheatmap(mat = as.matrix(test), show_colnames = F, show_rownames = F,
         annotation_col = subset.data.frame(celltype, row.names(celltype) %in% colnames(test)))
dev.off()


celltest <- subset.data.frame(celltype, row.names(celltype) %in% colnames(test))
celltest$YN <- 1
celltest <- celltest[order(celltest$`Cell type`),]

celltest -> test2
test2$cell <- row.names(test2)
as.data.frame(t(test)) -> test1
test1$cell <- row.names(test1)

test3 <- merge(test1, test2, by = "cell")
head(test3)
colnames(test3) <- c("cell","SE1","SE2","Cell","YN")

test3->test4

test4[test4 == 0] <- NA


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/top_AS_in_cells_boxplot.pdf")

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)

boxplot(SE1~Cell, test4, boxwex = 0.25, at = 1:6 - 0.2,col = "yellow", xaxt='n')
boxplot(SE2~Cell, test4, boxwex = 0.25, at = 1:6 + 0.2,col = "orange", add = T, xaxt='n')
axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, srt = 45,
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

table(test4[!is.na(test4$SE1),]$`Cell`)
table(test4[!is.na(test4$SE2),]$`Cell`)

axis(side = 3, at = 1:6 - .2, tick = F, labels = table(test4[!is.na(test4$SE1),]$`Cell`), line = -1)
axis(side = 3, at = 1:6 + .2, tick = F, labels = table(test4[!is.na(test4$SE2),]$`Cell`), line = -1)

legend("topright", inset=c(-0.2,.4), legend=c("SE1","SE2"), pch=19, cex = 1, bty = "n",
       col = c("yellow","orange"))

dev.off()







#### intron centric PSI for CD45 =================================================================================================================





load(file = "/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")

data.frame(do.call(rbind,strsplit(row.names(psi_sj_same_start_table), split="[:-]")), stringsAsFactors = F) -> sjinfo
sum(sjinfo$X1 == 1 & sjinfo$X2 >= 198638671 & sjinfo$X3 <= 198757283, na.rm=T)
# [1] 204

data.frame(do.call(rbind,strsplit(row.names(psi_sj_same_end_table), split="[:-]")), stringsAsFactors = F) -> sjinfo2
sum(sjinfo2$X1 == 1 & sjinfo2$X2 >= 198638671 & sjinfo2$X3 <= 198757283, na.rm=T)
# [1] 76


sj_CD45 <- rbind(psi_sj_same_start_table[sjinfo$X1 == 1 & sjinfo$X2 >= 198638671 & sjinfo$X3 <= 198757283,],
				psi_sj_same_end_table[sjinfo2$X1 == 1 & sjinfo2$X2 >= 198638671 & sjinfo2$X3 <= 198757283,])

sj_CD45[is.na(sj_CD45)] <- 0

summary(colSums(sj_CD45))
sum(colSums(sj_CD45) == 0)
# [1] 502

sj_CD45 <- sj_CD45[,row.names(celltype)] 


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_heatmap.pdf")

pheatmap(mat = as.matrix(sj_CD45), annotation_col = celltype, show_colnames = F, show_rownames = F)

pheatmap(mat = as.matrix(sj_CD45[apply(sj_CD45,1,sd)>=.1,]), annotation_col = celltype, show_colnames = F, show_rownames = F)

dev.off()



identical(colnames(sj_CD45), row.names(celltype))
# [1] TRUE

test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,
					SJ = sj_CD45[1,])

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_boxplot.pdf")

par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)

boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[1])
axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/All_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}

	
}

dev.off()




pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Meaningful_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0 & sd(test$SJ, na.rm = T) > 0 & sum(!is.na(test$SJ)) > 500){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}

	
}

dev.off()




pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Meaningful2_SJ_in_cells_boxplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
	print(i)
	test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
	test[test==0] <- NA

	if(sum(test$SJ, na.rm = T) > 0 & sd(test$SJ, na.rm = T) > 0 & sum(!is.na(test$SJ)) > 500 & median(test$SJ, na.rm = T) != 1){
		boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
	axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

	axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)

	}	
}

dev.off()




library(vioplot)
library(beanplot)
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Meaningful_SJ_in_cells_boxplot_and_vioplot.pdf")
par(mar=c(8.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(i in 1:nrow(sj_CD45)){
  print(i)
  test <- data.frame(cell = row.names(celltype), CellType = celltype$`Cell type`,SJ = sj_CD45[i,])
  test[test==0] <- NA

  if(sum(test$SJ, na.rm = T) > 0 & sd(test$SJ, na.rm = T) > 0 & sum(!is.na(test$SJ)) > 500){
    boxplot(SJ~CellType, test, boxwex = 0.25, col = "yellow", xaxt='n', main = row.names(sj_CD45)[i])
  axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(names(table(celltype$`Cell type`)),"\n", table(celltype$`Cell type`)))

  axis(side = 3, at = 1:6, tick = F, labels = table(test[test$SJ>0,]$`CellType`), line = -1)


  test2 <- test[!is.na(test$SJ),c("CellType", "SJ")]

  if(sum(lapply(split(test2, test2$CellType), function(x) sd(x$SJ))==0) == 0){
    print(i)
      vioplot(split(test2, test2$CellType)[[1]]$SJ,
        split(test2, test2$CellType)[[2]]$SJ,
        split(test2, test2$CellType)[[3]]$SJ,
        split(test2, test2$CellType)[[4]]$SJ,
        split(test2, test2$CellType)[[5]]$SJ,
        split(test2, test2$CellType)[[6]]$SJ, 
        names = F, main = row.names(sj_CD45)[i])
      axis(side = 1, at = 1:6,line = .3, tick = F, cex.axis=1, las = 2, 
     labels = paste(levels(test2$CellType),"\n", table(celltype$`Cell type`)))

  axis(side = 3, at = 1:6, labels = table(test[test$SJ>0,]$`CellType`), tick = F, las = 1)

}
}  
}

dev.off()





#### Seurat of intron-centric PSI ============================



setwd("/mnt/data5/BGI/UCB/tangchao/Seurat/CDC45")

library(Seurat)
apply(sj_CD45,1,sd) -> sj_CD45_sd
apply(sj_CD45,1,function(x) sum(x>0)) -> sj_CD45_mc
sj_CD45 <- sj_CD45[sj_CD45_sd!=0,]


row.names(sj_CD45) -> CD45_SJ
row.names(sj_CD45) <- paste("SJ",1:nrow(sj_CD45))

sj_CD45_2 <- t(t(sj_CD45)[!duplicated(t(sj_CD45)),])

ucb <- CreateSeuratObject( raw.data = sj_CD45_2, min.cells = 1, min.genes = 1, is.expr = 0, project="CD45")


# 3. Normalizing the data
#ucb <- NormalizeData(object = ucb, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Detection of variable genes across the single cells
pdf("Filter_genes.pdf",width = 8,height = 7)
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2, x.high.cutoff = 4, y.cutoff = -3)
length(x = ucb@var.genes)
## 2052 genes left
dev.off()

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
pdf("Select_pc.pdf",width = 8,height = 7)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:4)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-11

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = ucb, PCs = 1:4)
# good: PC1-15
PCElbowPlot(object = ucb)
dev.off()
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:4, resolution = 0.1, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:4, do.fast = TRUE)
pdf("tSNE.pdf",width = 8,height = 7)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()






#### Rtsne ========================================================================


library(Rtsne)
t(sj_CD45)[!duplicated(t(sj_CD45)),] -> test
tsne_out <- Rtsne(test, perplexity=560, theta = 1)
plot(tsne_out$Y, xlab = "t-SNE1", ylab = "t-SNE2", main = "UCB SE psi",pch = 20,
     cex = .5)

as.data.frame(tsne_out$Y) -> test2
test2$V3 <- celltype[row.names(test),]

plot(x = test2$V1, y = test2$V2, xlab = "t-SNE1", ylab = "t-SNE2", main = "CD45 Rtsne",pch = 20,
     cex = .8, col = as.numeric(test2$V3))
legend(x = 6.4, y = 1,legend = levels(test2$V3), pch = 20, col = 1:6, bty = "n", cex = .6)

sum(colSums(test>0)>=10)
test[,colSums(test>0)>=10] -> test3
test3[,apply(test3,2,function(x){sd(x[x>0])})>.3] ->test3
test3[!duplicated(test3),] -> test3
tsne_out2 <- Rtsne(test3, perplexity=300, theta = 1)
plot(tsne_out2$Y, xlab = "t-SNE1", ylab = "t-SNE2", main = "UCB SE psi",pch = 20,
     cex = .5)

as.data.frame(tsne_out2$Y) -> test4
test4$V3 <- celltype[row.names(test3),]

plot(x = test4$V1, y = test4$V2, xlab = "t-SNE1", ylab = "t-SNE2", main = "CD45 Rtsne",pch = 20,
     cex = .8, col = as.numeric(test4$V3))
legend(x = 5, y = 2,legend = levels(test4$V3), pch = 20, col = 1:6, bty = "n", cex = .6)










#### CD45 RABC =================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ1 = "1:198692374-198696711"
SJ2 = "1:198692374-198699563"
SJ3 = "1:198696910-198699563"
SJ4 = "1:198696910-198702386"
SJ5 = "1:198699705-198702386"
SJ6 = "1:198699705-198703297"
SJ7 = "1:198702531-198703297"
SJ8 = "1:198692374-198703297"
SJ9 = "1:198696910-198703297"
SJ10= "1:198692374-198702386"

te[.(c(SJ1,SJ2,SJ3,SJ4,SJ5,SJ6,SJ7,SJ8,SJ9,SJ10)),1:5]
as.data.frame(te[.(c(SJ1,SJ2,SJ3,SJ4,SJ5,SJ6,SJ7,SJ8,SJ9,SJ10))])[,-1] -> CD45_SJ

CD45_SJ[is.na(CD45_SJ)] <- 0
CD45_SJ[CD45_SJ<2] <- 0


CD45_R <- character()
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    R <- "CD45-"
  }else{
      if(CD45_SJ[8,i] > 0){
    R1 <- "RO"
  }else{
    R1 <- NULL
  }
  if(CD45_SJ[4,i] > 0){
    R2 <- "RAC"
  }else{
    R2 <- NULL
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    R3 <- "RBC"
  }else{
    R3 <- NULL
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    R4 <- "RAB"
  }else{
    R4 <- NULL
  }
  if(CD45_SJ[9,i] > 0){
    R5 <- "RA"
  }else{
    R5 <- NULL
  }
  if(CD45_SJ[10,i] > 0){
    R6 <- "RC"
  }else{
    R6 <-NULL
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    R7 <- "RB"
  }else{
    R7 <- NULL
  }
  if(CD45_SJ[2,i] < 0 & CD45_SJ[4,i] < 0 & CD45_SJ[6,i] < 0 & CD45_SJ[8,i] < 0 & CD45_SJ[4,i] < 0 & CD45_SJ[9,i] < 0 & CD45_SJ[10,i] < 0){
    R8 <- "RABC"
  }else{
    R8 <- NULL
  }
  else{
    R9 <- "Other"
  }
  R <- paste(R1,R2,R3,R4,R5,R6,R7,R8,R9,sep = "|")
  }
  CD45_R[i] <- R
}






CD45_R <- character()
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    R <- "CD45-"
  }else{
    R <- NULL
  }
  if(CD45_SJ[8,i] > 0){
    R1 <- "RO"
  }else{
    R1 <- NULL
  }
  if(CD45_SJ[4,i] > 0){
    R2 <- "RAC"
  }else{
    R2 <- NULL
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    R3 <- "RBC"
  }else{
    R3 <- NULL
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    R4 <- "RAB"
  }else{
    R4 <- NULL
  }
  if(CD45_SJ[9,i] > 0){
    R5 <- "RA"
  }else{
    R5 <- NULL
  }
  if(CD45_SJ[10,i] > 0){
    R6 <- "RC"
  }else{
    R6 <-NULL
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    R7 <- "RB"
  }else{
    R7 <- NULL
  }
  if(CD45_SJ[2,i] < 0 & CD45_SJ[4,i] < 0 & CD45_SJ[6,i] < 0 & CD45_SJ[8,i] < 0 & CD45_SJ[4,i] < 0 & CD45_SJ[9,i] < 0 & CD45_SJ[10,i] < 0){
    R8 <- "RABC"
  }else{
    R8 <- NULL
  }

  RE <- paste(c(R, R1,R2,R3,R4,R5,R6,R7,R8),collapse = "|")
  if(is.null(R) & is.null(R1) & is.null(R2) & is.null(R3) & is.null(R4) & is.null(R5) & is.null(R6) & is.null(R7) & is.null(R8)){
    RE <- "Other"
  }
  
  CD45_R[i] <- RE
}


CD45_R_tab <- CD45_SJ[1:9,]
row.names(CD45_R_tab) <- c("CD45-","RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_tab[CD45_R_tab>0] <- 0
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    CD45_R_tab["CD45-",i] <- 1
  }
  if(CD45_SJ[8,i] > 0){
    CD45_R_tab["RO",i] <- 1
  }
  if(CD45_SJ[4,i] > 0){
    CD45_R_tab["RAC",i] <- 1
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    CD45_R_tab["RBC",i] <- 1
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    CD45_R_tab["RAB",i] <- 1
  }
  if(CD45_SJ[9,i] > 0){
    CD45_R_tab["RA",i] <- 1
  }
  if(CD45_SJ[10,i] > 0){
    CD45_R_tab["RC",i] <- 1
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    CD45_R_tab["RB",i] <- 1
  }
  if(sum(CD45_SJ[,i]) > 0 & CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0 & CD45_SJ[8,i] == 0 & CD45_SJ[9,i] == 0 & CD45_SJ[10,i] == 0){
    CD45_R_tab["RABC",i] <- 1
  }
}



load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

colnames(CD45_R_tab) <- substr(colnames(CD45_R_tab),1,10)
CD45_R_tab[,row.names(tsne)] -> CD45_R_tab

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(tsne, pch = 20, cex = 0)
# CD45-
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[1,] > 0],], pch = 20, col = rainbow(1, start = .1, alpha = .4))
# RO
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[2,] > 0],], pch = 20, col = rainbow(1, start = .2, alpha = .4))
# RAC
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[3,] > 0],], pch = 20, col = rainbow(1, start = .3, alpha = .4))
# RBC
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[4,] > 0],], pch = 20, col = rainbow(1, start = .4, alpha = .4))
# RAB
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[5,] > 0],], pch = 20, col = rainbow(1, start = .5, alpha = .4))
# RA
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[6,] > 0],], pch = 20, col = rainbow(1, start = .6, alpha = .4))
# RC
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[7,] > 0],], pch = 20, col = rainbow(1, start = .7, alpha = .4))
# RB 
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[8,] > 0],], pch = 20, col = rainbow(1, start = .8, alpha = .4))
# RABC
points(tsne[colnames(CD45_R_tab)[CD45_R_tab[9,] > 0],], pch = 20, col = rainbow(1, start = .9, alpha = .4))
legend(x = 24,y = 36, legend = row.names(CD45_R_tab), pch = 20, cex = 1.2, col = rainbow(9,start = .1,end = .9), bty = "n")
rowSums(CD45_R_tab)

par(mfrow = c(3,3), mar = c(0.5,0.5,1,0.5))
for(i in 1:9){
  plot(tsne, pch = 20, cex = 0, main = row.names(CD45_R_tab)[i], xaxt="n",  yaxt="n")
  points(tsne[colnames(CD45_R_tab)[CD45_R_tab[i,] > 0],], pch = 20, col = rainbow(1, start = i/10))
}


table(colSums(CD45_R_tab))
x = barplot(table(colSums(CD45_R_tab)), ylim = c(0,3000))
text(x = x, y = table(colSums(CD45_R_tab))+100, labels = table(colSums(CD45_R_tab)))

CD45_R_tab[,colSums(CD45_R_tab)==3]

rowSums(CD45_R_tab[,colSums(CD45_R_tab)==2])
CD45_R_tab[,colSums(CD45_R_tab)==2][,which(CD45_R_tab[,colSums(CD45_R_tab)==2][2,]==0)]



#### Use the juction reads information also =======================================================================================


CD45_R_r_tab <- CD45_SJ[1:9,]
row.names(CD45_R_r_tab) <- c("CD45-","RO","RAC","RBC","RAB","RA","RC","RB","RABC")
CD45_R_r_tab[CD45_R_r_tab>0] <- 0
for (i in 1:ncol(CD45_SJ)) {
  if(sum(CD45_SJ[,i]) == 0){
    CD45_R_r_tab["CD45-",i] <- 1
  }
  if(CD45_SJ[8,i] > 0){
    CD45_R_r_tab["RO",i] <- CD45_SJ[8,i]
  }
  if(CD45_SJ[4,i] > 0){
    CD45_R_r_tab["RAC",i] <- CD45_SJ[4,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0){
    CD45_R_r_tab["RBC",i] <- CD45_SJ[2,i]
  }
  if(CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RAB",i] <- CD45_SJ[6,i]
  }
  if(CD45_SJ[9,i] > 0){
    CD45_R_r_tab["RA",i] <- CD45_SJ[9,i]
  }
  if(CD45_SJ[10,i] > 0){
    CD45_R_r_tab["RC",i] <- CD45_SJ[10,i]
  }
  if(CD45_SJ[2,i] > 0 & CD45_SJ[6,i] > 0){
    CD45_R_r_tab["RB",i] <- mean(CD45_SJ[2,i], CD45_SJ[6,i])
  }
  if(sum(CD45_SJ[,i]) > 0 & CD45_SJ[2,i] == 0 & CD45_SJ[4,i] == 0 & CD45_SJ[6,i] == 0 & CD45_SJ[8,i] == 0 & CD45_SJ[9,i] == 0 & CD45_SJ[10,i] == 0){
    CD45_R_r_tab["RABC",i] <- mean(CD45_SJ[c(1,3,5,7),i])
}
}


colnames(CD45_R_r_tab) <- substr(colnames(CD45_R_r_tab),1,10)
CD45_R_r_tab[,row.names(tsne)] -> CD45_R_r_tab

CD45_R_r_tab[,1:5]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==3]
CD45_R_r_tab[,colSums(CD45_R_r_tab>0)==2][,1:10]

apply(CD45_R_r_tab, 2, function(x) round(x/sum(x),2)) -> tang

tang[,colSums(tang>0)==2] -> test

hist(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), breaks = 20, main = "", xlab = "Fraction of minor subset")
plot(density(as.numeric(apply(test,2,function(x){ 1-max(x/sum(x)) })), adjust = 0.2 ), main = "", xlab = "Fraction of minor subset")
abline(v = 0.05, col = 2)


tang[tang < 0.05] <- 0
tang[,colSums(tang>0)==2][,1:10]
tang[tang>0] <- 1

table(colSums(tang))
x = barplot(table(colSums(tang)), ylim = c(0,3400))
text(x = x, y = table(colSums(tang))+200, labels = table(colSums(tang)))
par(mfrow = c(3,3), mar = c(0.5,0.5,1,0.5))
for(i in 1:9){
  plot(tsne, pch = 20, cex = 1, main = row.names(tang)[i], xaxt="n",  yaxt="n", col = gray(0.9))
  points(tsne[colnames(tang)[tang[i,] > 0],], pch = 20, cex = .2, col = rainbow(1, start = i/10))
}




dim(CD45_R_r_tab)
dim(tang)
identical(colnames(CD45_R_r_tab),colnames(tang))

CD45_R_r_tab -> chao
chao[tang==0] <- 0
apply(chao,1,summary)

chao[,colSums(chao>0)==2][,1:10]


ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log(as.numeric(chao[9,]+1))), alpha = 1, size = .8,
             position = position_jitter(.2, .2))+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)

ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log(as.numeric(chao[8,]+1))), alpha = 1, size = .8,
             position = position_jitter(.2, .2))+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)




ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao[9,]+1))), alpha = 1, size = .8,
             position = position_jitter(.2, .2))+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)

ggplot(data = tsne, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao[8,]+1))), alpha = 1, size = .8,
             position = position_jitter(.2, .2))+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao)[8])+
  theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
          panel.border = element_blank())

t(t(chao)[order(t(chao)[,1]),]) -> chao1
f1 <- ggplot(data = tsne[colnames(chao1),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao1[1,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle("R-")+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,2]),]) -> chao2
f2 <- ggplot(data = tsne[colnames(chao2),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao2[2,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao2)[2])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,3]),]) -> chao3
f3 <- ggplot(data = tsne[colnames(chao3),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao3[3,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao3)[3])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,4]),]) -> chao4
f4 <- ggplot(data = tsne[colnames(chao4),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao4[4,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao4)[4])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,5]),]) -> chao5
f5 <- ggplot(data = tsne[colnames(chao5),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao5[5,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao5)[5])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,6]),]) -> chao6
f6 <- ggplot(data = tsne[colnames(chao6),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao6[6,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao6)[6])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,7]),]) -> chao7
f7 <- ggplot(data = tsne[colnames(chao7),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao7[7,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao7)[7])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,8]),]) -> chao8
f8 <- ggplot(data = tsne[colnames(chao8),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao8[8,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao8)[8])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

t(t(chao)[order(t(chao)[,9]),]) -> chao9
f9 <- ggplot(data = tsne[colnames(chao9),], aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = log10(as.numeric(chao9[9,]+1))), size = .2)+
  scale_colour_gradient(high = 'red',low = 'Grey', guide = F)+
  ggtitle(row.names(chao9)[9])+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())




library(cowplot)
plot_grid(f1, f2, f3, f4, f5, f6, f7, f8, f9, ncol = 3)






