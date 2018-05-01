
#### Hypergeometric for cell type specific SJ:
load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")

library(RColorBrewer)
require(data.table)  # v1.6.6
require(gdata) 
FastRemoveMissingValues = function(DT) {
  # either of the following for loops
  
  # by name :
  for (j in names(DT))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

FastRemoveMissingValues(te)

te <- data.frame(as.data.frame(te)[,-1], row.names = te[[1]])

FastRemoveMissingValues = function(DT) {
  # either of the following for loops
  
  # by name :
  #for (j in names(DT))
  #  set(DT,which(DT[[j]]==1),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==1),j,0)
}

colnames(te) <- substr(colnames(te), 1, 10)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

te <- te[, row.names(Cell_type)]

rowSums(te) -> te_rowsum
te <- te[te_rowsum > 0, ]
dim(te)
# [1] 1784983   3039


Bcell_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="B Cells")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="B Cells")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


CD4_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="CD4+ T Cells")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="CD4+ T Cells")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


CD8_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="CD8+ T Cells")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="CD8+ T Cells")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


MK_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="Megakaryocytes")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="Megakaryocytes")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


Mono_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="Monocytes")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="Monocytes")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


NK_Sep <- t(apply(te, 1, function(x){
	k = sum(x[which(Cell_type$CellType=="NK cells")] > 0)
	D = sum(x>0)
	n = length(x) - sum(x>0)
	N = sum(Cell_type$CellType=="NK cells")
	pval = phyper(k, D, n, N, lower.tail=FALSE)
	if(k==0) {
    adj_pval <- pval
  	} else {
    adj_pval <- pval * k
  	}
  	enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))


Tcell_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")] > 0)
  D = sum(x>0)
  n = length(x) - sum(x>0)
  N = sum(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")
  pval = phyper(k, D, n, N, lower.tail=FALSE)
  if(k==0) {
    adj_pval <- pval
    } else {
    adj_pval <- pval * k
    }
    enrichment <- c(as.integer(k), D, n, N, pval, adj_pval); names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval") 
  return(enrichment)
}))



Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)

CD4_Sep <- as.data.frame(CD4_Sep)
CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)

CD8_Sep <- as.data.frame(CD8_Sep)
CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)

MK_Sep <- as.data.frame(MK_Sep)
MK_Sep$Cell_p <- MK_Sep$k/MK_Sep$N
MK_Sep$Bg_p <- MK_Sep$D/(MK_Sep$D + MK_Sep$n)
MK_Sep$Other_p <- (MK_Sep$D - MK_Sep$k)/(MK_Sep$D + MK_Sep$n - MK_Sep$N)

Mono_Sep <- as.data.frame(Mono_Sep)
Mono_Sep$Cell_p <- Mono_Sep$k/Mono_Sep$N
Mono_Sep$Bg_p <- Mono_Sep$D/(Mono_Sep$D + Mono_Sep$n)
Mono_Sep$Other_p <- (Mono_Sep$D - Mono_Sep$k)/(Mono_Sep$D + Mono_Sep$n - Mono_Sep$N)

NK_Sep <- as.data.frame(NK_Sep)
NK_Sep$Cell_p <- NK_Sep$k/NK_Sep$N
NK_Sep$Bg_p <- NK_Sep$D/(NK_Sep$D + NK_Sep$n)
NK_Sep$Other_p <- (NK_Sep$D - NK_Sep$k)/(NK_Sep$D + NK_Sep$n - NK_Sep$N)

Tcell_Sep <- as.data.frame(Tcell_Sep)
Tcell_Sep$Cell_p <- Tcell_Sep$k/Tcell_Sep$N
Tcell_Sep$Bg_p <- Tcell_Sep$D/(Tcell_Sep$D + Tcell_Sep$n)
Tcell_Sep$Other_p <- (Tcell_Sep$D - Tcell_Sep$k)/(Tcell_Sep$D + Tcell_Sep$n - Tcell_Sep$N)

save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, Tcell_Sep, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sep.RData")


head(Bcell_Sep[Bcell_Sep[,"adj_pval"] == 0 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"n"]!=0 & Bcell_Sep[,"k"]>10,], 40)

head(Bcell_Sep[Bcell_Sep[,"adj_pval"] == 0 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10,], 40)

head(Bcell_Sep[Bcell_Sep[,"adj_pval"] == 0 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10 & Bcell_Sep[,"k"]<607,], 40)

sum(Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10)


sum(Bcell_Sep[,"adj_pval"] == 0 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10)
# [1] 122
sum(CD4_Sep[,"adj_pval"] == 0 & CD4_Sep[,"D"]>0 & CD4_Sep[,"k"]>10)
# [1] 81
sum(CD8_Sep[,"adj_pval"] == 0 & CD8_Sep[,"D"]>0 & CD8_Sep[,"k"]>10)
# [1] 66


sum(Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10 & Bcell_Sep$Cell_p > 4*Bcell_Sep$Bg_p)
# [1] 934
sum(CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep[,"k"]>10 & CD4_Sep$Cell_p > 1.5*CD4_Sep$Bg_p)
# [1] 286
sum(CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep[,"k"]>10 & CD8_Sep$Cell_p > 2*CD8_Sep$Bg_p)
# [1] 621
sum(MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep[,"k"]>10 & MK_Sep$Cell_p > 2*MK_Sep$Bg_p)
# [1] 196
sum(Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep[,"k"]>10 & Mono_Sep$Cell_p > 4*Mono_Sep$Bg_p)
# [1] 614
sum(NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep[,"k"]>10 & NK_Sep$Cell_p > 4*NK_Sep$Bg_p)
# [1] 411

sum(CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep[,"k"]>10 & CD8_Sep$Cell_p > 4*CD8_Sep$Other_p)
# [1] 210
sum(CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep[,"k"]>10 & CD4_Sep$Cell_p > 4*CD4_Sep$Other_p)
# [1] 277

sum(Tcell_Sep[,"adj_pval"] < 0.01 & Tcell_Sep[,"D"]>0 & Tcell_Sep[,"k"]>10 & Tcell_Sep$Cell_p > 1.5*Tcell_Sep$Bg_p)
# [1] 0
sum(Tcell_Sep[,"adj_pval"] < 0.01 & Tcell_Sep[,"D"]>0 & Tcell_Sep[,"k"]>10 & Tcell_Sep$Cell_p > Tcell_Sep$Bg_p)
# [1] 3111
sum(Tcell_Sep[,"adj_pval"] < 0.01 & Tcell_Sep[,"D"]>0 & Tcell_Sep[,"k"]>1000 & Tcell_Sep$Cell_p > Tcell_Sep$Bg_p)
# [1] 360

sum(Tcell_Sep[,"adj_pval"] < 0.01 & Tcell_Sep[,"D"]>0 & Tcell_Sep[,"k"]>100 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p)
# [1] 720

sum(Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p) # 1129
sum(CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p) # 53
sum(CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p) # 140
sum(MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p) # 3041
sum(Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p) # 2400
sum(NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p) # 1079

sum(Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 4*Bcell_Sep$Bg_p) # 244
sum(CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 4*CD4_Sep$Bg_p) # 0
sum(CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 4*CD8_Sep$Bg_p) # 3
sum(MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 4*MK_Sep$Bg_p) # 2662
sum(Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 4*Mono_Sep$Bg_p) # 1563
sum(NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 4*NK_Sep$Bg_p) # 469



Bcell_j <- row.names(Bcell_Sep[Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep[,"k"]>10 & Bcell_Sep$Cell_p > 4*Bcell_Sep$Bg_p,])
CD4_j <- row.names(CD4_Sep[CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep[,"k"]>10 & CD4_Sep$Cell_p > 4*CD4_Sep$Other_p,])
CD8_j <- row.names(CD8_Sep[CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep[,"k"]>10 & CD8_Sep$Cell_p > 4*CD8_Sep$Other_p,])
MK_j <- row.names(MK_Sep[MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep[,"k"]>10 & MK_Sep$Cell_p > 2*MK_Sep$Bg_p,])
Mono_j <- row.names(Mono_Sep[Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep[,"k"]>10 & Mono_Sep$Cell_p > 4*Mono_Sep$Bg_p,])
NK_j <- row.names(NK_Sep[NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep[,"k"]>10 & NK_Sep$Cell_p > 4*NK_Sep$Bg_p,])

Tcell_j <- row.names(Tcell_Sep[Tcell_Sep[,"adj_pval"] < 0.01 & Tcell_Sep[,"D"]>0 & Tcell_Sep[,"k"]>100 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p,])

length(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j))
# [1] 2642
length(unique(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j)))
# [1] 2612

unique(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j))[unique(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j)) %in% row.names(te)] -> juncs

#### Lib size normalization

te_colsum <- colSums(te)

for(i in 1:ncol(te)){
  te[,i] <- (te[,i]/te_colsum[i])*mean(te_colsum)
}

te[juncs,row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 2612 3039
save(te_sub, Cell_type, sj, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_Normalized_for_Seurat.RData")

mat = log10(as.matrix(te[juncs,row.names(Cell_type)])+1)

library(RColorBrewer)
library(pheatmap)
breaksList = seq(0,5,0.5)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(20)[12:20])

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_Normalized_breaks.pdf")
pheatmap(mat,annotation_col = Cell_type, color = colors, breaks = breaksList,
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()





length(c(Bcell_j,CD8_j,MK_j,Mono_j,NK_j))
# [1] 2776
length(unique(c(Bcell_j,CD8_j,MK_j,Mono_j,NK_j)))
# [1] 2730

unique(c(Bcell_j,CD8_j,MK_j,Mono_j,NK_j))[unique(c(Bcell_j,CD8_j,MK_j,Mono_j,NK_j)) %in% row.names(te)] -> juncs

te[juncs,] -> reads_sub

Cell_type_sub <- Cell_type[Cell_type$CellType != "CD4+ T Cells", ]

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/test3.pdf")
pheatmap(log10(as.matrix(reads_sub[,row.names(Cell_type_sub)])+1),annotation_col = Cell_type_sub, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)



dim(reads_sub)
[1] 2730 3039
sum(row.names(reads_sub) %in% sj$sj)
[1] 2729
sum(row.names(reads_sub) %in% sj[sj$classification == "Novel",]$sj)
[1] 384
sum(row.names(reads_sub) %in% sj[sj$classification == "Annotated",]$sj)
[1] 2150
sum(row.names(reads_sub) %in% sj[sj$classification == "Unannotated",]$sj)
[1] 195
sum(row.names(reads_sub) %in% sj[sj$classification != "Annotated",]$sj)
[1] 579

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Annotated_Cell_Specific_SJ.pdf")
pheatmap(log10(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification == "Annotated",]$sj,row.names(Cell_type_sub)])+1),annotation_col = Cell_type_sub, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Novel_Cell_Specific_SJ.pdf")
pheatmap(log10(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification == "Novel",]$sj,row.names(Cell_type_sub)])+1),annotation_col = Cell_type_sub, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Unannotated_Cell_Specific_SJ.pdf")
pheatmap(log10(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification != "Annotated",]$sj,row.names(Cell_type_sub)])+1),annotation_col = Cell_type_sub, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()





unique(c(Bcell_j,Tcell_j,MK_j,Mono_j,NK_j))[unique(c(Bcell_j,Tcell_j,MK_j,Mono_j,NK_j)) %in% row.names(te)] -> juncs


te[juncs,] -> reads_sub

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/test4.pdf")
pheatmap(log10(as.matrix(reads_sub[,row.names(Cell_type)])+1),annotation_col = Cell_type, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()




#### Lib size normalization

te_colsum <- apply(te,2,sum)

for(i in 1:ncol(te)){
  te[,i] <- (te[,i]/te_colsum[i])*mean(te_colsum)
}


te[juncs,] -> reads_sub

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_SJ_Normalized.pdf")
pheatmap(log2(as.matrix(reads_sub[,row.names(Cell_type)])+1),annotation_col = Cell_type, 
    cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


library(RColorBrewer)
breaksList = seq(0,5,0.5)

colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(20)[12:20])

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_SJ_Normalized_breaks.pdf")
pheatmap(log10(as.matrix(reads_sub[,row.names(Cell_type)])+1),annotation_col = Cell_type, 
  color = colors, breaks = breaksList,
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()

save(reads_sub, Cell_type, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_SJ_Normalized_breaks.RData")



#### Annotated or Unannotated
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Annotated_Cell_Specific_SJ_Normalized_breaks.pdf")
pheatmap(log10(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification == "Annotated",]$sj,row.names(Cell_type)])+1),
  annotation_col = Cell_type, cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F, 
  color = colors, breaks = breaksList)
dev.off()



mat <- log2(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification != "Annotated",]$sj,row.names(Cell_type)])+1)
breaksList = seq(0,max(mat),1)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(20)[8:20])
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Unannotated_Cell_Specific_SJ_Normalized_breaks.pdf")
pheatmap(mat, annotation_col = Cell_type, cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F, 
  color = colors, breaks = breaksList)
dev.off()

mat <- log2(as.matrix(reads_sub[row.names(reads_sub) %in% sj[sj$classification == "Novel",]$sj,row.names(Cell_type)])+1)
breaksList = seq(0,max(mat)+1,1)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(20)[8:20])
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Novel_Cell_Specific_SJ_Normalized_breaks.pdf")
pheatmap(mat, annotation_col = Cell_type, cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F, 
  color = colors, breaks = breaksList)
dev.off()




#### For Seurat: -------

sum(Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p) # 1129
sum(CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p) # 53
sum(CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p) # 140
sum(MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p) # 3041
sum(Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p) # 2400
sum(NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p) # 1079

Bcell_j <- row.names(Bcell_Sep[Bcell_Sep[,"adj_pval"] < 0.01 & Bcell_Sep[,"D"]>0 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p,])
CD4_j <- row.names(CD4_Sep[CD4_Sep[,"adj_pval"] < 0.01 & CD4_Sep[,"D"]>0 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p,])
CD8_j <- row.names(CD8_Sep[CD8_Sep[,"adj_pval"] < 0.01 & CD8_Sep[,"D"]>0 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p,])
MK_j <- row.names(MK_Sep[MK_Sep[,"adj_pval"] < 0.01 & MK_Sep[,"D"]>0 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p,])
Mono_j <- row.names(Mono_Sep[Mono_Sep[,"adj_pval"] < 0.01 & Mono_Sep[,"D"]>0 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p,])
NK_j <- row.names(NK_Sep[NK_Sep[,"adj_pval"] < 0.01 & NK_Sep[,"D"]>0 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p,])

unique(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j))[unique(c(Bcell_j,CD4_j,CD8_j,MK_j,Mono_j,NK_j)) %in% row.names(te)] -> juncs
te[juncs,] -> reads_sub
dim(reads_sub)
# [1] 7253 3039

save(reads_sub, Cell_type, sj, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/For_Seurat.RData")






load("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_specific_SJ_Normalized_breaks.RData")

library(Seurat)
library(dplyr)
library(Matrix)

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = reads_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb, normalization.method = "LogNormalize")

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 4.5, y.cutoff = -5)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-11

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_Cell_specific_SJ.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()






load("/Users/tangchao/CloudStation/tangchao/project/UCB/All_Cell_specific_SJ_Normalized_for_Seurat.RData")

library(Seurat)
library(dplyr)
library(Matrix)

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = te_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb, normalization.method = "LogNormalize")

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 4.5, y.cutoff = -5)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-20

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_All_Cell_specific_SJ.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()



load("/Users/tangchao/CloudStation/tangchao/project/UCB/For_Seurat.RData")

library(Seurat)
library(dplyr)
library(Matrix)

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = reads_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb, normalization.method = "LogNormalize")

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2, x.high.cutoff = 4.5, y.cutoff = -5)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-20

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_All_Cell_specific_SJ2.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()


#### Unannotated ====

dim(reads_sub)
# [1] 7253 3039
sum(row.names(reads_sub) %in% sj[sj$classification != "Annotated", ]$sj)
reads_sub_unannotated <- reads_sub[row.names(reads_sub) %in% sj[sj$classification != "Annotated", ]$sj,]
dim(reads_sub_unannotated)
# [1]  857 3039

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = reads_sub_unannotated, min.cells = 1, min.genes = 1, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb, normalization.method = "LogNormalize")

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 7, y.cutoff = -5)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-18

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE, num.pc = 18)
JackStrawPlot(object = ucb, PCs = 1:18)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:18, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_All_Cell_specific_Unannotated_SJ.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()



#### Unannotated SJ intron-centric psi ====
load("/Users/tangchao/CloudStation/tangchao/project/UCB/psi_list_left_right_table.RData")
dim(reads_sub_unannotated)
# [1]  857 3039
dim(psi_sj_same_start_table)
# [1] 460327   3574
sum(row.names(reads_sub_unannotated) %in% row.names(psi_sj_same_start_table))
# [1] 731

psi <- psi_sj_same_start_table[row.names(reads_sub_unannotated)[row.names(reads_sub_unannotated) %in% row.names(psi_sj_same_start_table)], colnames(reads_sub_unannotated)]
dim(psi)
# [1]  731 3039
psi[is.na(psi)] <- 0

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi, min.cells = 1, min.genes = 1, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 7, y.cutoff = -8)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-18

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100, do.print = FALSE, num.pc = 20)
JackStrawPlot(object = ucb, PCs = 1:18)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:13, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:13, do.fast = F)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_All_Cell_specific_Unannotated_SJ_PSI.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()






load("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sep.RData")
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)
dim(te_sub)
sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])

table(with(sj, classification[sj %in% sj_Bcell]))
#  Annotated Unannotated       Novel
#       1042          60          22
table(with(sj, classification[sj %in% sj_Tcell]))
#  Annotated Unannotated       Novel
#        279           9          28
table(with(sj, classification[sj %in% sj_CD4]))
#  Annotated Unannotated       Novel
#         50           3           0
table(with(sj, classification[sj %in% sj_CD8]))
#  Annotated Unannotated       Novel
#        138           2           0
table(with(sj, classification[sj %in% sj_MK]))
#  Annotated Unannotated       Novel
#       2895          74          72
table(with(sj, classification[sj %in% sj_Mono]))
#  Annotated Unannotated       Novel
#       1807          45         547
table(with(sj, classification[sj %in% sj_NK]))
#  Annotated Unannotated       Novel
#       1012          30          37




load(file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_sep_PSI.RData")
psi_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
psi_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
psi_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
psi_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
psi_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
psi_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
psi_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])


sum(sj_Bcell %in% sj[sj$classification == "Novel",]$sj & !sj_Bcell %in% psi_Bcell)
[1] 16
sum(sj_Tcell %in% sj[sj$classification == "Novel",]$sj & !sj_Tcell %in% psi_Tcell)
[1] 28
sum(sj_MK %in% sj[sj$classification == "Novel",]$sj & !sj_MK %in% psi_MK)
[1] 67
sum(sj_Mono %in% sj[sj$classification == "Novel",]$sj & !sj_Mono %in% psi_Mono)
[1] 480
sum(sj_NK %in% sj[sj$classification == "Novel",]$sj & !sj_NK %in% psi_NK)
[1] 33


load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")

sj_Bcell_sub <- sj_Bcell[sj_Bcell %in% sj[sj$classification == "Novel",]$sj & !sj_Bcell %in% psi_Bcell]
sj_Tcell_sub <- sj_Tcell[sj_Tcell %in% sj[sj$classification == "Novel",]$sj & !sj_Tcell %in% psi_Tcell]
sj_MK_sub <- sj_MK[sj_MK %in% sj[sj$classification == "Novel",]$sj & !sj_MK %in% psi_MK]
sj_Mono_sub <- sj_Mono[sj_Mono %in% sj[sj$classification == "Novel",]$sj & !sj_Mono %in% psi_Mono]
sj_NK_sub <- sj_NK[sj_NK %in% sj[sj$classification == "Novel",]$sj & !sj_NK %in% psi_NK]

te_sub <- psi_sj_same_start_table[c(sj_Bcell_sub,sj_Tcell_sub,sj_MK_sub,sj_Mono_sub,sj_NK_sub)[c(sj_Bcell_sub,sj_Tcell_sub,sj_MK_sub,sj_Mono_sub,sj_NK_sub) %in% row.names(psi_sj_same_start_table)],]

te_sub[is.na(te_sub)] <- -0.2


Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_PSI_specific_SJ_PSI_heatmap.pdf")
pheatmap(te_sub[, row.names(Cell_type)], annotation_col = Cell_type, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



rowsum <- rowSums(te_sub[, row.names(Cell_type)] >= 0)      
tail(sort(rowsum),6)

te_sub[names(tail(sort(rowsum),4)), row.names(Cell_type)]

Bcell_Sep[names(tail(sort(rowsum),6)),1:5]
                       k  D    n   N         pval
22:22906554-22915824  56 61 2978 607 3.519789e-36
7:142788548-142801040  0  0 3039 607 0.000000e+00
5:77348709-77582954    4 46 2993 607 9.674885e-01
6:32554794-32758953    0  0 3039 607 0.000000e+00
6:32554794-32814423    0  0 3039 607 0.000000e+00
7:142792799-142802502  0  0 3039 607 0.000000e+00
Mono_Sep[names(tail(sort(rowsum),6)),1:5]
                      k  D    n  N      pval
22:22906554-22915824  0 61 2978 52 0.6547281
7:142788548-142801040 0  0 3039 52 0.0000000
5:77348709-77582954   0 46 2993 52 0.5506271
6:32554794-32758953   0  0 3039 52 0.0000000
6:32554794-32814423   0  0 3039 52 0.0000000
7:142792799-142802502 0  0 3039 52 0.0000000

apply(psi_sj_same_start_table[names(tail(sort(rowsum),6)),],1,function(x) sum(x>=0,na.rm=T))
 22:22906554-22915824 7:142788548-142801040   5:77348709-77582954
                  294                   436                   514
  6:32554794-32758953   6:32554794-32814423 7:142792799-142802502
                  754                   754                  1753
apply(psi_sj_same_start_table[names(tail(sort(rowsum),6)),],1,function(x) sum(x>0,na.rm=T))
 22:22906554-22915824 7:142788548-142801040   5:77348709-77582954
                   79                     8                    60
  6:32554794-32758953   6:32554794-32814423 7:142792799-142802502
                    2                    10                   247 
                     

apply(psi_sj_same_start_table[names(tail(sort(rowsum),6)), row.names(Cell_type[Cell_type$CellType=="B Cells",])],1,function(x) sum(x>=0,na.rm=T))
 22:22906554-22915824 7:142788548-142801040   5:77348709-77582954
                  209                    15                    73
  6:32554794-32758953   6:32554794-32814423 7:142792799-142802502
                  492                   492                   100
 
 apply(psi_sj_same_start_table[names(tail(sort(rowsum),6)), row.names(Cell_type[Cell_type$CellType=="B Cells",])],1,function(x) sum(x>0,na.rm=T))
 22:22906554-22915824 7:142788548-142801040   5:77348709-77582954
                   65                     0                     4
  6:32554794-32758953   6:32554794-32814423 7:142792799-142802502
                    2                     7                     4
                    
sums <- rowSums(te_sub[,row.names(Cell_type[Cell_type$CellType=="B Cells",])] >= 0)      

tail(sort(sums),2)
6:32554794-32758953 6:32554794-32814423
                492                 492
names(tail(sort(sums),2))
[1] "6:32554794-32758953" "6:32554794-32814423"



apply(te[c("6:32554794-32758953", "6:32554794-32814423"),],1,function(x) sum(x>=0,na.rm=T))

apply(te[names(tail(sort(rowsum),6)),],1,function(x) sum(x>=10,na.rm=T))
 22:22906554-22915824 7:142788548-142801040   5:77348709-77582954
                   71                     7                    51
  6:32554794-32758953   6:32554794-32814423 7:142792799-142802502
                    2                     7                   224

apply(te[names(tail(sort(rowsum),6)), row.names(Cell_type[Cell_type$CellType=="B Cells",])],1,function(x) sum(x>=0,na.rm=T))


