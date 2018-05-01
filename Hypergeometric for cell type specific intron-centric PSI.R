#### Hypergeometric for cell type specific intron-centric AS SJ:
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

te <- psi_sj_same_start_table[, row.names(Cell_type)]
dim(te)
# [1] 460327   3039
rm("psi_sj_same_start_table")
gc()

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

te <- as.data.frame(te)
FastRemoveMissingValues(te)


Bcell_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

Bcell_Sep <- as.data.frame(Bcell_Sep)
Bcell_Sep$Cell_p <- Bcell_Sep$k/Bcell_Sep$N
Bcell_Sep$Bg_p <- Bcell_Sep$D/(Bcell_Sep$D + Bcell_Sep$n)
Bcell_Sep$Other_p <- (Bcell_Sep$D - Bcell_Sep$k)/(Bcell_Sep$D + Bcell_Sep$n - Bcell_Sep$N)
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1, ])
# [1] 1463    9
dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
# [1] 1003    9

Bcell_Sep <- Bcell_Sep_formore


load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
te <- psi_sj_same_start_table[, row.names(Cell_type)]
dim(te)
# [1] 460327   3039
rm("psi_sj_same_start_table")
gc()
te <- as.data.frame(te)


Bcell_Sep <- t(apply(te, 1, function(x){
  x <- x[!is.na(x)]
  k = sum(x[names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",])] > 0.15)#
  #k = sum(x[which(Cell_type$CellType=="B Cells")] >= 0.15)

  D = sum(x >= 0.15)
  n = length(x) - D
  N = sum(names(x) %in% row.names(Cell_type[Cell_type$CellType=="B Cells",]))
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
Bcell_Sep -> Bcell_Sep2

hist(Bcell_Sep$pval)


te <- as.data.frame(te)
FastRemoveMissingValues(te)

CD4_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

CD4_Sep <- as.data.frame(CD4_Sep)
CD4_Sep$Cell_p <- CD4_Sep$k/CD4_Sep$N
CD4_Sep$Bg_p <- CD4_Sep$D/(CD4_Sep$D + CD4_Sep$n)
CD4_Sep$Other_p <- (CD4_Sep$D - CD4_Sep$k)/(CD4_Sep$D + CD4_Sep$n - CD4_Sep$N)

head(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ],20)
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1, ])
# [1] 539   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
# [1] 84  9


CD8_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

CD8_Sep <- as.data.frame(CD8_Sep)
CD8_Sep$Cell_p <- CD8_Sep$k/CD8_Sep$N
CD8_Sep$Bg_p <- CD8_Sep$D/(CD8_Sep$D + CD8_Sep$n)
CD8_Sep$Other_p <- (CD8_Sep$D - CD8_Sep$k)/(CD8_Sep$D + CD8_Sep$n - CD8_Sep$N)

head(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ],20)
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1, ])
# [1] 1293    9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
# [1] 142   9


MK_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Megakaryocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

MK_Sep <- as.data.frame(MK_Sep)
MK_Sep$Cell_p <- MK_Sep$k/MK_Sep$N
MK_Sep$Bg_p <- MK_Sep$D/(MK_Sep$D + MK_Sep$n)
MK_Sep$Other_p <- (MK_Sep$D - MK_Sep$k)/(MK_Sep$D + MK_Sep$n - MK_Sep$N)
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1, ])
# [1] 2729    9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
# [1] 2694    9



Mono_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="Monocytes")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

Mono_Sep <- as.data.frame(Mono_Sep)
Mono_Sep$Cell_p <- Mono_Sep$k/Mono_Sep$N
Mono_Sep$Bg_p <- Mono_Sep$D/(Mono_Sep$D + Mono_Sep$n)
Mono_Sep$Other_p <- (Mono_Sep$D - Mono_Sep$k)/(Mono_Sep$D + Mono_Sep$n - Mono_Sep$N)
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1, ])
# [1] 1701    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
# [1] 1647    9



NK_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="NK cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

NK_Sep <- as.data.frame(NK_Sep)
NK_Sep$Cell_p <- NK_Sep$k/NK_Sep$N
NK_Sep$Bg_p <- NK_Sep$D/(NK_Sep$D + NK_Sep$n)
NK_Sep$Other_p <- (NK_Sep$D - NK_Sep$k)/(NK_Sep$D + NK_Sep$n - NK_Sep$N)
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1, ])
# [1] 1143    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])
# [1] 947   9


save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sep_PSI.RData")


dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ]) # [1] 1003   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ]) # [1] 84 9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ]) # [1] 142  9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ]) # [1] 2694    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ]) # [1] 1647    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ]) # [1] 947   9

dim(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Bg_p, ]) # [1] 851   9
dim(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Bg_p, ]) # [1] 0 9
dim(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Bg_p, ]) # [1] 71  9
dim(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Bg_p, ]) # [1] 2694    9
dim(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Bg_p, ]) # [1] 1646    9
dim(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Bg_p, ]) # [1] 935   9



sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])


length(c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK))
# [1] 6517
length(unique(c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)))
# [1] 5916

unique(c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)) -> juncs

psi_sj_same_start_table[juncs,row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 5916 3039
save(te_sub, Cell_type, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_intron_centric_PSI.RData")

te_sub[is.na(te_sub)] <- -0.2

sj_anno <- matrix(rep(names(table(Cell_type$CellType)), c(length(sj_Bcell),length(sj_CD4),length(sj_CD8),length(sj_MK),length(sj_Mono),length(sj_NK))), ncol = 1)
row.names(sj_anno) <- c(sj_Bcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)
sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
dim(sj_anno)
# [1] 5916    1
identical(row.names(sj_anno), row.names(te_sub))
# [1] TRUE


library(RColorBrewer)
library(pheatmap)
breaksList = seq(-0.2,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(14)[4:14])
#colors = c(rep(brewer.pal(3,"Greens")[1], 4), brewer.pal(,"Reds")[4:9])

brewer.pal(6,"BrBG")
ann_colors = list(
  Individual = c(UCB1 = "#7570B3", UCB3 = "#E7298A", UCB4 = "#66A61E"),
  CellType = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E"),
  Cell = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- names(table(Cell_type$CellType))



pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_PSI_heatmap2.pdf")
pheatmap(te_sub,annotation_col = Cell_type, color = colors, breaks = breaksList, annotation_row = sj_anno, 
  annotation_colors = ann_colors, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()





#### Add T Cell Specific PSI ====

Tcell_Sep <- t(apply(te, 1, function(x){
  k = sum(x[which(Cell_type$CellType=="CD4+ T Cells" | Cell_type$CellType=="CD8+ T Cells")] >= 0.15)
  D = sum(x >= 0.15)
  n = length(x) - D
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

Tcell_Sep <- as.data.frame(Tcell_Sep)
Tcell_Sep$Cell_p <- Tcell_Sep$k/Tcell_Sep$N
Tcell_Sep$Bg_p <- Tcell_Sep$D/(Tcell_Sep$D + Tcell_Sep$n)
Tcell_Sep$Other_p <- (Tcell_Sep$D - Tcell_Sep$k)/(Tcell_Sep$D + Tcell_Sep$n - Tcell_Sep$N)

head(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ],20)
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1, ])
# [1] 1132    9
dim(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])
# [1] 524   9

save(Bcell_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep, Tcell_Sep, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_sep_PSI.RData")


sj_Tcell <- row.names(Tcell_Sep[Tcell_Sep$adj_pval < 0.01 & Tcell_Sep$Cell_p > .1 & Tcell_Sep$Cell_p > 2*Tcell_Sep$Other_p, ])


length(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK))
# [1] 7041
length(unique(c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)))
# [1] 6336

c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK) -> juncs

psi_sj_same_start_table[juncs,row.names(Cell_type)] -> te_sub
dim(te_sub)
# [1] 7041 3039
save(te_sub, Cell_type, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_intron_centric_PSI2.RData")

te_sub[is.na(te_sub)] <- -0.2

cells <- c("B Cells", "T Cells", "CD4+ T Cells", "CD8+ T Cells", "Megakaryocytes", "Monocytes", "NK cells")
sj_anno <- matrix(rep(cells, c(length(sj_Bcell),length(sj_Tcell),length(sj_CD4),length(sj_CD8),length(sj_MK),length(sj_Mono),length(sj_NK))), ncol = 1)
row.names(sj_anno) <- c(sj_Bcell, sj_Tcell, sj_CD4, sj_CD8, sj_MK, sj_Mono, sj_NK)
#sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
sj_anno <- as.data.frame(sj_anno)
sj_anno$Cell <- factor(sj_anno$Cell, levels = cells)
dim(sj_anno)
# [1] 7041    1
identical(row.names(sj_anno), row.names(te_sub))
# [1] TRUE


library(RColorBrewer)
library(pheatmap)
breaksList = seq(-0.2,1,0.1)
library(RColorBrewer)
colors = c(brewer.pal(3,"Greens")[1], colorRampPalette(brewer.pal(n = 9, name = "Reds"))(11))
#colors = c(rep(brewer.pal(3,"Greens")[1], 4), brewer.pal(,"Reds")[4:9])

brewer.pal(6,"BrBG")
ann_colors = list(
  Individual = c(UCB1 = "#7570B3", UCB3 = "#E7298A", UCB4 = "#66A61E"),
  CellType = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E"),
  Cell = c("#8C510A","#7570B3","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- c(names(table(Cell_type$CellType))[1],"T Cells", names(table(Cell_type$CellType))[-1])

cumsum(as.numeric(table(sj_anno$Cell)))[1:6]


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_PSI_heatmap3.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno$Cell)))[1:6],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()





#### Add Annotation information ====

load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj10.RData")
sj$classification <- NA
sj[sj$annotation==1,]$classification <- "Annotated"
sj[sj$novel=="Y",]$classification <- "Novel"
sj[sj$annotation==0 & sj$SCSpe == "N",]$classification <- "Unannotated"
sj$classification <- factor(sj$classification, levels = c("Annotated","Unannotated","Novel"))
table(sj$classification)

dim(te_sub)
# [1] 7041 3039
sum(row.names(te_sub) %in% sj$sj)
# [1] 7041

table(with(sj, classification[sj %in% row.names(te_sub)]))
#  Annotated Unannotated       Novel
#       6112         144          80



table(with(sj, classification[sj %in% sj_Bcell]))
#  Annotated Unannotated       Novel
#        958          39           6
table(with(sj, classification[sj %in% sj_Tcell]))
#  Annotated Unannotated       Novel
#        509          15           0
table(with(sj, classification[sj %in% sj_CD4]))
#  Annotated Unannotated       Novel
#         81           3           0
table(with(sj, classification[sj %in% sj_CD8]))
#  Annotated Unannotated       Novel
#        140           2           0
table(with(sj, classification[sj %in% sj_MK]))
#  Annotated Unannotated       Novel
#       2640          49           5
table(with(sj, classification[sj %in% sj_Mono]))
#  Annotated Unannotated       Novel
#       1549          29          69
table(with(sj, classification[sj %in% sj_NK]))
#  Annotated Unannotated       Novel
#        921          20           6

table(with(sj, classification[sj %in% sj_Bcell]))/sum(table(with(sj, classification[sj %in% sj_Bcell])))
#   Annotated Unannotated       Novel
# 0.955134596 0.038883350 0.005982054



sum(row.names(te_sub) %in% with(sj, sj[classification != "Annotated"]))


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_PSI_heatmap.pdf")
pheatmap(te_sub[row.names(te_sub) %in% with(sj, sj[classification != "Annotated"]),], 
  annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(no_NA).RData")
te <- te_table
rm("te_table")
gc()
#### Lib size normalization

te_colsum <- colSums(te)

for(i in 1:ncol(te)){
  te[,i] <- (te[,i]/te_colsum[i])*mean(te_colsum)
}
save(te, file = "/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")
te <- as.matrix(te)

te_sub2 <- te[row.names(te_sub[row.names(te_sub) %in% with(sj, sj[classification != "Annotated"]),]), row.names(Cell_type)]

te_sub2[is.na(te_sub2)] <- 0

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_heatmap.pdf")
pheatmap(log10(te_sub2+1), 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


test <- data.frame(R = te_sub2[1,], I = as.factor(Cell_type$Individual))
summary(aov(R ~ I, test))[[1]][[1,"Pr(>F)"]]

aovp <- vector()
for(i in 1:nrow(te_sub2)){
  test <- data.frame(R = te_sub2[i,], I = Cell_type$Individual)
  aovp[i] <- summary(aov(R ~ I, test))[[1]][[1,"Pr(>F)"]]
}

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Individual_specific_Unannotated_SJ_Reads_pvalue.pdf")
hist(aovp)
dev.off()



pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap.pdf")
pheatmap(log10(te_sub2+1)[aovp < 0.01, row.names(Cell_type[order(Cell_type$Individual),])], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()



aovp <- vector()
for(i in 1:nrow(te_sub2)){
  test <- data.frame(R = te_sub2[i,], I = Cell_type$Individual)
  test <- test[test$R>0, ]
  if(length(unique(test$I))>1){
       aovp[i] <- summary(aov(R ~ I, test))[[1]][[1,"Pr(>F)"]]
    }else{
       aovp[i] <- 0
    }
}

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Individual_specific_Unannotated_SJ_Reads_pvalue2.pdf")
hist(aovp)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap2.pdf")
pheatmap(log10(te_sub2+1)[aovp < 0.01, row.names(Cell_type[order(Cell_type$Individual),])], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


tabn <- vector()
for(i in 1:nrow(te_sub2)){
  test <- data.frame(R = as.character(te_sub2[i,]), I = Cell_type$Individual, stringsAsFactors = F)
  test <- test[test$R > 0, ]
       tabn[i] <- length(unique(test$I))
}


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap3.pdf")
pheatmap(log10(te_sub2+1)[tabn == 1, row.names(Cell_type[order(Cell_type$Individual),])], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap4.pdf")
pheatmap(log10(te_sub2+1)[tabn == 3, ], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap5.pdf")
pheatmap(log10(te_sub2+1)[tabn != 3, ], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_Unannotated_SJ_Reads_Indivadual_specific_heatmap6.pdf")
pheatmap(log10(te_sub2+1)[tabn >= 2, ], 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()

dim(te_sub2[tabn >= 2,])
#[1]  210 3039

#### sashimi plot ====================================================================================================================

novel_sj_reads <- rowSums(te_sub2[tabn >= 2,])
novel_sj_maxreads <- as.data.frame(apply(te_sub2[tabn >= 2,], 1, function(x) colnames(te_sub2)[which.max(x)]))
colnames(novel_sj_maxreads) <- "Cell"
novel_sj_maxreads$SJ <- row.names(te_sub2[tabn >= 2,])
save(novel_sj_maxreads, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.RData")

#### sashimi2
sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, tsv = "/mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/", outDir = "~/", min = 10, ann_height = 4){
    if(dir.exists(outDir) == FALSE) dir.create(outDir)
    Coordinates <- data.frame(do.call(rbind,strsplit(junction, split = "[:-]")),stringsAsFactors=F)
    Coordinates$X2 <- abs(as.numeric(Coordinates$X2) - left_expand)
    Coordinates$X3 <- abs(as.numeric(Coordinates$X3) + right_expand)
    window <- paste(Coordinates$X1, ":", Coordinates$X2, "-", Coordinates$X3, sep = "")
    if(length(cell) == 1){
        Dir_out <- paste(outDir, cell, "_", junction, sep = "")
        work <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi.sh ", window, cell, Dir_out, min, ann_height, "')")
        eval(parse(text = work))
    }else{
        Dir_out <- paste(outDir, junction, sep = "")
        tmpdir <- tempdir()
        work1 <- paste("system('", "cat ", paste(tsv,cell,".tsv",sep = "", collapse = " "), " > ", tmpdir,"/tmp.tsv", "')", sep = "")
        work2 <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi2.sh ", window, paste(tmpdir,"/tmp.tsv", sep = ""), Dir_out, min, ann_height, "')")
        work3 <- paste("system('", "rm",paste(tmpdir,"/tmp.tsv", "')", sep = ""))
        eval(parse(text = work1))
        eval(parse(text = work2))
        eval(parse(text = work3))
    }
}


novel_sj_maxreads$Novel <- 1
for (i in 1:nrow(novel_sj_maxreads)) {
  if(novel_sj_maxreads[i,]$SJ %in% sj[sj$classification=="Unannotated",]$sj){
    novel_sj_maxreads[i,]$Novel <- 0
  }
}
table(novel_sj_maxreads$Novel)
#   0   1
# 127  83

save(novel_sj_maxreads, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.RData")

for(i in 1:nrow(novel_sj_maxreads)){
  print(i)
  if(novel_sj_maxreads[i,"Novel"] == 1){
    sashimi(junction = novel_sj_maxreads[i,]$SJ, cell = novel_sj_maxreads[i,]$Cell, outDir = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/sashimi/Novel/", ann_height = 20)
  }else{
    sashimi(junction = novel_sj_maxreads[i,]$SJ, cell = novel_sj_maxreads[i,]$Cell, outDir = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/sashimi/Unannotated/", ann_height = 20)
  }
  
}





load("/mnt/data5/BGI/UCB/Rdata/CD45_samples.RData")

novel_CD45 <- te[c("1:198694141-198696711", "1:198694141-198699563", "1:198729766-198731616"), row.names(Cell_type)]

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/Unannotated_CD45_SJ.pdf")
pheatmap(log10(novel_CD45+1), 
  annotation_col = Cell_type, annotation_names_row = FALSE, 
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()







