#### BPSC for Cell Type Specific AS events
library(BPSC)
#################################################################################################################################################
#### SE #########################################################################################################################################
#################################################################################################################################################


load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
PSI <- SE_psi[, 13:ncol(SE_psi)]
row.names(PSI) <- SE_psi$loci

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

PSI <- PSI[, row.names(Cell_type)]
dim(PSI)
# [1] 17670  3039

rs <- rowSums(!is.na(PSI))
sum(rs == 0)
#[1] 53

PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0

group=c(rep(1,sum(Cell_type$CellType == "B Cells")),rep(2,sum(Cell_type$CellType != "B Cells")))

controlIds = colnames(PSI_tu)[which(group == 2)]

design = model.matrix(~group)
coef = 2
# Run BPglm for differential expression analysis
library(doParallel)
registerDoParallel(cores=16)

res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)

hist(res$PVAL, breaks = 20)

# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 1379
bpglm.DE.ids -> BCel_Sep



table(Cell_type$CellType)
group <- rep(c(1,2,1,1,1,1), table(Cell_type$CellType))
controlIds = colnames(PSI_tu)[which(group == 1)]
design = model.matrix(~group)
coef = 2
res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
hist(res$PVAL, breaks = 20)
# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 1153
bpglm.DE.ids -> CD4_Sep


group <- rep(c(1,1,2,1,1,1), table(Cell_type$CellType))
controlIds = colnames(PSI_tu)[which(group == 1)]
design = model.matrix(~group)
coef = 2
res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
hist(res$PVAL, breaks = 20)
# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 793
bpglm.DE.ids -> CD8_Sep


group <- rep(c(1,1,1,2,1,1), table(Cell_type$CellType))
controlIds = colnames(PSI_tu)[which(group == 1)]
design = model.matrix(~group)
coef = 2
res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
hist(res$PVAL, breaks = 20)
# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 152
bpglm.DE.ids -> MK_Sep


group <- rep(c(1,1,1,1,2,1), table(Cell_type$CellType))
controlIds = colnames(PSI_tu)[which(group == 1)]
design = model.matrix(~group)
coef = 2
res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
hist(res$PVAL, breaks = 20)
# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 376
bpglm.DE.ids -> Mono_Sep


group <- rep(c(1,1,1,1,1,2), table(Cell_type$CellType))
controlIds = colnames(PSI_tu)[which(group == 1)]
design = model.matrix(~group)
coef = 2
res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)
hist(res$PVAL, breaks = 20)
# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 295
bpglm.DE.ids -> NK_Sep



lapply(list(BCel_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep),length)

length(Reduce(union, list(BCel_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep)))



c(BCel_Sep, CD4_Sep, CD8_Sep, MK_Sep, Mono_Sep, NK_Sep) -> juncs
PSI[juncs, row.names(Cell_type)] -> te_sub
dim(te_sub)


te_sub[is.na(te_sub)] <- -0.2

cells <- c("B Cells", "CD4+ T Cells", "CD8+ T Cells", "Megakaryocytes", "Monocytes", "NK cells")
sj_anno <- matrix(rep(cells, c(length(BCel_Sep),length(CD4_Sep),length(CD8_Sep),length(MK_Sep),length(Mono_Sep),length(NK_Sep))), ncol = 1)
rownames(sj_anno) <- row.names(te_sub)
#sj_anno <- data.frame(sj_anno[juncs,])
colnames(sj_anno) <- "Cell"
sj_anno <- as.data.frame(sj_anno)
sj_anno$Cell <- factor(sj_anno$Cell, levels = cells)

dim(sj_anno)
# [1] 1215    1
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
  Cell = c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")
)
names(ann_colors[[2]]) <- names(table(Cell_type$CellType))
names(ann_colors[[3]]) <- c(names(table(Cell_type$CellType)))

cumsum(as.numeric(table(sj_anno$Cell)))[1:5]


library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/All_Cell_specific_SE_heatmap_BPSC.pdf")
pheatmap(te_sub, annotation_col = Cell_type, color = colors, breaks = breaksList, 
  annotation_names_row = FALSE, gaps_row = cumsum(as.numeric(table(sj_anno$Cell)))[1:5],
  gaps_col = cumsum(as.numeric(table(Cell_type$CellType)))[1:5],
  cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F)
dev.off()







PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0
PSI_tu <- PSI_tu[, row.names(Cell_type[Cell_type$CellType %in% c("B Cells", "CD8+ T Cells"),])]

group = c(rep(1,sum(Cell_type$CellType == "B Cells")),rep(2,sum(Cell_type$CellType == "CD8+ T Cells")))

controlIds = colnames(PSI_tu)[which(group == 1)]

design = model.matrix(~group)
coef = 2
# Run BPglm for differential expression analysis
library(doParallel)
registerDoParallel(cores=16)

res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)

hist(res$PVAL, breaks = 20)

# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 1313
bpglm.DE.ids

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Cell_specific_SE_heatmap_BPSC_CD8.Vs.BCel.pdf")
pheatmap(PSI_tu[bpglm.DE.ids,], cluster_cols = F, cluster_rows = T, show_colnames = F, show_rownames = F)
dev.off()





PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 0
PSI_tu <- PSI_tu[, row.names(Cell_type[Cell_type$CellType %in% c("CD4+ T Cells", "CD8+ T Cells"),])]

group = c(rep(1,sum(Cell_type$CellType == "CD4+ T Cells")),rep(2,sum(Cell_type$CellType == "CD8+ T Cells")))

controlIds = colnames(PSI_tu)[which(group == 1)]

design = model.matrix(~group)
coef = 2
# Run BPglm for differential expression analysis
library(doParallel)
registerDoParallel(cores=16)

res = BPglm(data = as.matrix(PSI_tu), controlIds = controlIds, design = design, coef = coef, estIntPar = TRUE, useParallel = TRUE)

hist(res$PVAL, breaks = 20)

# Summarize the resutls
ss = summary(res)

fdr = p.adjust(res$PVAL, method = "BH")
bpglm.DE.ids = row.names(PSI_tu)[which(fdr <= 0.01)]
length(bpglm.DE.ids)# Print the indices of the true DE genes:
# [1] 414
bpglm.DE.ids


PSI_tu[bpglm.DE.ids,] -> plot_psi

Cell_type[colnames(plot_psi),] -> Cell_type_plot
Cell_type_plot$CellType <- as.factor(as.character(Cell_type_plot$CellType))

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Cell_specific_SE_heatmap_BPSC_CD8.Vs.CD4.pdf")
pheatmap(plot_psi, cluster_cols = F, cluster_rows = T, show_colnames = F, show_rownames = F, annotation_col = Cell_type_plot)
dev.off()






