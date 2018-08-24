load(file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_intron_centric_PSI2.RData")
load("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sep_PSI.RData")
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])

union(sj_CD4, sj_CD8) -> juncs

te_sub[juncs,row.names(Cell_type[grepl(pattern="T", Cell_type$CellType),])] -> tcel_sub
dim(tcel_sub)
[1]  226 2257

library(Seurat)
library(dplyr)
library(Matrix)

psi_sub <- tcel_sub[!duplicated(tcel_sub), ]
psi_sub[is.na(tcel_sub)] <- 0
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/FindVariableGenes.pdf")
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 6, y.cutoff = -5)
length(x = ucb@var.genes)
dev.off()

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/PCA.pdf")
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-11
dev.off()

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100)
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/JackStraw.pdf")
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10
dev.off()

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:7, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:7, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/Tcell_tSNE_7s_2.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()

dev.off()


# 10. Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
library(dplyr)
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
ucb.markers %>% group_by(cluster) %>% top_n(2, avg_logFC) -> top2
write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers.txt",sep = "\t",quote = F)

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_FeaturePlot_top2.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = top2$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_FeaturePlot_top10.pdf",width = 20,height = 60)
FeaturePlot(object = ucb, features.plot = top10$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_RidgePlot_top2.pdf",width = 20,height = 20)
print(RidgePlot(object = ucb, features.plot = top2$gene))
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_RidgePlot_top10.pdf",width = 20,height = 60)
print(RidgePlot(object = ucb, features.plot = top10$gene))
dev.off()


pdf('/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_VlnPlot_top2.pdf',width = 20,height = 20)
print(VlnPlot(object = ucb, features.plot = top2$gene))
dev.off()
pdf('/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell/New_markers_VlnPlot_top10.pdf',width = 20,height = 60)
print(VlnPlot(object = ucb, features.plot = top10$gene))
dev.off()


