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
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- JackStraw(object = ucb, num.replicate = 100)
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
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_1.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()



ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:4, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
ucb <- RunTSNE(object = ucb, dims.use = 1:4, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_4s.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()




ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:7, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
ucb <- RunTSNE(object = ucb, dims.use = 1:7, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()






#### Second =====================================================================================================================================



# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 3.0, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- JackStraw(object = ucb, num.replicate = 100)
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
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_2.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()










#### Third ======================================================================================================================================




# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 3.0, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- JackStraw(object = ucb, num.replicate = 100)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_3.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()




















#### Fourth =====================================================================================================================================




# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 2.0, x.high.cutoff = 6, y.cutoff = -3)
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
ucb <- JackStraw(object = ucb, num.replicate = 100, prop.freq = 0.05, num.pc = 20)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_4.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()












#### Fifth =====================================================================================================================================


sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])


te_sub[sj_CD8,row.names(Cell_type[grepl(pattern="T", Cell_type$CellType),])] -> cd8_sub
dim(cd8_sub)
[1]  142 2257

library(Seurat)
library(dplyr)
library(Matrix)

psi_sub <- cd8_sub[!duplicated(cd8_sub), ]
psi_sub[is.na(cd8_sub)] <- 0
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- JackStraw(object = ucb, num.replicate = 100, prop.freq = 0.01, num.pc = 20)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:4, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:4, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_CD8_for_CD48_2.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()





ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:4, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
ucb <- RunTSNE(object = ucb, dims.use = 1:7, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()







































#################################################################################################################################################
#### Find markers for correlation with gene expression
#################################################################################################################################################


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
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- JackStraw(object = ucb, num.replicate = 100)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

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


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()

dev.off()
samplegroup = sub(".[0-9]+$", "", names(ucb@ident))
names(samplegroup) = names(ucb@ident)
samplegroup = as.factor(samplegroup)
ucb@ident = as.factor(samplegroup)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5)





# . Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers <- ucb.markers[ucb.markers$p_val_adj<0.05,]
ucb.markers %>% group_by(cluster) %>% top_n(10) -> top10
ucb.markers %>% group_by(cluster) %>% top_n(2) -> top2
#write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.txt",sep = "\t",quote = F)
#pdf("New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)
#dev.off()

#pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.pdf",width = 20,height = 20)
#FeaturePlot(object = ucb, features.plot = top2$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
#dev.off()

table(ucb.markers$cluster)
 0  1  2  3  4  5  6  7
16 17 10 14  8  8 26 13

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster0.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==0]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster1.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==1]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster2.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==2]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster3.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==3]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster4.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==4]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster5.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==5]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster6.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==6]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_New_markers_cluster7.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==7]), cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_Seurat_final.Rda")




load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")

tsne_ge <- as.data.frame(Combine@dr$tsne@cell.embeddings)
dim(tsne_ge)
tsne_ge$cluster <- 0
names(ucb@ident)[ucb@ident==7]
tsne_ge[names(ucb@ident)[ucb@ident==7],]$cluster <- 1
tsne_ge$cluster <- as.factor(tsne_ge$cluster)

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_cluster7_marker_CRTAM_in_gene_expression.pdf",width = 20,height = 20)
FeaturePlot(object = Combine, features.plot = "CRTAM", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T, pch.use = 20)
ggplot(tsne_ge, aes(x = tSNE_1, y = tSNE_2, colour = cluster))+
  geom_point(size = 2)
dev.off()



#### CD8 in tSNE

psi_sub[sj_CD8, ]

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_CD8SJ_in_tSNE.pdf")
for (i in 1:nrow(psi_sub[sj_CD8, ])){
  tmp_data <- data.frame(tsne, PSI = as.numeric(psi_sub[sj_CD8, ][i,row.names(tsne)]))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(psi_sub[sj_CD8, ])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


#### CD8 in bulk

load('/mnt/data1/projects/UCB/data_ucb/BP_ucb_bulk_SJ/psi_same_start_end_SJ.RData')
load("/mnt/data1/projects/UCB/results_ucb/tangchao/bulk_gene_and_psi/RData/bulk_gene_Seurat.RData") ## just for the sample information
library(data.table)

psi_sj_same_start[sj_CD8,] -> bulk_CD8_psi
psi_sj_same_start[paste("chr",sj_CD8,sep=""),] -> bulk_CD8_psi
dim(bulk_CD8_psi)
[1] 142  59
data_plot[grep("CD[48]-",data_plot$cellname),]
row.names(data_plot[grep("CD[48]-",data_plot$cellname),])
#[1] "S000RD15" "S0018A13"

bulk_CD8_psi[,row.names(data_plot[grep("CD[48]-",data_plot$cellname),])] -> bulk_cd48_CD8_psi
colnames(bulk_cd48_CD8_psi) <- paste(colnames(bulk_cd48_CD8_psi),c("CD8","CD4","CD4","CD8"),sep="_")
library(RColorBrewer)
pheatmap(bulk_cd48_CD8_psi, show_rownames=F, color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100))

bulk_cd48_CD8_psi[rowSums(is.na(bulk_cd48_CD8_psi))>0,]



bulk_na_sj <- substr(row.names(bulk_cd48_CD8_psi[rowSums(is.na(bulk_cd48_CD8_psi))>0,]),4,100)

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_bulk_CD8_sep_SJ_in_SC_tSNE.pdf")
for (i in 1:nrow(psi_sub[bulk_na_sj, ])){
  tmp_data <- data.frame(tsne, PSI = as.numeric(psi_sub[bulk_na_sj, ][i,row.names(tsne)]))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(psi_sub[bulk_na_sj, ])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()

save(bulk_na_sj, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/bulk_cd8_sep_sj.RData")




#### Marker SJs' Host gene expression

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
gene_exp <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/UCB_rsem_TPM_mat.txt", sep = "\t", header=T, row.names=1)

cluster3_host_gene <- unique(sj[sj$sj %in% with(ucb.markers,gene[cluster==3]),"SJ_start_exon"])

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_Cluster3_hostgene_in_tSNE.pdf")
for(i in 1:length(cluster3_host_gene)){
  tmp_data <- data.frame(tsne, PSI = log2(as.numeric(gene_exp[cluster3_host_gene,][i,row.names(tsne)])+1))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  #tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(gene_exp[cluster3_host_gene,])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


cluster4_host_gene <- unique(sj[sj$sj %in% with(ucb.markers,gene[cluster==4]),"SJ_start_exon"])
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_Cluster4_hostgene_in_tSNE.pdf")
for(i in 1:length(cluster4_host_gene)){
  tmp_data <- data.frame(tsne, PSI = log2(as.numeric(gene_exp[cluster4_host_gene,][i,row.names(tsne)])+1))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  #tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(gene_exp[cluster4_host_gene,])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


cluster5_host_gene <- unique(sj[sj$sj %in% with(ucb.markers,gene[cluster==5]),"SJ_start_exon"])
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_Cluster5_hostgene_in_tSNE.pdf")
for(i in 1:length(cluster5_host_gene)){
  tmp_data <- data.frame(tsne, PSI = log2(as.numeric(gene_exp[cluster5_host_gene,][i,row.names(tsne)])+1))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  #tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(gene_exp[cluster5_host_gene,])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


cluster6_host_gene <- unique(sj[sj$sj %in% with(ucb.markers,gene[cluster==6]),"SJ_start_exon"])
cluster6_host_gene <- unlist(strsplit(cluster6_host_gene,"\\|"))
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_Cluster6_hostgene_in_tSNE.pdf")
for(i in 1:length(cluster6_host_gene)){
  tmp_data <- data.frame(tsne, PSI = log2(as.numeric(gene_exp[cluster6_host_gene,][i,row.names(tsne)])+1))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  #tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(gene_exp[cluster6_host_gene,])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()


cluster7_host_gene <- unique(sj[sj$sj %in% with(ucb.markers,gene[cluster==7]),"SJ_start_exon"])
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_tSNE_7s_2_Cluster7_hostgene_in_tSNE.pdf")
for(i in 1:length(cluster7_host_gene)){
  tmp_data <- data.frame(tsne, PSI = log2(as.numeric(gene_exp[cluster7_host_gene,][i,row.names(tsne)])+1))
  tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
  #tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
  #tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
  print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(colour = PSI), size = 2)+
      scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
      ggtitle(row.names(gene_exp[cluster7_host_gene,])[i])+
      theme(panel.background = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            axis.text= element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank())
)
}#statements
dev.off()

save.image(file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/Tcell_Seurat.RData")





















