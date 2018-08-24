load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

psi_sj_same_start_table <- psi_sj_same_start_table[,row.names(Cell_type)]

psi_sj_same_start_table[is.na(psi_sj_same_start_table)] <- 0
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sj_same_start_table, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = 4.5, y.cutoff = -5)
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
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/tSNE_1.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()




#### Second =====================================================================================================================================

ucb <- CreateSeuratObject(raw.data = psi_sj_same_start_table, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

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
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/tSNE_2.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()






#### Third ======================================================================================================================================

ucb <- CreateSeuratObject(raw.data = psi_sj_same_start_table, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

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
# good: PC1-11

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:6, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:6, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/tSNE_3.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()










