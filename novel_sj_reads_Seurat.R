load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(no_NA).RData")
Basic_stat <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Basic_stat.txt", header = T, row.names = 1, stringsAsFactors=F)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Basic_stat <- Basic_stat[row.names(Cell_type),]
te_table <- te_table[, row.names(Cell_type)]

for(i in 1:ncol(te_table)){
	te_table[i] <- te_table[i]/Basic_stat[i, "Mapped_reads"]*mean(Basic_stat[, "Mapped_reads"])
}

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
length(with(sj,sj[novel=="Y"]))
#[1] 613530

novel_sj <- te_table[with(sj,sj[novel=="Y"]), ]
dim(novel_sj)
#[1] 613530   3039
novel_sj_rs <- rowSums(novel_sj)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#     0.00      1.81      4.81     79.58     24.72 136768.09
novel_sj_rs_c <- rowSums(novel_sj>0)
summary(novel_sj_rs_c)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.000    1.000    1.000    1.924    1.000 2165.000
sum(novel_sj_rs_c>=10)
#[1] 12412
novel_sj2 <- novel_sj[novel_sj_rs_c>=10,]
novel_sj2_sd <- apply(novel_sj2, 1, function(x) sd(x[x>0]))
summary(novel_sj2_sd)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.179    0.686    1.358   29.561    4.518 4250.377

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = novel_sj2, min.cells = 10, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_Detect_variable_genes.pdf")
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 5, y.cutoff = 0.5)
length(x = ucb@var.genes)
dev.off()
# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_PCA.pdf")
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
dev.off()
ucb = ProjectPCA(object = ucb, do.print = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_PCHeatmap.pdf")
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()
# good: PC1-11

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100)

pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_JackStrawPlot.pdf")
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10
dev.off()

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.2, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_TSNEPlot.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


samplegroup = sub(".[0-9]+$", "", names(ucb@ident))
names(samplegroup) = names(ucb@ident)
samplegroup = as.factor(samplegroup)
ucb@ident = as.factor(samplegroup)
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_tSNE_sample.pdf")
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5)
dev.off()


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- merge(tsne, Cell_type, by = 0, all.x = T)
#tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
library(RColorBrewer)
cols_cells <- brewer.pal(8, "Set1")[c(2:5,7,8)]
#types <- c("CD8", "CD4", "NK", "B", "Mono", "Mega")
tsne$CellType <- factor(tsne$CellType, levels=c("CD8+ T Cells", "CD4+ T Cells", "NK cells", "B Cells", "Monocytes", "Megakaryocytes"))

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()+
  scale_color_manual(values = cols_cells)




pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/novel_sj_seurat_TSNEPlot_gene_x0_y0.5.pdf")


ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- merge(tsne, Cell_type, by = 0, all.x = T)
#tsne$CellType <- as.factor(Cell_type$CellType)
tsne$CellType <- factor(tsne$CellType, levels=c("CD8+ T Cells", "CD4+ T Cells", "NK cells", "B Cells", "Monocytes", "Megakaryocytes"))

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()+
  scale_color_manual(values = cols_cells)+
  ggtitle("PC1-10")


ucb <- RunTSNE(object = ucb, dims.use = 1:12, do.fast = TRUE)

tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- merge(tsne, Cell_type, by = 0, all.x = T)
#tsne$CellType <- as.factor(Cell_type$CellType)
tsne$CellType <- factor(tsne$CellType, levels=c("CD8+ T Cells", "CD4+ T Cells", "NK cells", "B Cells", "Monocytes", "Megakaryocytes"))

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()+
  scale_color_manual(values = cols_cells)+
  ggtitle("PC1-12")



ucb <- RunTSNE(object = ucb, dims.use = 1:15, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- merge(tsne, Cell_type, by = 0, all.x = T)
#tsne$CellType <- as.factor(Cell_type$CellType)
tsne$CellType <- factor(tsne$CellType, levels=c("CD8+ T Cells", "CD4+ T Cells", "NK cells", "B Cells", "Monocytes", "Megakaryocytes"))

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()+
  scale_color_manual(values = cols_cells)+
  ggtitle("PC1-15")


ucb <- RunTSNE(object = ucb, dims.use = 1:20, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- merge(tsne, Cell_type, by = 0, all.x = T)
#tsne$CellType <- as.factor(Cell_type$CellType)
tsne$CellType <- factor(tsne$CellType, levels=c("CD8+ T Cells", "CD4+ T Cells", "NK cells", "B Cells", "Monocytes", "Megakaryocytes"))

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()+
  scale_color_manual(values = cols_cells)+
  ggtitle("PC1-20")

dev.off()



ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)

# 10. Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
library(dplyr)
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
ucb.markers %>% group_by(cluster) %>% top_n(2, avg_logFC) -> top2
write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers.txt",sep = "\t",quote = F)
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers_FeaturePlot_top2.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = top2$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers_FeaturePlot_top10.pdf",width = 20,height = 60)
FeaturePlot(object = ucb, features.plot = top10$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers_RidgePlot_top2.pdf",width = 20,height = 20)
print(RidgePlot(object = ucb, features.plot = top2$gene))
dev.off()


pdf('/mnt/data5/BGI/UCB/tangchao/Seurat/SJ/New_markers_VlnPlot_top2.pdf',18,12)
print(VlnPlot(object = ucb, features.plot = top2$gene))
dev.off()


customize_Seurat_FeaturePlot <- function(p, alpha.use = 1, gradient.use = c("yellow", "red"), expression.threshold = 0, is.log1p.transformed = F) {
  
  #### Main function ####
  main_function <- function(p = p, alpha.use = alpha.use, gradient.use = gradient.use, expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed) {
    
    # Order data by gene expresion level
    p$data <- p$data[order(p$data$gene),]
    
    # Define lower limit of gene expression level
    if (isTRUE(is.log1p.transformed)) {
      expression.threshold <- expression.threshold
    } else {
      expression.threshold <- log1p(expression.threshold)
    }
    
    # Compute maximum value in gene expression
    max.exp <- max(p$data$gene)
    
    # Fill points using the gene expression levels
    p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p <- p + scale_colour_gradientn(colours = gradient.use, guide = F, limits = c(expression.threshold, max.exp), na.value = "grey") +
      scale_fill_gradientn(colours = gradient.use, name = expression(atop(Expression, (log))), limits = c(expression.threshold, max.exp), na.value = "grey") +
      scale_alpha_continuous(range = alpha.use, guide = F)
  }
  
  #### Execution of main function ####
  # Apply main function on all features
  p <- lapply(X = p, alpha.use = alpha.use, gradient.use = gradient.use, 
              expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(p))))
}


p <- FeaturePlot(object = ucb, features.plot = top2$gene, cols.use = c("gray97", "blue"),pch.use = 16, no.legend = F, do.return = T)
print(customize_Seurat_FeaturePlot(p,expression.threshold = 0, gradient.use = c("gray97", "blue")))












































