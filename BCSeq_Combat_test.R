library(data.table)
load(file = "/mnt/data5/BGI/UCB/tangchao/BCseq/featureCounts_raw_merge_ensembl.RData")

mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/BCseq/featurecount_ensembl")
filenames=list.files(path=mypath, pattern="epCount$", full.names=TRUE)
tmp<-fread(filenames[1], header=TRUE, select = 1:6, sep = "\t")

all_feature <- cbind(tmp, te)
setkey(all_feature, "Geneid")

bcseq_result <- fread("/mnt/data5/BGI/UCB/tangchao/BCseq/bcseq_result_ensembl", sep = " ", header = T)
dim(bcseq_result)

bcseq_result[, gid := unique(all_feature[[1]])]
test <- data.frame(as.data.frame(bcseq_result)[,-1], row.names = bcseq_result[[1]])
test <- test[test$mark>0,]

all_bc <- test[,grep("\\.1", colnames(test))]

colnames(all_bc) <- colnames(all_feature)[grep("UCB", colnames(all_feature))]

save(all_bc, file = "/mnt/data5/BGI/UCB/tangchao/BCseq/all_bc.RData")

rs <- log10(rowSums(all_bc))
bc_tu <- all_bc[rs >= 2, ]


library("sva")

#RNAseq - parametric ComBat

#exprs <- data.frame(log2(1+all_bc))


#Define center as batch
batch = read.table( "/mnt/data5/BGI/UCB/ExpMat_NewID/batch.rmSuspB.txt", header=F, stringsAsFactors=F, sep = "\t", row.names = 1)
sample_batch <- batch[colnames(all_bc), ]

combat <- ComBat(dat=bc_tu, batch=sample_batch, mod=NULL, par.prior=FALSE, prior.plots=FALSE)

save(combat, file = "/mnt/data5/BGI/UCB/tangchao/BCseq/High_BCSeq_combat_nonpar.RData")


rs <- rowSums(combat)

pdf("~/rs.hist.pdf")
hist(log10(rs),breaks=100)
dev.off()


bc_tu <- combat[log10(rs) >= 2, ]

library(Seurat)
ucb <- CreateSeuratObject(raw.data = bc_tu, min.cells = 3, min.genes=500, is.expr = 1, project="UCB_cluster")
# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0)
length(x = ucb@var.genes)
# [1] 5071

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
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

dim(Cell_type)
tsne <- tsne[row.names(Cell_type), ]
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/BCseq/Seurat/enseembl_tsne4.pdf")
pdf("~/enseembl_tsne4.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()



















