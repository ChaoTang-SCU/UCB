
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")

gtf <- read.table(pfgtf, head = F, sep = "\t", stringsAsFactors = F)
exon <- subset(gtf,gtf$V3=="exon")
exon$gene_id <- substr(exon$V9,9,23)
exon <- unique(exon)
colnames(exon) <- c("GeneID","Chr","Start","End","Strand")
write.table(exon, "/mnt/data5/BGI/UCB/tangchao/BCseq/hg38_Ensembl_exon.txt", sep = "\t", row.names = F, quote = F)


for f in `awk '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt` \
	do featureCounts -a /mnt/data5/BGI/UCB/tangchao/BCseq/hg38_Ensembl_exon.txt \
					 -o /mnt/data5/BGI/UCB/tangchao/BCseq/featurecount_ensembl/$f.epCount \
					 /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam \
					 -F SAF -f --read2pos 5 -O -T 10 \
done


#### R
library("data.table")
mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/BCseq/featurecount_ensembl")

library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="epCount$", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x, header=TRUE, select = 7, sep = "\t")
    names(tmp) <- strsplit(names(tmp), "/")[[1]][9]
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}
system.time(te<-multmerge(mypath))

save(te, file = "/mnt/data5/BGI/UCB/tangchao/BCseq/featureCounts_raw_merge_ensembl.RData")

filenames=list.files(path=mypath, pattern="epCount$", full.names=TRUE)
tmp<-fread(filenames[1], header=TRUE, select = 1:6, sep = "\t")

all_feature <- cbind(tmp, te)
setkey(all_feature, "Geneid")

gene_length <- all_feature[,.(gene_length = sum(Length)), by=Geneid]
gene_str <- all_feature[,.(gene_str = min(Start)), by=Geneid]
gene_chr <- all_feature[,.(gene_chr = unique(Chr)), by=Geneid]

all_mean <- all_feature[, lapply(.SD, mean), by=Geneid, .SDcols = colnames(all_feature)[grep("UCB", colnames(all_feature))]]
all_var <- all_feature[, lapply(.SD, var), by=Geneid, .SDcols = colnames(all_feature)[grep("UCB", colnames(all_feature))]]
names(all_var) <- c("Geneid",paste(colnames(all_var)[grep("UCB", colnames(all_var))], "_var", sep = ""))
names(all_mean) <- c("Geneid",paste(colnames(all_mean)[grep("UCB", colnames(all_mean))], "_mean", sep = ""))
all_mean_var <- merge(all_mean, all_var, by = "Geneid")
all_mean_var <- all_mean_var[,order(names(all_mean_var)), with = F]

m1 <- merge(gene_chr, gene_str, by = "Geneid")
m2 <- merge(m1, gene_length, by = "Geneid")
m2[, name := m2[[1]]]

m3 <- merge(m2, all_mean_var, by = "Geneid")
gene_cnt_summary <- as.data.frame(m3)

write.table(gene_cnt_summary, "/mnt/data5/BGI/UCB/tangchao/BCseq/gene_cnt_summary_ensembl", sep = "\t", row.names = F, quote = F)

system("/mnt/data4/software/BCseq/BCseq gene_cnt_summary_ensembl bcseq_result_ensembl 3039 58051 1 1-3039")

bcseq_result <- fread("/mnt/data5/BGI/UCB/tangchao/BCseq/bcseq_result_ensembl", sep = " ", header = T)
dim(bcseq_result)

bcseq_result[, gid := m3[[1]]]
test <- data.frame(as.data.frame(bcseq_result)[,-1], row.names = bcseq_result[[1]])
test <- test[test$mark>0,]

all_bc <- test[,grep("\\.1", colnames(test))]

colnames(all_bc) <- colnames(all_feature)[grep("UCB", colnames(all_feature))]

rs <- log10(rowSums(all_bc))
bc_tu <- all_bc[rs >= 2, ]

library(Seurat)
ucb <- CreateSeuratObject(raw.data = bc_tu, min.cells = 3, min.genes=500, is.expr = 1, project="UCB_cluster")
# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0)
length(x = ucb@var.genes)
# [1] 5953

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
ucb <- RunTSNE(object = ucb, dims.use = 1:22, do.fast = TRUE)
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
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()





# *** 0. Load the data ***
mat_all = all_bc

## Convert GeneID to GeneName in matrix
gname = read.table( "/mnt/data5/BGI/UCB/ExpMat_NewID/Gene_info.txt",header=T,row.names=1,sep = "\t")
gname = gname[!grepl("ERCC-",rownames(gname)), ]
#gname = gname[grepl("protein_coding|^TR_|^IG_",gname$Biotype), ]
gname = gname[grepl("protein_coding",gname$Biotype), ]
gname = gname[!duplicated(gname$GeneName), ]

flag = rownames(mat_all) %in% rownames(gname)
mat_all = mat_all[flag,]
rownames(mat_all) <- gname[rownames(mat_all),1]


## Split matrix by each batch
batch = read.table( "/mnt/data5/BGI/UCB/ExpMat_NewID/batch.rmSuspB.txt", header=F, stringsAsFactors=F, sep = "\t")
batch <- batch[batch$V1 %in% colnames(mat_all), ]
b1 = batch[batch$V2==1, "V1"]
b2 = batch[batch$V2==2, "V1"] 
b3 = batch[batch$V2==3, "V1"] 
b4 = batch[batch$V2==4, "V1"] 
mat_b1 = mat_all[, b1]
mat_b2 = mat_all[, b2]
mat_b3 = mat_all[, b3]
mat_b4 = mat_all[, b4]

# *** 1. Setup the Seurat objects ***
s1 = CreateSeuratObject( raw.data=mat_b1, min.cells = 3, min.genes=500, is.expr=1, project="Batch1")
s2 = CreateSeuratObject( raw.data=mat_b2, min.cells = 3, min.genes=500, is.expr=1, project="Batch2")
s3 = CreateSeuratObject( raw.data=mat_b3, min.cells = 3, min.genes=500, is.expr=1, project="Batch3")
s4 = CreateSeuratObject( raw.data=mat_b4, min.cells = 3, min.genes=500, is.expr=1, project="Batch4")
s1@meta.data$batch = "UCB1"
s2@meta.data$batch = "UCB3"
s3@meta.data$batch = "UCB3S"
s4@meta.data$batch = "UCB4"

## Normalization and Scaling 
s1 = NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 1000000)
s2 = NormalizeData(s2, normalization.method = "LogNormalize", scale.factor = 1000000)
s3 = NormalizeData(s3, normalization.method = "LogNormalize", scale.factor = 1000000)
s4 = NormalizeData(s4, normalization.method = "LogNormalize", scale.factor = 1000000)
s1 = ScaleData(s1)
s2 = ScaleData(s2)
s3 = ScaleData(s3)
s4 = ScaleData(s4)

## Determine genes to use for CCA, must be highly variable in at least 2 datasets
s1 = FindVariableGenes(s1, do.plot = F)
s2 = FindVariableGenes(s2, do.plot = F)
s3 = FindVariableGenes(s3, do.plot = F)
s4 = FindVariableGenes(s4, do.plot = F)
ob.list = list( s1,s2,s3,s4 )
genes.use = c()
for (i in 1:length(ob.list)) {
  genes.use = c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 2000))
}
genes.use = names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use = genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}
length(genes.use) #1767 genes


# *** 2. Perform a canonical correlation analysis (CCA) ***
Combine = RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 30)
pdf("Raw_CC.pdf",18,8)
p1 = DimPlot(object = Combine, reduction.use = "cca", group.by = "batch", pt.size = 0.5, do.return = TRUE)
p2 = VlnPlot(object = Combine, features.plot = "CC1", group.by = "batch", do.return = TRUE)
plot_grid(p1, p2)
dev.off()

## CC Selection
pdf("CC_selection1.pdf",10,8)
MetageneBicorPlot(Combine, grouping.var = "batch", dims.eval = 1:30)  #1-20
dev.off()
pdf("CC_selection2.pdf",10,20)
DimHeatmap(object = Combine, reduction.type = "cca", cells.use = 500, dim.use = 1:30, do.balanced = TRUE)  #1-20
dev.off()

## Run rare non-overlapping filtering
Combine = CalcVarExpRatio(object = Combine, reduction.type = "pca", grouping.var = "batch", dims.use = 1:20)
Combine = SubsetData(object = Combine, subset.name = "var.ratio.pca", accept.low = 0.5)


# *** 3. Align the CCA subspaces ***
Combine = AlignSubspace(object = Combine, reduction.type = "cca", grouping.var = "batch", dims.align = 1:20)
pdf("Align_CC.pdf",10,8)
p1 = VlnPlot(object = Combine, features.plot = "ACC1", group.by = "batch", do.return = TRUE)
p2 = VlnPlot(object = Combine, features.plot = "ACC2", group.by = "batch", do.return = TRUE)
plot_grid(p1, p2)
dev.off()


# *** 4. t-SNE and Clustering ***
Combine = FindClusters(Combine, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = T, resolution = 0.8)
Combine = RunTSNE(Combine, reduction.use = "cca.aligned", dims.use = 1:20)
pdf("tSNE.pdf",10,8)
TSNEPlot(Combine, do.return = T, pt.size = 1.5, group.by = "batch")
TSNEPlot(Combine, do.label = T, do.return = T, pt.size = 1.5)
dev.off()

save(Combine, file = "save.Rda")



tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)
dim(tsne)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)
Cell_type <- Cell_type[row.names(tsne),]
dim(Cell_type)

tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()















library("sva")

#RNAseq - parametric ComBat

all_bc

exprs <- data.frame(log2(1+all_bc))


#Define center as batch
batch = read.table( "/mnt/data5/BGI/UCB/ExpMat_NewID/batch.rmSuspB.txt", header=F, stringsAsFactors=F, sep = "\t", row.names = 1)
sample_batch <- batch[colnames(exprs), ]

combat <- ComBat(dat=exprs, batch=sample_batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

rs <- rowSums(combat)

exprs_tu <- (2^combat)-1
rs <- rowSums(exprs_tu)

library(Seurat)
ucb <- CreateSeuratObject(raw.data = exprs_tu, min.cells = 3, min.genes=500, is.expr = 1, project="UCB_cluster")
# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.24, x.high.cutoff = 8, y.cutoff = 0)
length(x = ucb@var.genes)
# [1] 2954

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
#TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
#TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

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
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()
















batch = read.table( "/mnt/data5/BGI/UCB/ExpMat_NewID/batch.rmSuspB.txt", header=F, stringsAsFactors=F, sep = "\t", row.names = 1)
sample_batch <- batch[colnames(bc_tu), ]

combat <- ComBat(dat=bc_tu, batch=sample_batch, mod=NULL, par.prior=FALSE, prior.plots=FALSE, mean.only=T)

rs <- rowSums(combat)

library(Seurat)
ucb <- CreateSeuratObject(raw.data = combat, min.cells = 3, min.genes=500, is.expr = 1, project="UCB_cluster")
# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.24, x.high.cutoff = 8, y.cutoff = 0)
length(x = ucb@var.genes)
# [1] 2439
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0)
length(x = ucb@var.genes)
# [1] 5750

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
#TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
#TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

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
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()



































