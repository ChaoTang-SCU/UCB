#### feature counts

for f in `awk '{print $1}' /mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt` \
	do featureCounts -a /mnt/data4/software/subread-1.6.2-Linux-x86_64/annotation/hg38_RefSeq_exon.txt \
					 -o /mnt/data5/BGI/UCB/tangchao/BCseq/featurecount/$f.epCount \
					 /mnt/data5/BGI/UCB/tangchao/data/STAR/$f/$f\_Aligned.sortedByCoord.out.bam \
					 -F SAF -f --read2pos 5 -O -T 10 \
done


#### R
library("data.table")
mypath<-file.path("/mnt/data5/BGI/UCB/tangchao/BCseq/featurecount")

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

save(te, file = "/mnt/data5/BGI/UCB/tangchao/BCseq/featureCounts_raw_merge.RData")

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



library(clusterProfiler)

keytypes(org.Hs.eg.db)
#gene_enid <- mapIds(org.Hs.eg.db, keys = as.character(all_mean_var[[1]]), column = "ENSEMBL", keytype = "ENTREZID")
gene_id <- bitr(all_mean_var[[1]], fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
colnames(gene_id) <- c("Geneid", "gene_name")


gene_chr <- as.data.frame(gene_chr)
gene_chr <- gene_chr[!duplicated(gene_chr$Geneid),]
gene_str <- as.data.frame(gene_str)
gene_length <- as.data.frame(gene_length)

m1 <- merge(gene_chr, gene_str, by = "Geneid")
m2 <- merge(m1, gene_length, by = "Geneid")
m3 <- merge(m2, gene_id, by = "Geneid")

gene_cnt_summary <- merge(m3, as.data.frame(all_mean_var), by = "Geneid")

write.table(gene_cnt_summary, "/mnt/data5/BGI/UCB/tangchao/BCseq/gene_cnt_summary", sep = "\t", row.names = F, quote = F)



system("/mnt/data4/software/BCseq/BCseq gene_cnt_summary bcseq_result 3039 27805 1 1-3039")


bcseq_result <- read.table("/mnt/data5/BGI/UCB/tangchao/BCseq/bcseq_result", sep = " ", header = T, row.names = 1)
dim(bcseq_result)
#[1] 27805  9123

test <- bcseq_result[bcseq_result$mark >0, ]

all_bc <- test[,grep("\\.1", colnames(test))]
colnames(all_bc) <- colnames(all_feature)[grep("UCB", colnames(all_feature))]


library(Seurat)

ucb <- CreateSeuratObject(raw.data = all_bc, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = .8)
length(x = ucb@var.genes)
# [1] 3017

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
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:7, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:7, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne <- tsne[row.names(Cell_type), ]
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

#pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_Cell_specific_SJ.pdf")
pdf("/mnt/data5/BGI/UCB/tangchao/BCseq/Seurat/tsne3.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()






















