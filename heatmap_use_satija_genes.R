
setwd("/mnt/data5/BGI/UCB/tangchao/Seurat/Satija")
refe <- read.csv("referencce/inline-supplementary-material-2.csv", sep = ",",)
TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", sep = "\t", header = T, row.names = 1)

library("AnnotationDbi")
library("org.Hs.eg.db")
refe$gene_id <- mapIds(org.Hs.eg.db,
                     keys=as.character(refe$gene),
                     column="ENSEMBL",
                     keytype="SYMBOL",
                     multiVals="first")


TPM_tu <- TPM[unique(as.character(na.omit(refe$gene_id))), ]
dim(TPM_tu)
# [1]  349 3574
identical(as.character(refe[!is.na(refe$gene_id) & !duplicated(refe$gene_id), ]$gene_id), row.names(TPM_tu))
# [1] TRUE
row.names(TPM_tu) <- as.character(refe[!is.na(refe$gene_id) & !duplicated(refe$gene_id), ]$gene)


library(pheatmap)

pdf("Satija1.pdf")
pheatmap(as.matrix(TPM_tu), show_rownames = F, show_colnames = F)
dev.off()

pdf("Satija2.pdf")
pheatmap(log(as.matrix(TPM_tu)+1), show_rownames = F, show_colnames = F)
dev.off()

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", sep = "\t", header = F, row.names = 1) -> celltype
colnames(celltype) <- "CellType"
celltype$sample <- substr(row.names(celltype),1,4)

celltype <- celltype[order(celltype$CellType), ]

pdf("Satija3.pdf")
pheatmap(log(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, , cluster_cols = F, annotation_col = celltype)
dev.off()



hts <- t(apply(TPM_tu,1,function(x)(x-mean(x))/sd(x)))
max(hts)
min(hts)
hts[hts>(5)] = 5
#
pdf("Satija4.pdf")
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = celltype)
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype)
dev.off()


mean(apply(TPM_tu,1,mean))
# [1] 487.8391

sum(apply(TPM_tu,1,mean) > median(apply(TPM_tu,1,mean)))

TPM_tu_high <- TPM_tu[apply(TPM_tu,1,mean) >= median(apply(TPM_tu,1,mean)), ]
TPM_tu_low <- TPM_tu[apply(TPM_tu,1,mean) <= median(apply(TPM_tu,1,mean)), ]

pdf("Satija5.pdf")
pheatmap(log(as.matrix(TPM_tu_high[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, , cluster_cols = F, annotation_col = celltype)
pheatmap(log(as.matrix(TPM_tu_low[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, , cluster_cols = F, annotation_col = celltype)
dev.off()



geneType <- na.omit(refe[,4:5])
geneType <- geneType[!duplicated(geneType$gene_id),]
row.names(geneType) <- as.character(geneType$gene_id)
geneType <- subset.data.frame(geneType, select = 1)

pdf("Satija6.pdf")
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F)
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F, cluster_rows = F)
dev.off()


pdf("Satija7.pdf")
pheatmap(log10(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F)
pheatmap(log10(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F, cluster_rows = F)
dev.off()




#### Add Unknown

setwd("/mnt/data5/BGI/UCB/tangchao/Seurat/Satija")
refe <- read.csv("referencce/inline-supplementary-material-2.csv", sep = ",",)
TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", sep = "\t", header = T, row.names = 1)

read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", sep = "\t", header = F, row.names = 1) -> celltype
colnames(celltype) <- "CellType"
celltype$sample <- substr(row.names(celltype),1,4)

celltype <- celltype[order(celltype$CellType), ]

colnames(TPM)[!colnames(TPM) %in% row.names(celltype)]
unknown <- data.frame(CellType = "Unknown", sample = substr(colnames(TPM)[!colnames(TPM) %in% row.names(celltype)],1,4),row.names = colnames(TPM)[!colnames(TPM) %in% row.names(celltype)])
celltype <- rbind(celltype, unknown)
celltype$sample <- as.factor(celltype$sample)

library("AnnotationDbi")
library("org.Hs.eg.db")
refe$gene_id <- mapIds(org.Hs.eg.db,
                     keys=as.character(refe$gene),
                     column="ENSEMBL",
                     keytype="SYMBOL",
                     multiVals="first")


TPM_tu <- TPM[unique(as.character(na.omit(refe$gene_id))), ]
dim(TPM_tu)
# [1]  349 3574
identical(as.character(refe[!is.na(refe$gene_id) & !duplicated(refe$gene_id), ]$gene_id), row.names(TPM_tu))
# [1] TRUE
row.names(TPM_tu) <- as.character(refe[!is.na(refe$gene_id) & !duplicated(refe$gene_id), ]$gene_id)

hts <- t(apply(TPM_tu,1,function(x)(x-mean(x))/sd(x)))
max(hts)
min(hts)
hts[hts>(5)] = 5


geneType <- na.omit(refe[,4:5])
geneType <- geneType[!duplicated(geneType$gene_id),]
row.names(geneType) <- as.character(geneType$gene_id)
geneType <- subset.data.frame(geneType, select = 1)

library(pheatmap)


pdf("Satija8.pdf")
pheatmap(log10(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F)
pheatmap(log10(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_rows = F)
pheatmap(log10(as.matrix(TPM_tu[,row.names(celltype)])+1), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F, cluster_rows = F)

pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F)
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_rows = F)
pheatmap(as.matrix(hts[,row.names(celltype)]), show_rownames = F, show_colnames = F, annotation_col = celltype, annotation_row = geneType, cluster_cols = F, cluster_rows = F)
dev.off()




