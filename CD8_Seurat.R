
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

cd8_psi <- psi_sj_same_start_table[,row.names(Cell_type[Cell_type$CellType == "CD4+ T Cells",])]
dim(cd8_psi)
#[1] 460327   1678
rm("psi_sj_same_start_table")
rm("psi_sj_same_end_table")
gc()
cd8_psi <- cd8_psi[rowSums(!is.na(cd8_psi))>0,]
dim(cd8_psi)
#[1] 238282    579
cd8_psi <- cd8_psi[rowSums(cd8_psi, na.rm = T)>0,]
dim(cd8_psi)
#[1] 187369    579
sd8 <- apply(cd8_psi,1,function(x) sd(x, na.rm = T))
sd8[is.na(sd4)] <- 0
cd8_psi <- cd8_psi[sd8>0,]
dim(cd8_psi)
#[1] 115292    579



library(Seurat)
library(dplyr)
library(Matrix)

cd8_psi[is.na(cd8_psi)] <- 0
save(cd8_psi, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/cd8_psi.RData")
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = cd8_psi, min.cells = 10, min.genes=10, is.expr = 0, project="CD8_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 2.3, y.cutoff = -5)
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
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.3, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_tSNE_1.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_Seurat_result.RData")

# . Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers <- ucb.markers[ucb.markers$p_val_adj<0.05,]
ucb.markers %>% group_by(cluster) %>% top_n(10) -> top10
table(ucb.markers$cluster)
#  0   1   2   3
# 61  13 406  34
#write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.txt",sep = "\t",quote = F)
#pdf("New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
#DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene
# [1] "13:95371214-95619324"
ucb.markers[ucb.markers$gene %in% sj[sj$annotation==0,"sj"],]$gene
#[1] "7:142792081-142801943" "2:127076707-127081831" "13:95371214-95619324"

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_tSNE_New_markers_cluster0.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==0])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_tSNE_New_markers_cluster1.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==1])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_tSNE_New_markers_novel.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T, pt.size = 3)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_tSNE_New_markers_unannotated.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$annotation==0,"sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T, pt.size = 3)
dev.off()



cd8_psi["13:95371214-95619324",1:5]
cd8_psi["7:142792081-142801943",1:5]
cd8_psi["2:127076707-127081831",1:5]

cd8_psi["13:95371214-95619324",cd8_psi["13:95371214-95619324",]>0]
table(substr(names(cd8_psi["13:95371214-95619324",cd8_psi["13:95371214-95619324",]>0]),1,4))
#UCB1 UCB3 UCB4
#   7   72   12
table(substr(names(cd8_psi["7:142792081-142801943",cd8_psi["7:142792081-142801943",]>0]),1,4))
#UCB1 UCB3 UCB4
#   4   11  145
table(substr(names(cd8_psi["2:127076707-127081831",cd8_psi["2:127076707-127081831",]>0]),1,4))
#UCB1 UCB3 UCB4
#  18   68   48

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")
te <- te[c("13:95371214-95619324","7:142792081-142801943","2:127076707-127081831"),colnames(cd8_psi)]
colnames(te[,order(te[1,], decreasing=T)])[1:10]
# [1] "UCB3.02009" "UCB3.02052" "UCB3.01447" "UCB3.02162" "UCB4.00489"
# [6] "UCB3.02120" "UCB3.01894" "UCB3.01720" "UCB3.01839" "UCB3.01869"
colnames(te[,order(te[2,], decreasing=T)])[1:10]
# [1] "UCB4.00844" "UCB4.00970" "UCB4.00817" "UCB4.01024" "UCB4.00098"
# [6] "UCB4.00650" "UCB4.00297" "UCB4.00172" "UCB4.00279" "UCB4.01210"
colnames(te[,order(te[3,], decreasing=T)])[1:10]
# [1] "UCB3.01882" "UCB1.00089" "UCB3.01562" "UCB1.00092" "UCB4.00967"
# [6] "UCB4.00172" "UCB3.02009" "UCB1.00087" "UCB4.00371" "UCB4.00725"

cd8_psi["13:95371214-95619324", colnames(te[,order(te[1,], decreasing=T)])[1:10]]
#UCB3.02009 UCB3.02052 UCB3.01447 UCB3.02162 UCB4.00489 UCB3.02120 UCB3.01894
# 0.3943662  0.5238095  0.2676056  0.2040816  1.0000000  0.5000000  0.2911392
#UCB3.01720 UCB3.01839 UCB3.01869
# 0.0000000  0.5000000  0.0000000


cd8_psi["7:142792081-142801943", colnames(te[,order(te[2,], decreasing=T)])[1:10]]
#UCB4.00844 UCB4.00970 UCB4.00817 UCB4.01024 UCB4.00098 UCB4.00650 UCB4.00297
# 0.3490338  0.3812500  0.4352941  0.2961672  0.3130435  0.2267846  0.2838710
#UCB4.00172 UCB4.00279 UCB4.01210
# 0.2319224  0.2748299  0.0747992


cd8_psi["2:127076707-127081831", colnames(te[,order(te[3,], decreasing=T)])[1:10]]
#UCB3.01882 UCB1.00089 UCB3.01562 UCB1.00092 UCB4.00967 UCB4.00172 UCB3.02009
# 0.9851995  0.9704142  0.9654179  0.9587302  0.9705882  1.0000000  1.0000000
#UCB1.00087 UCB4.00371 UCB4.00725
# 1.0000000  0.9598930  1.0000000



sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, tsv = "/mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/", outDir = "~/", min = 10, ann_height = 4){
  if(dir.exists(outDir) == FALSE) dir.create(outDir)
  dir.create("/home/tangchao/test2")
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
    tmpdir <- "/home/tangchao/test2"
    work1 <- paste("system('", "cat ", paste(tsv,cell,".tsv",sep = "", collapse = " "), " > ", tmpdir,"/tmp.tsv", "')", sep = "")
    work2 <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi2.sh ", window, paste(tmpdir,"/tmp.tsv", sep = ""), Dir_out, min, ann_height, "')")
    work3 <- paste("system('", "rm",paste(tmpdir,"/tmp.tsv", "')", sep = ""))
    eval(parse(text = work1))
    eval(parse(text = work2))
    eval(parse(text = work3))
  }
}

sashimi(junction = "13:95371214-95619324", cell = colnames(te[,order(te[1,], decreasing=T)])[1:10], 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 1)

sashimi(junction = "7:142792081-142801943", cell = colnames(te[,order(te[2,], decreasing=T)])[1:10], 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 10)

sashimi(junction = "2:127076707-127081831", cell = colnames(te[,order(te[3,], decreasing=T)])[1:10], 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 10)



sashimi(junction = "13:95371214-95619324", cell = c('UCB3.02009','UCB3.02052','UCB3.01447','UCB3.02162','UCB4.00489','UCB3.02120','UCB3.01894','UCB3.01720','UCB3.01839','UCB3.01869'), 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 1)

sashimi(junction = "7:142792081-142801943", cell = c("UCB4.00844", "UCB4.00970", "UCB4.00817", "UCB4.01024", "UCB4.00098", "UCB4.00650", "UCB4.00297", "UCB4.00172", "UCB4.00279", "UCB4.01210"), 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 1)

sashimi(junction = "2:127076707-127081831", cell = c("UCB3.01882", "UCB1.00089", "UCB3.01562", "UCB1.00092", "UCB4.00967", "UCB4.00172", "UCB3.02009", "UCB1.00087", "UCB4.00371", "UCB4.00725"), 
        outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 1)


#### Only 2:127076707-127081831 has host gene:
# ENST00000409400

TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", header = T, row.names = 1, stringsAsFactors = F)
TPM["ENSG00000136717",]

tmp_data <- data.frame(tsne, GE = log10(as.numeric(TPM["ENSG00000136717",row.names(tsne)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD8_Seurat/CD8_marker_annotated_hostgene_2_127076707_127081831.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
    geom_point(aes(colour = GE), size = 2)+
    scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
    ggtitle("ENSG00000136717")+
    theme(panel.background = element_blank())
dev.off()








