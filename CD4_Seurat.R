
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

cd4_psi <- psi_sj_same_start_table[,row.names(Cell_type[Cell_type$CellType == "CD4+ T Cells",])]
dim(cd4_psi)
#[1] 460327   1678
rm("psi_sj_same_start_table")
rm("psi_sj_same_end_table")
gc()
cd4_psi <- cd4_psi[rowSums(!is.na(cd4_psi))>0,]
dim(cd4_psi)
#[1] 320874   1678
cd4_psi <- cd4_psi[rowSums(cd4_psi, na.rm = T)>0,]
dim(cd4_psi)
#[1] 286661   1678
sd4 <- apply(cd4_psi,1,function(x) sd(x, na.rm = T))
sd4[is.na(sd4)] <- 0
cd4_psi <- cd4_psi[sd4>0,]
dim(cd4_psi)
#[1] 87141   1678



library(Seurat)
library(dplyr)
library(Matrix)

cd4_psi[is.na(cd4_psi)] <- 0
save(cd4_psi, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/cd4_psi.RData")
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = cd4_psi, min.cells = 10, min.genes=10, is.expr = 0, project="CD4_cluster")

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


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_1.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


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
#[1] "7:142792081-142801943" "6:31269526-31355106"   "6:31270086-31355316"
#[4] "13:95371214-95619324"  "2:97118518-97142639"   "2:97209710-97211545"
#[7] "2:95899337-95904760"

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_cluster0.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==0])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_cluster1.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==1])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_cluster2.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==2])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_cluster3.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==3])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_novel.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T, pt.size = 3)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/CD4_tSNE_New_markers_unannotated.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$annotation==0,"sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T, pt.size = 3)
dev.off()



cd4_psi["13:95371214-95619324",1:5]
cd4_psi["7:142792081-142801943",1:5]

cd4_psi["13:95371214-95619324",cd4_psi["13:95371214-95619324",]>0]
table(substr(names(cd4_psi["13:95371214-95619324",cd4_psi["13:95371214-95619324",]>0]),1,4))
#UCB1 UCB3 UCB4
#   7  141   48
table(substr(names(cd4_psi["7:142792081-142801943",cd4_psi["7:142792081-142801943",]>0]),1,4))
#UCB1 UCB3 UCB4
#  15   28  460

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te_no_NA_and_normalized.RData")
te <- te[c("13:95371214-95619324","7:142792081-142801943"),colnames(cd4_psi)]
colnames(te[,order(te[1,], decreasing=T)])[1:10]
# [1] "UCB3.02123" "UCB3.02023" "UCB3.01942" "UCB3.01714" "UCB3.02205"
# [6] "UCB3.01980" "UCB3.02114" "UCB3.01551" "UCB3.01791" "UCB3.02196"
colnames(te[,order(te[2,], decreasing=T)])[1:10]
# [1] "UCB4.00843" "UCB4.01290" "UCB4.00499" "UCB4.00868" "UCB4.01228"
# [6] "UCB4.01213" "UCB4.00159" "UCB4.01109" "UCB4.00116" "UCB4.00647"

cd4_psi["13:95371214-95619324", colnames(te[,order(te[1,], decreasing=T)])[1:10]]
#UCB3.02123 UCB3.02023 UCB3.01942 UCB3.01714 UCB3.02205 UCB3.01980 UCB3.02114
# 0.0000000  0.2571429  0.4264706  0.2888889  0.5000000  0.4177215  0.5925926
#UCB3.01551 UCB3.01791 UCB3.02196
# 0.2500000  0.4583333  0.5833333

cd4_psi["7:142792081-142801943", colnames(te[,order(te[2,], decreasing=T)])[1:10]]
#UCB4.00843 UCB4.01290 UCB4.00499 UCB4.00868 UCB4.01228 UCB4.01213 UCB4.00159
#0.46860465 0.27365591 0.29719189 0.37096774 0.06995004 0.42465753 0.28744186
#UCB4.01109 UCB4.00116 UCB4.00647
#0.36483516 0.31319555 0.12773109


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
    	outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 1)

sashimi(junction = "7:142792081-142801943", cell = colnames(te[,order(te[2,], decreasing=T)])[1:10], 
    	outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/CD4_Seurat/", ann_height = 10, left_expand = 100, right_expand = 100, min = 10)












