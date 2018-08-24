load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

nsj <- sj[sj$annotation == 0, "sj"]
psi_tu <- psi_sj_same_end_table[nsj[nsj %in% row.names(psi_sj_same_end_table)], row.names(Cell_type)]
psi_tu <- psi_tu[rowSums(!is.na(psi_tu))>=10, ]


library(Seurat)
library(dplyr)
library(Matrix)

psi_tu[is.na(psi_tu)] <- 0

# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_tu, min.cells = 1, min.genes=10, is.expr = 0, project="CD8_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = .23, x.high.cutoff = 6, y.cutoff = -5)
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
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:20, resolution = 0.7, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)


tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$Individual <- as.factor(Cell_type[row.names(tsne),]$Individual)
tsne$CellType <- as.factor(Cell_type[row.names(tsne),]$CellType)

library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/tSNE_1.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/Seurat_result.RData")


# . Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers <- ucb.markers[ucb.markers$p_val_adj<0.05,]
ucb.markers %>% group_by(cluster) %>% top_n(10) -> top10
table(ucb.markers$cluster)
# 0  1  2  3  4  5  6  7  8  9 10 11 12
# 4  1  7  6 10  5  4  8  3  3  4  3 20

load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
markers_info <- merge(x = ucb.markers, y = sj[sj$sj %in% ucb.markers$gene, c(1,12,42:45)], by.x = "gene", by.y = "sj", all.x = T)

save(markers_info, file = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/markers_info.RData")

#write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.txt",sep = "\t",quote = F)
#pdf("New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
#DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)

ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene
# [1] "13:95371214-95619324"
ucb.markers[ucb.markers$gene %in% sj[sj$annotation==0,"sj"],]$gene
#[1] "7:142792081-142801943" "2:127076707-127081831" "13:95371214-95619324"

pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/Novel_markers_feature.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = T)
FeaturePlot(object = ucb, features.plot = ucb.markers[ucb.markers$gene %in% sj[sj$novel=="Y","sj"],]$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()


markers_info[markers_info$novel == "Y", ]

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")

tsne_ge <- as.data.frame(Combine@dr$tsne@cell.embeddings)
dim(tsne_ge)


TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", header = T, row.names = 1, stringsAsFactors = F)
TPM["ENSG00000116824",]

tmp_data <- data.frame(tsne_ge, GE = log2(as.numeric(TPM["ENSG00000116962",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000116962_1:236055179-236295737.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("1:236055179-236295737")+
 		theme(panel.background = element_blank())
dev.off()



tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM["ENSG00000090382",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000090382_12:69116358-69353570.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000090382")+
 		theme(panel.background = element_blank())
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000090382_Feature.pdf", height = 6,width = 6)
FeaturePlot(object = Combine, features.plot = "LYZ", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()



tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM["ENSG00000257764",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000257764_12:69116358-69353570.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000257764_AC020656.1")+
 		theme(panel.background = element_blank())
dev.off()



tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM["ENSG00000121766",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000121766_1:31132361-31355110.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000121766_ZCCHC17")+
 		theme(panel.background = element_blank())
FeaturePlot(object = Combine, features.plot = "ZCCHC17", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()



tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM[" ENSG00000205913",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000205913_16:2747589-2893732.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle(" ENSG00000205913_SRRM2-AS1")+
 		theme(panel.background = element_blank())
FeaturePlot(object = Combine, features.plot = "SRRM2-AS1", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()


tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM[" ENSG00000162076",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000162076_16:2747589-2893732.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle(" ENSG00000162076_FLYWCH2")+
 		theme(panel.background = element_blank())
FeaturePlot(object = Combine, features.plot = "FLYWCH2", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()

#ENSG00000230795

tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM["ENSG00000230795",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000230795.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000230795_HLA-K")+
 		theme(panel.background = element_blank())
#FeaturePlot(object = Combine, features.plot = "HLA-K", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()


tmp_data <- data.frame(tsne_ge, GE = log10(as.numeric(TPM["ENSG00000230795",row.names(tsne_ge)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/ENSG00000230795.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000230795_HLA-K")+
 		theme(panel.background = element_blank())
#FeaturePlot(object = Combine, features.plot = "HLA-K", cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()



load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(no_NA).RData")
# 1:31132361−31355110
# 6:29927351−29934184
# 22:22906554−22919750
# 12:69146556−69353736
# 1:26479778−26520853
# 8:7916460−8097116
# 6:37186872−37286458

te_table <- te_table[c("1:26479778−26520853","8:7916460−8097116","6:37186872−37286458"),row.names(Cell_type)]

plot_tab <- t(apply(te_table,1,function(x) head(names(x)[order(x, decreasing=T)])))


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


for(i in 1:nrow(plot_tab)){
  print(i)
  sashimi(junction = row.names(plot_tab)[i], cell = plot_tab[i,], 
    	  outDir = "/mnt/data5/BGI/UCB/tangchao/PSI_Seurat/All_unannotated_sj/", ann_height = 20, left_expand = 100, right_expand = 100, min = 10)
}




