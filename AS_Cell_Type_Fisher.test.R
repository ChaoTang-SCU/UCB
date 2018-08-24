load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
PSI <- SE_psi[, 13:ncol(SE_psi)]
row.names(PSI) <- SE_psi$loci
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)
PSI <- PSI[, row.names(Cell_type)]
dim(PSI)
# [1] 17670  3039
rs <- rowSums(!is.na(PSI))
sum(rs == 0)
#[1] 53
PSI_tu <- PSI[rs>=1, ]
PSI_Cells <- t(apply(PSI_tu, 1, function(x) table(na.omit(cbind(Cell_type, x))$CellType)))
PSI_Cellp <- t(apply(PSI_Cells, 1, function(x) x/table(Cell_type$CellType)))

fisher.test(matrix(c(31,109,43,2,4,10,607,1678,579,24,52,99),6,2))$p.value

fisher.test(matrix(c(as.numeric(PSI_Cells[9,]), as.numeric(table(Cell_type$CellType))-as.numeric(PSI_Cells[9,])),6,2))

#psi_p <- apply(PSI_Cells, 1, function(x) {fisher.test(matrix(c(as.numeric(x), as.numeric(table(Cell_type$CellType))-as.numeric(x)),6,2), simulate.p.value = TRUE, B = 1e5)$p.value})

psi_p <- vector()
for(i in 1:nrow(PSI_Cells)){
	print(paste(i, "of", nrow(PSI_Cells)))
	tab <- matrix(c(as.numeric(PSI_Cells[i,]), as.numeric(table(Cell_type$CellType))-as.numeric(PSI_Cells[i,])),6,2)
	if(sum(tab[,2])==0){
		psi_p[i] <- NA
	}else{
		psi_p[i] <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
	}	
}
save(psi_p, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/psi_p.RData")
hist(psi_p)
psi_p[is.na(psi_p)] <- 1

psi_p2 <- vector()
for(i in 1:nrow(PSI_Cells)){
	print(paste(i, "of", nrow(PSI_Cells)))
	tab <- matrix(c(as.numeric(PSI_Cells[i,1:3]), as.numeric(table(Cell_type$CellType))[1:3]-as.numeric(PSI_Cells[i,1:3])),3,2)
	if(sum(tab[,2])==0 | sum(tab[,1])==0){
		psi_p2[i] <- NA
	}else{
		psi_p2[i] <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
	}	
}
save(psi_p2, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/psi_p2.RData")

PSI_Cellp <- as.data.frame(PSI_Cellp)
PSI_Cellp$psi_p2 <- psi_p2

row.names(PSI_Cellp[head(order(PSI_Cellp$psi_p2),10),])
dim(PSI_tu[row.names(PSI_Cellp[head(order(PSI_Cellp$psi_p2),10),]),])
#[1]   10 3039

PSI_sub <- PSI_tu[row.names(PSI_Cellp[head(order(PSI_Cellp$psi_p2),100),]),]
tsne_tu <- merge(tsne, t(PSI_sub), by = 0, all.x=T)
tsne_tu <- data.frame(tsne_tu[,-1], row.names = tsne_tu[,1])

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Fisher_test_top100_featureplot.pdf", width=9, height=6)
for(i in 3:ncol(tsne_tu)){
  tsne_tu[rev(order(tsne_tu[,i], decreasing=T)),] -> CD45_R_r_tab1
  print(ggplot(data = CD45_R_r_tab1, aes(x = tSNE_1, y = tSNE_2))+
  geom_point(aes(colour = as.numeric(CD45_R_r_tab1[,i])), size = 1)+
  #scale_colour_gradient(high = 'red',low = 'Grey')+
  scale_colour_gradient(low = "red", high = "red",na.value = "grey50")+
  ggtitle(colnames(CD45_R_r_tab1)[i])+
  theme(legend.title = element_blank(),
      	axis.title = element_blank(),
        axis.text= element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()))
}
dev.off()



PSI_Cellp$CSI <- apply(PSI_Cellp[,1:3], 1, function(x){(max(x) - min(x))/min(x)})
PSI_Cellp$delta <- apply(PSI_Cellp[,1:3], 1, function(x){max(x) - min(x)})

plot(na.omit(PSI_Cellp)$psi_p2, na.omit(PSI_Cellp)$CSI, ylim = c(0,2), pch = 19, cex=.2)
plot(na.omit(PSI_Cellp)$psi_p2, na.omit(PSI_Cellp)$delta, ylim = c(0,2), pch = 19, cex=.2)

hist(PSI_Cellp[which(is.infinite(PSI_Cellp$CSI)), "psi_p2"])

test <- PSI_Cellp[PSI_Cellp$psi_p2 < 0.01, ]
test <- na.omit(test)

plot(test$psi_p2, test$CSI, ylim = c(0,2), pch = 19, cex=.2)
plot(test$psi_p2, test$delta, ylim = c(0,.05), pch = 19, cex=.2)

test <- test[test$CSI>1, ]
plot(test$CSI, test$delta, pch = 19, cex=.2)



PSI_sub <- PSI[row.names(PSI_Cellp[PSI_Cellp$psi_p2 > 0.01,]),]
sum(rowSums(!is.na(PSI_sub))==0)
#[1] 134
PSI_tu <- PSI_sub[rowSums(!is.na(PSI_sub))!=0, ]

PSI_tu[is.na(PSI_tu)] <- -1

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 1, min.genes=1, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = -0.1, x.high.cutoff = .8, y.cutoff = 0)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

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
tsne <- tsne[row.names(Cell_type),]
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Remove_Fisher_test_Significant_SE_Highly_variable_SE_Seurat_tSNE_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()


















PSI_sub <- PSI[row.names(PSI_Cellp[PSI_Cellp$psi_p2 < 0.05,]),]
sum(rowSums(!is.na(PSI_sub))==0)
#[1] 134
PSI_tu <- PSI_sub[rowSums(!is.na(PSI_sub))!=0, ]

PSI_tu[is.na(PSI_tu)] <- 0

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 1, min.genes=1, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = -0.1, x.high.cutoff = 4, y.cutoff = -10)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

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
tsne <- tsne[row.names(Cell_type),]
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Cell_Sep/Fisher_test_Significant_SE_Seurat_tSNE_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()








