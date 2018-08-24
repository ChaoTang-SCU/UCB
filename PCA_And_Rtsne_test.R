#### DSU

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


sashimi(junction = "15:44715702-44717606", cell = c("UCB4.01387", "UCB4.01433"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Sashimi/")






load("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")


A3SS_psi[A3SS_psi$as2 == "1:11029782-11029901_1:11029782-11030178", c("UCB4.01449", "UCB4.01416", "UCB4.01463", "UCB4.01486", "UCB1.00084")]
#         UCB4.01449 UCB4.01416 UCB4.01463 UCB4.01486 UCB1.00084
# name216  0.9480519  0.9861023          1        NaN          0

sashimi(junction = "1:11029782-11030178", cell = c("UCB4.01449", "UCB4.01416", "UCB4.01463", "UCB4.01486", "UCB1.00084"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Sashimi/")



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_PSI_new.RData")


A5SS_psi[A5SS_psi$as2 == "1:11087903-11088122_1:11087911-11088122", c("UCB4.00727", "UCB4.00239", "UCB1.00006")]
#         UCB4.00727 UCB4.00239 UCB1.00006
# name232   0.772093          1        NaN
sashimi(junction = "1:11087903-11088122", cell = c("UCB4.00727", "UCB4.00239"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/Sashimi/")





load("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")

AFE_psi[23,which(!is.na(AFE_psi[23,]))]

AFE_psi[AFE_psi$as2 == "1:67426286-67429435_1:67426286-67429987", c("UCB4.01496", "UCB4.00842")]
#         UCB4.01496 UCB4.00842
#name1390  0.0172144          0
sashimi(junction = "1:67426286-67429987", cell = c("UCB4.01496", "UCB4.00842"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Sashimi/")



load("/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi_new.RData")
ALE_psi[23,which(!is.na(ALE_psi[23,]))]

ALE_psi[ALE_psi$as2 == "15:43793084-43793686_15:43793084-43793709", c("UCB4.01420", "UCB1.00040", "UCB1.00083", "UCB4.01490")]
#            UCB4.01420 UCB1.00040 UCB1.00083 UCB4.01490
#		     0.2590847  0.2229796          0        NaN

sashimi(junction = "15:43793084-43793709", cell = c("UCB4.01420", "UCB1.00040", "UCB1.00083", "UCB4.01490"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/Sashimi/")



load("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")
MXE_psi[23,which(!is.na(MXE_psi[23,]))]

MXE_psi[MXE_psi$as2 == "12:55938053-55938510_12:55938053-55938637", c("UCB4.00983", "UCB4.00891", "UCB4.00932", "UCB1.00003")]
#       UCB4.00983 UCB4.00891 UCB4.00932 UCB1.00003
#        0.3937729          0          1        NaN
sashimi(junction = "12:55938053-55938914", cell = c("UCB4.00983", "UCB4.00891", "UCB4.00932", "UCB1.00003"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/")


as.numeric(which(rowSums(is.na(MXE_psi[,19:20]))==1))
 [1]   2  11  29  35 265 266 268 290 291 292 293 294 295 296 297 298 299 300 301
[20] 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320
[39] 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339
[58] 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358
[77] 359 364 368 378 379 387 389 391 393 394 566 568 572 602 603 604

as.numeric(which(rowSums(is.na(MXE_psi[,19:20]))==2))
 [1]   1   3  21  22  24  25  26  27  28  30  31  34  36  37  38 267 361 362 363
[20] 367 369 370 371 375 376 377 380 381 382 383 384 385 386 388 390 395 396 397
[39] 398 399 400 401 402 403 405 406 564 565 567 569 570 575 576 578 579 580 581
[58] 582 583 584 601 605 606

MXE_psi[265,which(!is.na(MXE_psi[265,]))]
sashimi(junction = "12:62391887-62393052", cell = c("UCB3.01558", "UCB3.01587", "UCB3.01648", "UCB1.00001"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/")
sashimi(junction = "12:62391887-62393052", cell = c("UCB3.01558", "UCB3.00684", "UCB3.01587", "UCB3.01648"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/")



MXE_psi[22,which(!is.na(MXE_psi[22,]))]
MXE_psi[30,which(!is.na(MXE_psi[30,]))]
MXE_psi[31,which(!is.na(MXE_psi[31,]))]
MXE_psi[34,which(!is.na(MXE_psi[34,]))]


sashimi(junction = "10:114792673-114801303", cell = c("UCB3.00523", "UCB3.01212"), outDir = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/Sashimi/")









#### SE


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
PSI_tu[is.na(PSI_tu)] <- 100
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 0, min.genes = 0, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = -1.5, x.high.cutoff = 2, y.cutoff = -3)
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
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/Seurat_100_PC10_tSNE.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()


ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
ucb <- RunTSNE(object = ucb, dims.use = 1:20, do.fast = TRUE)
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/Seurat_100_PC20_tSNE.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()



ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/SE_Seurat_NA_to_100.RData")



#### Rtsne


library(Rtsne)
library(ggplot2)
Rtsne_out <- Rtsne(t(PSI_tu[ucb@var.genes, ]), perplexity=500, theta = 1)

tsne <- as.data.frame(Rtsne_out$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/Rtsne.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



#### ICA
library(ica)
imod2 <- icaimax(t(PSI_tu[ucb@var.genes, ]),2)

se_ica <- as.data.frame(imod2$Y)
colnames(se_ica) <- c("IC_1", "IC_2")
se_ica$CellType <- as.factor(Cell_type$CellType)
se_ica$Individual <- as.factor(Cell_type$Individual)


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/ICA.pdf")
ggplot(se_ica, aes(x = IC_1, y = IC_2, colour = CellType))+
  geom_point()
ggplot(se_ica, aes(x = IC_1, y = IC_2, colour = Individual))+
  geom_point()
dev.off()



#### ICA + Rtsne
imod2 <- icaimax(t(PSI_tu[ucb@var.genes, ]),20)

ica_Rtsne_se <- Rtsne(imod2$Y, perplexity=500, theta = 1)
tsne <- as.data.frame(ica_Rtsne_se$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/ICA_Rtsne.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

ica_Rtsne_se2 <- Rtsne(imod2$Y, perplexity=1000, theta=1, pca=TRUE)
tsne <- as.data.frame(ica_Rtsne_se2$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/ICA_Rtsne_2.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()






#### A3SS =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")
PSI <- A3SS_psi[, 15:ncol(A3SS_psi)]
row.names(PSI) <- A3SS_psi$as2
rs <- rowSums(!is.na(PSI))
PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- 100
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

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
tsne <- tsne[row.names(Cell_type),]
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Seurat/tSNE_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Seurat/SE_Seurat.RData")




library(Rtsne)
Rtsne_out <- Rtsne(t(PSI_tu[ucb@var.genes, row.names(Cell_type)]), perplexity=500, theta = 1)

tsne <- as.data.frame(Rtsne_out$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/Seurat/Rtsne.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()














#### AFE ========================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")
PSI <- AFE_psi[, 12:ncol(AFE_psi)]
row.names(PSI) <- AFE_psi$as2
rs <- rowSums(!is.na(PSI))
PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- 100
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)


library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

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
tsne <- tsne[row.names(Cell_type),]
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Seurat/tSNE_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Seurat/AFE_Seurat.RData")





#### ICA


afe_ica <- icaimax(t(PSI_tu[ucb@var.genes, ]),20)

afe_icas <- as.data.frame(afe_ica$Y)
colnames(afe_icas) <- paste("IC", 1:20, sep="")
afe_icas <- afe_icas[row.names(Cell_type), ]
afe_icas$CellType <- as.factor(Cell_type$CellType)
afe_icas$Individual <- as.factor(Cell_type$Individual)


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Seurat/ICA_1.pdf")
ggplot(afe_icas, aes(x = IC1, y = IC2, colour = CellType))+
  geom_point()
ggplot(afe_icas, aes(x = IC1, y = IC2, colour = Individual))+
  geom_point()
dev.off()


#### Rtsne

library(Rtsne)
Rtsne_out <- Rtsne(t(PSI_tu[ucb@var.genes, row.names(Cell_type)]), perplexity=500, theta = 1)

tsne <- as.data.frame(Rtsne_out$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Seurat/Rtsne.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


#### ICA + Rtsne

ica_Rtsne <- Rtsne(afe_icas[,1:20], perplexity=500, theta = 1)
tsne <- as.data.frame(ica_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/Seurat/ICA_Rtsne.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()






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


PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- 0
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

library(Seurat)
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = PSI_tu, min.cells = 2, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 1, y.cutoff = -3)
length(x = ucb@var.genes)
ucb@var.genes -> hvgene0

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
identical(row.names(tsne), row.names(Cell_type))
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/tSNE_0.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
dev.off()












PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- -.1
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

psd <- apply(PSI_tu, 1, sd)
pmean <- rowMeans(PSI_tu)

sum(pmean > -.098 & psd > 0.05)
# [1] 10932

PSI_tu <- PSI[row.names(PSI_tu)[pmean > -.098 & psd > 0.05], ]
PSI_tu[is.na(PSI_tu)] <- 100

pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)

#barplot(pca_psi$eig[,2],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")

dim(pca_psi$svd$U)

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:20], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()








PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- -.1
rs1 <- rowSums(PSI_tu>0&PSI_tu<1)

pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:10], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_2_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()




















#################################################################################################################################################
######################################################## - 0.1 ##################################################################################
#################################################################################################################################################


PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- -0.1
pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)
pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:20], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:10], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:5], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(PSI_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca_psi$svd$U
row.names(pcs) <- colnames(PSI_tu)
colnames(pcs) <- paste("PC",1:100,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)


pca_Rtsne <- Rtsne(t(combat[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



pca_Rtsne <- Rtsne(t(combat1[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[2:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc2_10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[3:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc3_10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[2:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc2_20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[3:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-0.1_combat_pc3_20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



#################################################################################################################################################
######################################################## 100 ####################################################################################
#################################################################################################################################################






PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- 100
pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)
pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:20], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:10], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()

pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:5], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(PSI_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca_psi$svd$U
row.names(pcs) <- colnames(PSI_tu)
colnames(pcs) <- paste("PC",1:100,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)

pca_Rtsne <- Rtsne(t(combat[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_combat1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()









#################################################################################################################################################
######################################################## - 1 ####################################################################################
#################################################################################################################################################




library(ggplot2)
PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- -1
pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)
pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:20], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:10], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:5], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(PSI_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca_psi$svd$U
row.names(pcs) <- colnames(PSI_tu)
colnames(pcs) <- paste("PC",1:100,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)

pca_Rtsne <- Rtsne(t(combat[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_-1_combat1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()







#################################################################################################################################################
######################################################## 0 ######################################################################################
#################################################################################################################################################



library(Rtsne)
library(ggplot2)
PSI_tu <- PSI[rs>=10, ]
PSI_tu[is.na(PSI_tu)] <- 0
pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)
pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:20], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:10], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(pca_psi$svd$U[,1:5], perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(PSI_tu),]
#Define center as batch
#mod <- model.matrix(~as.factor(sample.info$Cell_type), data=exprs)
pcs <- pca_psi$svd$U
row.names(pcs) <- colnames(PSI_tu)
colnames(pcs) <- paste("PC",1:100,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = NULL, mean.only = T)

pca_Rtsne <- Rtsne(t(combat[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat1_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat1_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:5,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_0_combat1_pc5.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()




#### -0.1 and remove PC1
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

PSI_tu <- PSI[rs>=1, ]
PSI_tu[is.na(PSI_tu)] <- 100
pca_psi <- FactoMineR::PCA(t(PSI_tu),ncp = 100,graph = F)

#### Combat
library(sva)
sample_batch <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/batch.list", header = F, row.names = 1, stringsAsFactors = F)
sample_batch <- sample_batch[colnames(PSI_tu),]
#Define center as batch
mod <- model.matrix(~as.factor(Cell_type$CellType), data=PSI_tu)
pcs <- pca_psi$svd$U
row.names(pcs) <- colnames(PSI_tu)
colnames(pcs) <- paste("PC",1:100,seep="")

combat <- sva::ComBat(dat=t(pcs), batch=sample_batch, par.prior=TRUE, prior.plots=FALSE, mod = NULL)
combat1 <- sva::ComBat(dat=t(pcs), batch=sample_batch, prior.plots=FALSE, par.prior=TRUE, mod = mod)


library(Rtsne)
pca_Rtsne <- Rtsne(t(combat[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_min1_combat_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_min1_combat_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:10,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_min1_combatmod_pc10.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()


pca_Rtsne <- Rtsne(t(combat1[1:20,]), perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/Seurat/PCA_Rtsne_100_min1_combatmod_pc20.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()








#### 
#load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
load("/mnt/nfs_nas/Lab/People/tangchao/Project/UCB/ICA/All_Unanno_SJ/data/sj12.RData")
#load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
load("/mnt/nfs_nas/Lab/People/tangchao/Project/UCB/ICA/All_Unanno_SJ/data/psi_list_left_right_table.RData")
#Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
Cell_type <- read.table("/mnt/nfs_nas/Lab/People/tangchao/Project/UCB/ICA/All_Unanno_SJ/data/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)

colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

nsj <- sj[sj$annotation == 0, "sj"]
psi_tu <- psi_sj_same_end_table[nsj[nsj %in% row.names(psi_sj_same_end_table)], row.names(Cell_type)]
psi_tu <- psi_tu[rowSums(!is.na(psi_tu))>=10, ]

psi_tu[is.na(psi_tu)] <- 100



#### Rtsne



library(Rtsne)
library(ggplot2)
tsne_out_ucb <- Rtsne(t(psi_tu), perplexity=500, theta = 1)
plot(tsne_out_ucb$Y, xlab = "t-SNE1", ylab = "t-SNE2", main = "UCB SE psi",pch = 20,
     cex = .5)

tsne <- as.data.frame(tsne_out_ucb$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/nfs_nas/Lab/People/tangchao/Project/UCB/ICA/All_Unanno_SJ/Figure/Rtsne_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()



#### PCA + Rtsne



pca_psi <- FactoMineR::PCA(t(psi_tu),ncp = 20,graph = F)

barplot(pca_psi$eig[,2],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")

dim(pca_psi$svd$U)

pca_Rtsne <- Rtsne(pca_psi$svd$U, perplexity=1000, theta = 1)
tsne <- as.data.frame(pca_Rtsne$Y)
colnames(tsne) <- c("tSNE_1", "tSNE_2")
tsne$CellType <- as.factor(Cell_type$CellType)
tsne$Individual <- as.factor(Cell_type$Individual)

pdf("/mnt/nfs_nas/Lab/People/tangchao/Project/UCB/ICA/All_Unanno_SJ/Figure/PCA_Rtsne_1.pdf")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = Individual))+
  geom_point()
dev.off()









































