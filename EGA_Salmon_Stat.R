#### EGA Sample information
all <- list.files("/mnt/raid61/Science_2014_chen/EGA_metadata", recursive = T, pattern="Analysis_Sample_meta_info", full.names=T)

all_list <- list()
for(i in 1:length(all)){
	all_list[[i]] <- read.table(all[i],sep="\t", stringsAsFactors=F)
}
all_tab <- do.call(rbind, all_list)

test <- suppressWarnings(data.frame(do.call(rbind, strsplit(all_tab$V3, split = "[=;]")), stringsAsFactors=F))
valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}

gtf_infor <- data.frame(BIOMATERIAL_TYPE = as.vector(apply(test,1,function(x){valueFind(x,type = "BIOMATERIAL_TYPE")})),            
                        MOLECULE = as.vector(apply(test,1,function(x){valueFind(x,type = "MOLECULE")})),
                        DISEASE = as.vector(apply(test,1,function(x){valueFind(x,type = "DISEASE")})),            
                        BIOMATERIAL_PROVIDER = as.vector(apply(test,1,function(x){valueFind(x,type = "BIOMATERIAL_PROVIDER")})),            
                        CELL_TYPE = as.vector(apply(test,1,function(x){valueFind(x,type = "CELL_TYPE")})),
                        MARKERS = as.vector(apply(test,1,function(x){valueFind(x,type = "MARKERS")})),            
                        TISSUE_TYPE = as.vector(apply(test,1,function(x){valueFind(x,type = "TISSUE_TYPE")})),            
                        donor_id = as.vector(apply(test,1,function(x){valueFind(x,type = "donor_id")})),
                        DONOR_AGE = as.vector(apply(test,1,function(x){valueFind(x,type = "DONOR_AGE")})),            
                        DONOR_SEX = as.vector(apply(test,1,function(x){valueFind(x,type = "DONOR_SEX")})),
                        stringsAsFactors=F)

all_info <- cbind(all_tab[,1:2], gtf_infor)

#"activated_macrophage", "monocyte", "MK", "B_CELL", "CD4T", "CD8T", "CLP", "DC", "NK", "ECUVp", "ECUVr", "EB", "HSC", "macrophage_inf", "macrophage", "neutrophil"
all_info$CELL_TYPE_ABB <- as.character(factor(all_info$CELL_TYPE, labels=c("activated_macrophage", "monocyte", "MK", "B_CELL", "CD4T", "CD8T", "CLP", "DC", "NK", "ECUVp", "ECUVr", "EB", "HSC", "macrophage_inf", "macrophage", "neutrophil")))
colnames(all_info)[1:2] <- c("SAMEA", "ID")
write.table(all_info, "/mnt/raid61/Science_2014_chen/STAR_output/All_Sample_Information.txt", sep = "\t", quote = F, row.names = F)


#### EGA CD45 Salmon ----------------------------------------------------------------------------------------------------------------------------


all_salmon_quant <- list.files("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_bulk/EGA_CDD45_Salmon/transcripts_quant", recursive=T, pattern="quant.sf", full.names=T)

my_list <- list()
for (i in 1:length(all_salmon_quant)) {
	my_list[[i]] <- read.table(all_salmon_quant[[i]], sep = "\t", header = T)[5]
}
cd45_tab <- do.call(cbind, my_list)

colnames(cd45_tab) <- do.call(rbind, strsplit(all_salmon_quant, "/"))[,11]
row.names(cd45_tab) <- c("RA", "RAB", "RABC", "RAC", "RB", "RBC", "RC", "RO")
cd45_tab <- cd45_tab[, order(substr(colnames(cd45_tab),10,100))]
cd45_tab <- t(cd45_tab)


all_info <- read.table("/mnt/raid61/Science_2014_chen/STAR_output/All_Sample_Information.txt", sep = "\t", header = T, stringsAsFactors=F)
TISSUE_TYPE <- data.frame(TISSUE_TYPE=all_info[,"TISSUE_TYPE"], row.names=paste(all_info$ID, all_info$CELL_TYPE_ABB, sep = "_"), stringsAsFactors = F)
TISSUE_TYPE <- data.frame(TISSUE_TYPE = TISSUE_TYPE[row.names(cd45_tab),],row.names=row.names(cd45_tab))


pdf("~/EGA_CD45_Cov.pdf", width=8, height=16)
pheatmap(log10(cd45_tab+1), annotation_row = TISSUE_TYPE)
dev.off()



cd45_tab -> cd45_tab_p
cd45_tab_p[cd45_tab_p<50] <- 0
t(apply(cd45_tab_p,1,function(x)x/sum(x))) -> cd45_tab_p
cd45_tab_p[is.na(cd45_tab_p)] <- 0

pdf("~/EGA_CD45_isoform_ratio.pdf", width=8, height=16)
pheatmap(cd45_tab_p, annotation_row = TISSUE_TYPE)
dev.off()


Cord_Blood <- paste(all_info[all_info$TISSUE_TYPE == "cord blood","ID"],all_info[all_info$TISSUE_TYPE == "cord blood","CELL_TYPE_ABB"],sep="_")


pdf("~/EGA_Cord_Blood_CD45_Cov.pdf", width=8, height=16)
pheatmap(log10(cd45_tab[Cord_Blood,]+1))
dev.off()

pdf("~/EGA_Cord_Blood_CD45_isoform_ratio.pdf", width=8, height=16)
pheatmap(cd45_tab_p[Cord_Blood,])
dev.off()



#### PCA 

cd45_tab_p_sub <- cd45_tab_p[rowSums(cd45_tab_p)!=0,]

pca_mf <- FactoMineR::PCA(cd45_tab_p_sub,ncp = 100,graph = F)

barplot(pca_mf$eig[,2][1:10],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")

dim(pca_mf$svd$U)

pca_mf_result <- data.frame(pca_mf$svd$U, cell = substr(row.names(cd45_tab_p_sub), 10,50))

library(ggplot2)
pdf("~/PCA_of_EGA_CD45_isoform_ratio.pdf.pdf")
#ggplot(pca_mf_result,aes(x=X1,y=X2,col=cell))+
#  geom_point(size=2)+ #Size and alpha just for fun
#  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
#  theme_classic()+
#  #guides(colour=FALSE)+
#  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
#  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))+
#  geom_text(label=pca_mf_result$cell,colour="black",size=3)+
#  theme(legend.position = "")

ggplot(pca_mf_result,aes(x=X1,y=X2,col=cell,shape=cell))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))+
  scale_shape_manual(values=1:nlevels(pca_mf_result$cell))
  #geom_text(label=pca_mf_result$cell,colour="black",size=3)+
  #theme(legend.position = "")
dev.off()



test <- cd45_tab_p[Cord_Blood,]
cd45_tab_p_sub <- test[rowSums(test)!=0,]
pca_mf <- FactoMineR::PCA(cd45_tab_p_sub,ncp = 100,graph = F)
pca_mf_result <- data.frame(pca_mf$svd$U, cell = substr(row.names(cd45_tab_p_sub), 10,50))

pdf("~/PCA_of_EGA_Cord_Blood_CD45_isoform_ratio.pdf.pdf")
ggplot(pca_mf_result,aes(x=X1,y=X2,col=cell,shape=cell))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))+
  scale_shape_manual(values=1:nlevels(pca_mf_result$cell))
  #geom_text(label=pca_mf_result$cell,colour="black",size=3)+
  #theme(legend.position = "")
dev.off()




#### UCB ==============================================================================

load("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
cov_tu <- apply(cov_tu,2,function(x)x/sum(x))
cov_tu[is.na(cov_tu)] <- 0
cov_tu <- t(cov_tu)

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

Cell_type <- Cell_type[row.names(cov_tu),]

pdf("~/UCB_CD45_isoform_ratio.pdf", width=8, height=16)
pheatmap(cov_tu, annotation_row = Cell_type, show_rownames = F)
dev.off()

library(reshape2)
cell10 <- as.character(melt(mapply(function(x) {sample(row.names(Cell_type[Cell_type$CellType == x, ]),10)}, levels(Cell_type$CellType)))[,3])

pdf("~/UCB_CD45_isoform_ratio_subset.pdf", width=8, height=8)
pheatmap(cov_tu[cell10,], annotation_row = Cell_type[cell10,], show_rownames = F)
dev.off()


cov_tu[cell10,] -> test
Cell_type[cell10,] -> test2

row.names(test) <- as.character(paste(Cell_type[cell10,"CellType"]))

pdf("~/UCB_CD45_isoform_ratio_subset2.pdf", width=8, height=8)
pheatmap(test, show_rownames = T)
dev.off()


#### PCA
test <- cov_tu[rowSums(cov_tu)>0, ]
pca_mf <- FactoMineR::PCA(test,ncp = 100,graph = F)
barplot(pca_mf$eig[,2][1:10],xlab = "PC1 to PC10",ylab = "Percentage of variance(%)")
dim(pca_mf$svd$U)

pca_mf_result <- data.frame(pca_mf$svd$U, cell = as.character(Cell_type[row.names(test), "CellType"]))

pca_mf_result <- pca_mf_result[order(pca_mf_result$cell),]
pdf("~/PCA_of_UCB_CD45_isoform_ratio.pdf.pdf")
ggplot(pca_mf_result,aes(x=X1,y=X2,col=cell))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))
  #geom_text(label=pca_mf_result$cell,colour="black",size=3)+
  #theme(legend.position = "")
ggplot(pca_mf_result,aes(x=X1,y=X2,col=cell,shape=cell))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))+
  scale_shape_manual(values=1:nlevels(pca_mf_result$cell))
  #geom_text(label=pca_mf_result$cell,colour="black",size=3)+
  #theme(legend.position = "")
ggplot(pca_mf_result,aes(x=X1,y=X2,colour=cell))+
  geom_point(size=2)+ #Size and alpha just for fun
  #scale_color_manual(values = c(col.m,col.m.cross,col.n,col.n.cross,col.t,col.t.cross))+ #your colors here
  theme_classic()+
  #guides(colour=FALSE)+
  xlab(paste("PC1(",round(pca_mf$eig[,2][1],2),"%)",sep = ""))+
  ylab(paste("PC2(",round(pca_mf$eig[,2][2],2),"%)",sep = ""))+
  scale_colour_manual(values=c("grey", "grey", "grey", "red", "blue", "green"))
dev.off()


#### Density

test <- cov_tu[rowSums(cov_tu)>0, ]
test <- cbind(test, Cell_type[row.names(test), ])
pdf("~/Density_of_UCB_CD45_RABC_ratio.pdf")
ggplot(test, aes(RABC, colour = CellType))+
  geom_density()
dev.off()

tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RABC, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RABC")+
  ylab(label = "Percentage of Cells")


tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RABC, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

p1 <- ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RABC")+
  ylab(label = "Percentage of Cells")


tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RB, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

p2 <- ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RB")+
  ylab(label = "Percentage of Cells")


tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RAB, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

p3 <- ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RAB")+
  ylab(label = "Percentage of Cells")


tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RBC, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

p4 <- ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RBC")+
  ylab(label = "Percentage of Cells")


tab <- melt(mapply(function(x){table(cut(test[test$CellType == x,]$RO, seq(0,1,0.05), include.lowest = T))/nrow(test)},levels(test$CellType)))
colnames(tab) <- c("Interval", "CellType", "Percentage")
tab$Interval <- as.numeric(tab$Interval)/20

p5 <- ggplot(tab, aes(x = Interval, y = Percentage, colour = CellType))+
  geom_line()+
  xlab(label = "Percentage of Isoform RO")+
  ylab(label = "Percentage of Cells")


pdf("~/Percentage_of_UCB_CD45_isoform_ratio.pdf")
p1
p2
p3
p4
p5
dev.off()










