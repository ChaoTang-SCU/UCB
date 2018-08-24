#### intron-centric PSI Vs. Salmon
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"
psi_sj_same_start_table[c(SJ34,SJ35,SJ36,SJ37), ] -> CD45_PSI
row.names(CD45_PSI) <- c("SJ34","SJ35","SJ36","SJ37")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

load(file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
CD45_PSI <- CD45_PSI[, colnames(cov_tu)]
data.frame(colSums(cov_tu[1:4,])/colSums(cov_tu)) -> RA
colnames(RA) <- "Salmon"
RA$PSI <- CD45_PSI["SJ34",]
RA <- na.omit(RA)
cbind(RA, Cell_type[row.names(RA),]) -> RA

library(ggplot2)
summary(lm(RA$Salmon~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Salmon),2), "\n", 
			"p = ", signif(cor.test(RA$PSI, RA$Salmon)$p.value,2), "\n", 
			"R2 = ", "0.97", "\n",
      "n = ", nrow(RA), sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RA+.pdf")
ggplot(data = RA, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()


#### RB+
data.frame(colSums(cov_tu[5:6,])/colSums(cov_tu)) -> RB
colnames(RB) <- "Salmon"
RB$PSI <- CD45_PSI["SJ35",]
RB <- na.omit(RB)
cbind(RB, Cell_type[row.names(RB),]) -> RB

summary(lm(RB$Salmon~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Salmon),2), "\n", 
			"p = ", signif(cor.test(RB$PSI, RB$Salmon)$p.value,2), "\n", 
			"R2 = ", "0.98", "\n",
      "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RB+.pdf")
ggplot(data = RB, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()



#### RO
t(data.frame(cov_tu[8,]/colSums(cov_tu))) -> RO
colnames(RO) <- "Salmon"
RO <- as.data.frame(RO)
RO$PSI <- as.numeric(CD45_PSI["SJ37",])
RO <- na.omit(RO)
cbind(RO, Cell_type[row.names(RO),]) -> RO

summary(lm(RO$Salmon~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Salmon),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Salmon)$p.value,2), "\n", 
            "R2 = ", "0.97", "\n",
            "n = ", nrow(RO), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RO+.pdf")
ggplot(data = RO, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()






#### intron-centric PSI Vs. Stringtie



library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"t_data.ctab",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 11, sep = "\t")
    setnames(tmp, "cov",strsplit(x,"/")[[1]][10])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[cov_merge<10] <- 0

cov_merge[,row.names(tsne)] -> cov_tu


#### CD45: use exon expression as cutoff and splice junction for validate
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "median",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_exp/")
system.time(median_merge<-multmerge(path))
median_merge <- as.data.frame(median_merge)
colnames(median_merge) <- substr(colnames(median_merge),1,10)

colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB3", replacement="UCB4")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB1", replacement="UCB3")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB5", replacement="UCB1")

row.names(median_merge) <- paste("E", 3:7, sep = "")
median_merge <- t(median_merge)
median_merge <- data.frame(median_merge)
median_merge$Exon_Sum <- as.numeric(rowSums(median_merge))

median_merge$difE3E7 <- median_merge$E3/(median_merge$E3+median_merge$E7)
median_merge[is.na(median_merge)] <- 0
exon_exp_tu <- median_merge[log10(median_merge$Exon_Sum)>=1.3 & median_merge$difE3E7>=.4 & median_merge$difE3E7<=.75, ]


sum(colnames(cov_tu) %in% row.names(exon_exp_tu))
cov_tu[, !colnames(cov_tu) %in% row.names(exon_exp_tu)] <- 0

#### validate use splice junction

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"

as.data.frame(te[.(c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37))])[,-1] -> CD45_SJ
colnames(CD45_SJ) <- substr(colnames(CD45_SJ),1,10)
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

CD45_SJ <- CD45_SJ[,colnames(cov_tu)]

cov_tu["RA", apply(CD45_SJ[c("SJ34","SJ47"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAB", apply(CD45_SJ[c("SJ34","SJ45","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RABC", apply(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAC", apply(CD45_SJ[c("SJ34","SJ46","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RB", apply(CD45_SJ[c("SJ35","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RBC", apply(CD45_SJ[c("SJ35","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RC", apply(CD45_SJ[c("SJ36","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RO", apply(CD45_SJ[c("SJ37"),], 2, function(x) sum(is.na(x))!=0)] <- 0

save(cov_tu, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/stringtie_cov_tu.RData")



#### RA+


data.frame(colSums(cov_tu[1:4,])/colSums(cov_tu)) -> RA
colnames(RA) <- "Stringtie"
RA$PSI <- CD45_PSI["SJ34",]
RA <- na.omit(RA)
cbind(RA, Cell_type[row.names(RA),]) -> RA

library(ggplot2)
summary(lm(RA$Stringtie~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RA$PSI, RA$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.92",  "\n",
            "n = ", nrow(RA), sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RA+.pdf")
ggplot(data = RA, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()



#### RB+

data.frame(colSums(cov_tu[5:6,])/colSums(cov_tu)) -> RB
colnames(RB) <- "Stringtie"
RB$PSI <- CD45_PSI["SJ35",]
RB <- na.omit(RB)
cbind(RB, Cell_type[row.names(RB),]) -> RB

summary(lm(RB$Stringtie~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RB$PSI, RB$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.92", "\n",
            "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RB+.pdf")
ggplot(data = RB, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()




#### RO
t(data.frame(cov_tu[8,]/colSums(cov_tu))) -> RO
colnames(RO) <- "Stringtie"
RO <- as.data.frame(RO)
RO$PSI <- as.numeric(CD45_PSI["SJ37",])
RO <- na.omit(RO)
cbind(RO, Cell_type[row.names(RO),]) -> RO

summary(lm(RO$Stringtie~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.69",  "\n",
            "n = ", nrow(RO), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RO+.pdf")
ggplot(data = RO, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()


RO[rowSums(RO[1:2])>0&rowSums(RO[1:2])<2,] -> RO
summary(lm(RO$Stringtie~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.33", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RO+_2.pdf")
ggplot(data = RO, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()





#### intron-centric PSI Vs. RT-PCR


test <- read.xls("/mnt/data5/BGI/UCB/tangchao/CD45/CD45 transcripts order.xlsx", sheet = 2)


#### RA+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[1:4,3:5])*0.01))) -> RA
RA$PSI <- as.numeric(CD45_PSI["SJ34",row.names(RA)])
cbind(RA, Cell_type[row.names(RA),]) -> RA
RA <- na.omit(RA)

summary(lm(RA$RT.PCR~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$PSI, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.92", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RA+.pdf")
ggplot(data = RA, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$RT.PCR~RA$PSI))[1], slope = coef(lm(RA$RT.PCR~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$RT.PCR~RA$PSI))[1], slope = coef(lm(RA$RT.PCR~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()






#### RB+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[5:6,3:5])*0.01))) -> RB
RB$PSI <- as.numeric(CD45_PSI["SJ35",row.names(RB)])
cbind(RB, Cell_type[row.names(RB),]) -> RB
RB <- na.omit(RB)

summary(lm(RB$RT.PCR~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$PSI, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.90", sep = "")

grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RB+.pdf")
ggplot(data = RB, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$RT.PCR~RB$PSI))[1], slope = coef(lm(RB$RT.PCR~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$RT.PCR~RB$PSI))[1], slope = coef(lm(RB$RT.PCR~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()






#### R0+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[8,3:5])*0.01))) -> RO
RO$PSI <- as.numeric(CD45_PSI["SJ37",row.names(RO)])
cbind(RO, Cell_type[row.names(RO),]) -> RO
RO <- na.omit(RO)

summary(lm(RO$RT.PCR~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.66", sep = "")

grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RO+.pdf")
ggplot(data = RO, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$RT.PCR~RO$PSI))[1], slope = coef(lm(RO$RT.PCR~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$RT.PCR~RO$PSI))[1], slope = coef(lm(RO$RT.PCR~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()




#### intron-centric PSI Vs. Salmon
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"
psi_sj_same_start_table[c(SJ34,SJ35,SJ36,SJ37), ] -> CD45_PSI
row.names(CD45_PSI) <- c("SJ34","SJ35","SJ36","SJ37")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

load(file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/cov_tu.RData")
CD45_PSI <- CD45_PSI[, colnames(cov_tu)]
data.frame(colSums(cov_tu[1:4,])/colSums(cov_tu)) -> RA
colnames(RA) <- "Salmon"
RA$PSI <- CD45_PSI["SJ34",]
RA <- na.omit(RA)
cbind(RA, Cell_type[row.names(RA),]) -> RA

RA <- RA[rowSums(RA[1:2])>0,]
library(ggplot2)
summary(lm(RA$Salmon~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Salmon),2), "\n", 
      "p = ", signif(cor.test(RA$PSI, RA$Salmon)$p.value,2), "\n", 
      "R2 = ", "0.92", "\n", 
      "n = ", nrow(RA), sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RA+_no00.pdf")
ggplot(data = RA, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()



RA <- RA[rowSums(RA[1:2])>0 & rowSums(RA[1:2])<2,]
library(ggplot2)
summary(lm(RA$Salmon~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Salmon),2), "\n", 
      "p = ", signif(cor.test(RA$PSI, RA$Salmon)$p.value,2), "\n", 
      "R2 = ", "0.82", "\n", 
      "n = ", nrow(RA), sep = "")

grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RA+_no00_no11.pdf")
ggplot(data = RA, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$PSI))[1], slope = coef(lm(RA$Salmon~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()



#### RB+
data.frame(colSums(cov_tu[5:6,])/colSums(cov_tu)) -> RB
colnames(RB) <- "Salmon"
RB$PSI <- CD45_PSI["SJ35",]
RB <- na.omit(RB)
cbind(RB, Cell_type[row.names(RB),]) -> RB

RB <- RB[rowSums(RB[1:2])>0,]

summary(lm(RB$Salmon~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Salmon),2), "\n", 
      "p = ", signif(cor.test(RB$PSI, RB$Salmon)$p.value,2), "\n", 
      "R2 = ", "0.94", "\n", 
      "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RB+_no00.pdf")
ggplot(data = RB, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()


RB <- RB[rowSums(RB[1:2])>0 & rowSums(RB[1:2])<2,]

summary(lm(RB$Salmon~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Salmon),2), "\n", 
      "p = ", signif(cor.test(RB$PSI, RB$Salmon)$p.value,2), "\n", 
      "R2 = ", "0.87", "\n", 
      "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RB+_no00_no11.pdf")
ggplot(data = RB, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$PSI))[1], slope = coef(lm(RB$Salmon~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()


#### RO
t(data.frame(cov_tu[8,]/colSums(cov_tu))) -> RO
colnames(RO) <- "Salmon"
RO <- as.data.frame(RO)
RO$PSI <- as.numeric(CD45_PSI["SJ37",])
RO <- na.omit(RO)
cbind(RO, Cell_type[row.names(RO),]) -> RO

RO <- RO[rowSums(RO[1:2])>0,]

summary(lm(RO$Salmon~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Salmon),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Salmon)$p.value,2), "\n", 
            "R2 = ", "0.95", "\n",
            "n = ", nrow(RO), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RO+_no00.pdf")
ggplot(data = RO, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()



RO <- RO[rowSums(RO[1:2])>0 & rowSums(RO[1:2])<2,]

summary(lm(RO$Salmon~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Salmon),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Salmon)$p.value,2), "\n", 
            "R2 = ", "0.88", "\n",
            "n = ", nrow(RO), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Salmon_RO+_no00_no11.pdf")
ggplot(data = RO, aes(x = PSI, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$PSI))[1], slope = coef(lm(RO$Salmon~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()




#### intron-centric PSI Vs. Stringtie



library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(paste(filenames,"t_data.ctab",sep = "/"), function(x){
    tmp<-fread(x, header=TRUE, select = 11, sep = "\t")
    setnames(tmp, "cov",strsplit(x,"/")[[1]][10])
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/result")

system.time(cov_merge<-multmerge(path))

cov_merge <- as.data.frame(cov_merge)
row.names(cov_merge) <- c("RA","RAB","RABC","RAC","RB","RBC","RC","RO")

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(Combine@dr$tsne@cell.embeddings)

cov_merge[cov_merge<10] <- 0

cov_merge[,row.names(tsne)] -> cov_tu


#### CD45: use exon expression as cutoff and splice junction for validate
library(data.table)
multmerge = function(mypath){
  #need to change the pattern accordingly !!!
  filenames=list.files(path=mypath, pattern="txt", full.names=TRUE)
  datalist = lapply(filenames, function(x){
    tmp<-fread(x,header=TRUE, select = 4, sep = "\t")
    setnames(tmp, "median",tail(strsplit(x,"/")[[1]],n=1))
    return(tmp)})
  Reduce(function(x,y) {cbind(x,y)}, datalist)
}

path<-file.path("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v4/Exon37_exp/")
system.time(median_merge<-multmerge(path))
median_merge <- as.data.frame(median_merge)
colnames(median_merge) <- substr(colnames(median_merge),1,10)

colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB3", replacement="UCB4")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB1", replacement="UCB3")
colnames(median_merge) <- gsub(colnames(median_merge), pattern="UCB5", replacement="UCB1")

row.names(median_merge) <- paste("E", 3:7, sep = "")
median_merge <- t(median_merge)
median_merge <- data.frame(median_merge)
median_merge$Exon_Sum <- as.numeric(rowSums(median_merge))

median_merge$difE3E7 <- median_merge$E3/(median_merge$E3+median_merge$E7)
median_merge[is.na(median_merge)] <- 0
exon_exp_tu <- median_merge[log10(median_merge$Exon_Sum)>=1.3 & median_merge$difE3E7>=.4 & median_merge$difE3E7<=.75, ]


sum(colnames(cov_tu) %in% row.names(exon_exp_tu))
cov_tu[, !colnames(cov_tu) %in% row.names(exon_exp_tu)] <- 0

#### validate use splice junction

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te.RData")
library(data.table)
SJ34 = "1:198692374-198696711"
SJ35 = "1:198692374-198699563"
SJ45 = "1:198696910-198699563"
SJ46 = "1:198696910-198702386"
SJ56 = "1:198699705-198702386"
SJ57 = "1:198699705-198703297"
SJ67 = "1:198702531-198703297"
SJ47 = "1:198696910-198703297"
SJ36 = "1:198692374-198702386"
SJ37 = "1:198692374-198703297"

as.data.frame(te[.(c(SJ34,SJ35,SJ45,SJ46,SJ56,SJ57,SJ67,SJ47,SJ36,SJ37))])[,-1] -> CD45_SJ
colnames(CD45_SJ) <- substr(colnames(CD45_SJ),1,10)
row.names(CD45_SJ) <- c("SJ34","SJ35","SJ45","SJ46","SJ56","SJ57","SJ67","SJ47","SJ36","SJ37")

CD45_SJ <- CD45_SJ[,colnames(cov_tu)]

cov_tu["RA", apply(CD45_SJ[c("SJ34","SJ47"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAB", apply(CD45_SJ[c("SJ34","SJ45","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RABC", apply(CD45_SJ[c("SJ34","SJ45","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RAC", apply(CD45_SJ[c("SJ34","SJ46","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RB", apply(CD45_SJ[c("SJ35","SJ57"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RBC", apply(CD45_SJ[c("SJ35","SJ56","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RC", apply(CD45_SJ[c("SJ36","SJ67"),], 2, function(x) sum(is.na(x))!=0)] <- 0
cov_tu["RO", apply(CD45_SJ[c("SJ37"),], 2, function(x) sum(is.na(x))!=0)] <- 0

save(cov_tu, file = "/mnt/data5/BGI/UCB/tangchao/CD45/CD45_v5/stringtie_cov_tu.RData")



#### RA+


data.frame(colSums(cov_tu[1:4,])/colSums(cov_tu)) -> RA
colnames(RA) <- "Stringtie"
RA$PSI <- CD45_PSI["SJ34",]
RA <- na.omit(RA)
cbind(RA, Cell_type[row.names(RA),]) -> RA

RA <- RA[rowSums(RA[1:2])>0,]
library(ggplot2)
summary(lm(RA$Stringtie~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RA$PSI, RA$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.85", "\n", 
            "n = ", nrow(RA), sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RA+_no00.pdf")
ggplot(data = RA, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()


RA <- RA[rowSums(RA[1:2])>0 & rowSums(RA[1:2])<2,]
library(ggplot2)
summary(lm(RA$Stringtie~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RA$PSI, RA$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.77", "\n", 
            "n = ", nrow(RA), sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RA+_no00_no11.pdf")
ggplot(data = RA, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$PSI))[1], slope = coef(lm(RA$Stringtie~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()

 

#### RB+

data.frame(colSums(cov_tu[5:6,])/colSums(cov_tu)) -> RB
colnames(RB) <- "Stringtie"
RB$PSI <- CD45_PSI["SJ35",]
RB <- na.omit(RB)
cbind(RB, Cell_type[row.names(RB),]) -> RB

RB <- RB[rowSums(RB[1:2])>0,]

summary(lm(RB$Stringtie~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RB$PSI, RB$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.84", "\n", 
            "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RB+_no00.pdf")
ggplot(data = RB, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()



RB <- RB[rowSums(RB[1:2])>0 & rowSums(RB[1:2])<2,]

summary(lm(RB$Stringtie~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RB$PSI, RB$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.72", "\n", 
            "n = ", nrow(RB), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RB+_no00_no11.pdf")
ggplot(data = RB, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$PSI))[1], slope = coef(lm(RB$Stringtie~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()




#### RO
t(data.frame(cov_tu[8,]/colSums(cov_tu))) -> RO
colnames(RO) <- "Stringtie"
RO <- as.data.frame(RO)
RO$PSI <- as.numeric(CD45_PSI["SJ37",])
RO <- na.omit(RO)
cbind(RO, Cell_type[row.names(RO),]) -> RO

RO <- RO[rowSums(RO[1:2])>0,]

summary(lm(RO$Stringtie~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.66", "\n", 
            "n = ", nrow(RO), sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RO+_no00.pdf")
ggplot(data = RO, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()


RO[rowSums(RO[1:2])>0&rowSums(RO[1:2])<2,] -> RO
summary(lm(RO$Stringtie~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$Stringtie),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$Stringtie)$p.value,2), "\n", 
            "R2 = ", "0.32", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_Stringtie_RO+_no00_no11.pdf")
ggplot(data = RO, aes(x = PSI, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$PSI))[1], slope = coef(lm(RO$Stringtie~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()





#### intron-centric PSI Vs. RT-PCR


test <- read.xls("/mnt/data5/BGI/UCB/tangchao/CD45/CD45 transcripts order.xlsx", sheet = 2)


#### RA+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[1:4,3:5])*0.01))) -> RA
RA$PSI <- as.numeric(CD45_PSI["SJ34",row.names(RA)])
cbind(RA, Cell_type[row.names(RA),]) -> RA
RA <- na.omit(RA)

summary(lm(RA$RT.PCR~RA$PSI))
lab = paste("cor = ", round(cor(RA$PSI, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$PSI, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.92", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RA+.pdf")
ggplot(data = RA, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$RT.PCR~RA$PSI))[1], slope = coef(lm(RA$RT.PCR~RA$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$RT.PCR~RA$PSI))[1], slope = coef(lm(RA$RT.PCR~RA$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()






#### RB+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[5:6,3:5])*0.01))) -> RB
RB$PSI <- as.numeric(CD45_PSI["SJ35",row.names(RB)])
cbind(RB, Cell_type[row.names(RB),]) -> RB
RB <- na.omit(RB)

summary(lm(RB$RT.PCR~RB$PSI))
lab = paste("cor = ", round(cor(RB$PSI, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$PSI, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.90", sep = "")

grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RB+.pdf")
ggplot(data = RB, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$RT.PCR~RB$PSI))[1], slope = coef(lm(RB$RT.PCR~RB$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$RT.PCR~RB$PSI))[1], slope = coef(lm(RB$RT.PCR~RB$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()






#### R0+
as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[8,3:5])*0.01))) -> RO
RO$PSI <- as.numeric(CD45_PSI["SJ37",row.names(RO)])
cbind(RO, Cell_type[row.names(RO),]) -> RO
RO <- na.omit(RO)

summary(lm(RO$RT.PCR~RO$PSI))
lab = paste("cor = ", round(cor(RO$PSI, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$PSI, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.66", sep = "")

grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/PSI_Vs_RT.PCR_RO+.pdf")
ggplot(data = RO, aes(x = PSI, y = RT.PCR, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$RT.PCR~RO$PSI))[1], slope = coef(lm(RO$RT.PCR~RO$PSI))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = PSI, y = RT.PCR))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$RT.PCR~RO$PSI))[1], slope = coef(lm(RO$RT.PCR~RO$PSI))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()





#### Salmon Vs. RT-PCR



library(gdata)
test <- read.xls("/mnt/data5/BGI/UCB/tangchao/CD45/CD45 transcripts order.xlsx", sheet = 2)


#### RA+

as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[1:4,3:5])*0.01))) -> RA

cbind(RA, Cell_type[row.names(RA),]) -> RA
RA <- na.omit(RA)

summary(lm(RA$Salmon~RA$RT.PCR))
lab = paste("cor = ", round(cor(RA$Salmon, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$Salmon, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.87", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))


pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RA+.pdf")
ggplot(data = RA, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$RT.PCR))[1], slope = coef(lm(RA$Salmon~RA$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$RT.PCR))[1], slope = coef(lm(RA$Salmon~RA$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()


RA <- RA[rowSums(RA[,c("Salmon", "RT.PCR")])>0,]

summary(lm(RA$Salmon~RA$RT.PCR))
lab = paste("cor = ", round(cor(RA$Salmon, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$Salmon, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.80", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RA+_no00.pdf")
ggplot(data = RA, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$RT.PCR))[1], slope = coef(lm(RA$Salmon~RA$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Salmon~RA$RT.PCR))[1], slope = coef(lm(RA$Salmon~RA$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()




#### RB+


as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[5:6,3:5])*0.01))) -> RB

cbind(RB, Cell_type[row.names(RB),]) -> RB
RB <- na.omit(RB)

summary(lm(RB$Salmon~RB$RT.PCR))
lab = paste("cor = ", round(cor(RB$Salmon, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$Salmon, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.88", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RB+.pdf")
ggplot(data = RB, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$RT.PCR))[1], slope = coef(lm(RB$Salmon~RB$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$RT.PCR))[1], slope = coef(lm(RB$Salmon~RB$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()


RB <- RB[rowSums(RB[,c("Salmon", "RT.PCR")])>0,]

summary(lm(RB$Salmon~RB$RT.PCR))
lab = paste("cor = ", round(cor(RB$Salmon, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$Salmon, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.79", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RB+_no00.pdf")
ggplot(data = RB, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$RT.PCR))[1], slope = coef(lm(RB$Salmon~RB$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Salmon~RB$RT.PCR))[1], slope = coef(lm(RB$Salmon~RB$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()




#### RO


as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[8,3:5])*0.01))) -> RO

cbind(RO, Cell_type[row.names(RO),]) -> RO
RO <- na.omit(RO)

summary(lm(RO$Salmon~RO$RT.PCR))
lab = paste("cor = ", round(cor(RO$Salmon, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$Salmon, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.55", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RO+.pdf")
ggplot(data = RO, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$RT.PCR))[1], slope = coef(lm(RO$Salmon~RO$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$RT.PCR))[1], slope = coef(lm(RO$Salmon~RO$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()


RO <- RO[rowSums(RO[,c("Salmon", "RT.PCR")])>0,]

summary(lm(RO$Salmon~RO$RT.PCR))
lab = paste("cor = ", round(cor(RO$Salmon, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$Salmon, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.03", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Salmon_RO+_no00.pdf")
ggplot(data = RO, aes(x = RT.PCR, y = Salmon, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$RT.PCR))[1], slope = coef(lm(RO$Salmon~RO$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = RT.PCR, y = Salmon))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Salmon~RO$RT.PCR))[1], slope = coef(lm(RO$Salmon~RO$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()


#### Stringtie Vs. RT-PCR


#### RA+

as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[1:4,3:5])*0.01))) -> RA

cbind(RA, Cell_type[row.names(RA),]) -> RA
RA <- na.omit(RA)

summary(lm(RA$Stringtie~RA$RT.PCR))
lab = paste("cor = ", round(cor(RA$Stringtie, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$Stringtie, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.84", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RA+.pdf")
ggplot(data = RA, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$RT.PCR))[1], slope = coef(lm(RA$Stringtie~RA$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$RT.PCR))[1], slope = coef(lm(RA$Stringtie~RA$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()


RA <- RA[rowSums(RA[,c("Stringtie", "RT.PCR")])>0,]

summary(lm(RA$Stringtie~RA$RT.PCR))
lab = paste("cor = ", round(cor(RA$Stringtie, RA$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RA$Stringtie, RA$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.78", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RA+_no00.pdf")
ggplot(data = RA, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$RT.PCR))[1], slope = coef(lm(RA$Stringtie~RA$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")

ggplot(data = RA, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RA$Stringtie~RA$RT.PCR))[1], slope = coef(lm(RA$Stringtie~RA$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RA$CellType))+
  ggtitle("RA+")
dev.off()


#### RB+


as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[5:6,3:5])*0.01))) -> RB

cbind(RB, Cell_type[row.names(RB),]) -> RB
RB <- na.omit(RB)

summary(lm(RB$Stringtie~RB$RT.PCR))
lab = paste("cor = ", round(cor(RB$Stringtie, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$Stringtie, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.75", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RB+.pdf")
ggplot(data = RB, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$RT.PCR))[1], slope = coef(lm(RB$Stringtie~RB$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$RT.PCR))[1], slope = coef(lm(RB$Stringtie~RB$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()


RB <- RB[rowSums(RB[,c("Stringtie", "RT.PCR")])>0,]

summary(lm(RB$Stringtie~RB$RT.PCR))
lab = paste("cor = ", round(cor(RB$Stringtie, RB$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RB$Stringtie, RB$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.70", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RB+_no00.pdf")
ggplot(data = RB, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$RT.PCR))[1], slope = coef(lm(RB$Stringtie~RB$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")

ggplot(data = RB, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RB$Stringtie~RB$RT.PCR))[1], slope = coef(lm(RB$Stringtie~RB$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RB$CellType))+
  ggtitle("RB+")
dev.off()




#### RO


as.data.frame(do.call(rbind,lapply(split(test, test$Cell), function(x) colSums(x[8,3:5])*0.01))) -> RO

cbind(RO, Cell_type[row.names(RO),]) -> RO
RO <- na.omit(RO)

summary(lm(RO$Stringtie~RO$RT.PCR))
lab = paste("cor = ", round(cor(RO$Stringtie, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$Stringtie, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.08", sep = "")

library(grid)
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RO+.pdf")
ggplot(data = RO, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$RT.PCR))[1], slope = coef(lm(RO$Stringtie~RO$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$RT.PCR))[1], slope = coef(lm(RO$Stringtie~RO$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()


RO <- RO[rowSums(RO[,c("Stringtie", "RT.PCR")])>0,]

summary(lm(RO$Stringtie~RO$RT.PCR))
lab = paste("cor = ", round(cor(RO$Stringtie, RO$RT.PCR),2), "\n", 
            "p = ", signif(cor.test(RO$Stringtie, RO$RT.PCR)$p.value,2), "\n", 
            "R2 = ", "0.00", sep = "")
grob <- grobTree(textGrob(lab, x=0.05,  y=0.9, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))

pdf("/mnt/data5/BGI/UCB/tangchao/CD45/CD45_Stat3/RT.PCR_Vs_Stringtie_RO+_no00.pdf")
ggplot(data = RO, aes(x = RT.PCR, y = Stringtie, colour = CellType, shape = CellType))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$RT.PCR))[1], slope = coef(lm(RO$Stringtie~RO$RT.PCR))[2])+
  annotation_custom(grob = grob)+
  scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")

ggplot(data = RO, aes(x = RT.PCR, y = Stringtie))+
  geom_point()+
  geom_abline(intercept = coef(lm(RO$Stringtie~RO$RT.PCR))[1], slope = coef(lm(RO$Stringtie~RO$RT.PCR))[2], col = "red")+
  annotation_custom(grob = grob)+
  #scale_shape_manual(values=1:nlevels(RO$CellType))+
  ggtitle("RO")
dev.off()











